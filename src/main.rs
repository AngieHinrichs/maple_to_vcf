use anyhow::{Context, Result};
use bio::io::fasta;
use clap::Parser;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(author, version, about = "Convert MAPLE files to multi-sample VCF", long_about = None)]
struct Args {
    /// Directory containing input MAPLE files
    #[arg(short = 'i', long = "input-directory")]
    input_directory: PathBuf,

    /// File containing IDs to exclude from conversion
    #[arg(short = 'e', long = "exclude-file")]
    exclude_file: Option<PathBuf>,

    /// Minimum minor allele frequency to include a variant (default: 0.0)
    #[arg(short = 'm', long = "min-maf", default_value = "0.0")]
    min_maf: f64,

    /// Minimum proportion of non-N genotypes to include a variant (default: 0.5)
    #[arg(short = 'n', long = "min-non-N", default_value = "0.5")]
    min_non_n: f64,

    /// Reference FASTA file
    #[arg(short = 'r', long = "ref-fasta")]
    ref_fasta: PathBuf,

    /// Output multi-sample VCF file
    #[arg(short = 'o', long = "output")]
    output: PathBuf,
}

#[derive(Debug, Clone)]
struct Snv {
    position: usize,
    alt_base: char,
}

#[derive(Debug, Clone)]
struct VcfRecord {
    alt: Vec<char>,
    genotypes: Vec<String>,
}

impl VcfRecord {
    fn new(sample_count: usize) -> Self {
        VcfRecord {
            alt: Vec::new(),
            genotypes: vec!["0".to_string(); sample_count],
        }
    }
}

fn parse_maple_file(
    file_path: &Path,
    exclude_ids: &HashSet<String>,
) -> Result<Option<(String, Vec<Snv>)>> {
    let file = File::open(file_path)
        .with_context(|| format!("Failed to open file: {}", file_path.display()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Read first line to get sample ID
    let first_line = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("Empty file"))??;

    if !first_line.starts_with('>') {
        eprintln!(
            "File {} does not start with '>', skipping.",
            file_path.display()
        );
        return Ok(None);
    }

    let sample_id = first_line[1..].trim().to_string();

    if exclude_ids.contains(&sample_id) {
        return Ok(None);
    }

    let mut snvs = Vec::new();

    for line_result in lines {
        let line = line_result?;
        let words: Vec<&str> = line.trim().split('\t').collect();

        if words.is_empty() {
            continue;
        }

        match words[0] {
            "N" | "n" if words.len() == 3 => {
                let start: usize = words[1].parse()?;
                let length: usize = words[2].parse()?;
                for pos in start..start + length {
                    snvs.push(Snv {
                        position: pos,
                        alt_base: 'N',
                    });
                }
            }
            "-" if words.len() == 3 => {
                // Deletion, skip
                continue;
            }
            "A" | "C" | "G" | "T" | "a" | "c" | "g" | "t" if words.len() == 2 => {
                let position: usize = words[1].parse()?;
                let alt_base = words[0].chars().next().unwrap().to_ascii_uppercase();
                snvs.push(Snv {
                    position,
                    alt_base,
                });
            }
            _ => {
                eprintln!(
                    "Invalid MAPLE file format in {}, skipping.",
                    file_path.display()
                );
                return Ok(None);
            }
        }
    }

    Ok(Some((sample_id, snvs)))
}

fn parse_maple_files(
    maple_directory: &Path,
    exclude_ids: &HashSet<String>,
) -> Result<HashMap<String, Vec<Snv>>> {
    let mut maple_data = HashMap::new();

    for entry in fs::read_dir(maple_directory)? {
        let entry = entry?;
        let file_path = entry.path();

        if file_path.is_file() {
            match parse_maple_file(&file_path, exclude_ids) {
                Ok(Some((id, snvs))) => {
                    maple_data.insert(id, snvs);
                }
                Ok(None) => {
                    // File was skipped
                }
                Err(e) => {
                    eprintln!("Error parsing {}: {}", file_path.display(), e);
                }
            }
        }
    }

    Ok(maple_data)
}

fn build_vcf_records(
    maple_data: &HashMap<String, Vec<Snv>>,
    sample_ids: &[String],
) -> HashMap<usize, VcfRecord> {
    let sample_count = sample_ids.len();
    let mut vcf: HashMap<usize, VcfRecord> = HashMap::new();

    for sample_id in sample_ids {
        let snvs = &maple_data[sample_id];
        let ix = sample_ids.iter().position(|id| id == sample_id).unwrap();

        for snv in snvs {
            let record = vcf
                .entry(snv.position)
                .or_insert_with(|| VcfRecord::new(sample_count));

            if snv.alt_base == 'N' {
                // Don't add N as an alt base, just set genotype to "."
                record.genotypes[ix] = ".".to_string();
            } else {
                // Add alt base if not already present
                if !record.alt.contains(&snv.alt_base) {
                    record.alt.push(snv.alt_base);
                }
                // Set genotype to 1 + index of alt_base in alt list
                let alt_index = record.alt.iter().position(|&c| c == snv.alt_base).unwrap();
                record.genotypes[ix] = (alt_index + 1).to_string();
            }
        }
    }

    vcf
}

fn make_info_field(record: &VcfRecord, sample_count: usize) -> (String, f64, usize) {
    let n_count = record.genotypes.iter().filter(|&gt| gt == ".").count();
    let mut alt_counts = vec![0; record.alt.len()];

    for gt in &record.genotypes {
        if gt != "." && gt != "0" {
            if let Ok(idx) = gt.parse::<usize>() {
                if idx > 0 && idx <= alt_counts.len() {
                    alt_counts[idx - 1] += 1;
                }
            }
        }
    }

    let non_n_count = sample_count - n_count;
    let alt_freqs: Vec<f64> = if non_n_count > 0 {
        alt_counts
            .iter()
            .map(|&count| count as f64 / non_n_count as f64)
            .collect()
    } else {
        vec![0.0; alt_counts.len()]
    };

    let (maf, info_str) = if alt_counts.is_empty() {
        let info = format!("AC=0;AN={};AF=0", non_n_count);
        (0.0, info)
    } else {
        let minor_alt_freqs: Vec<f64> = alt_freqs
            .iter()
            .map(|&freq| if freq <= 0.5 { freq } else { 1.0 - freq })
            .collect();
        let maf = minor_alt_freqs
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);

        let ac_str = alt_counts
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let af_str = alt_freqs
            .iter()
            .map(|f| format!("{:.6}", f))
            .collect::<Vec<_>>()
            .join(",");
        let info = format!("AC={};AN={};AF={}", ac_str, non_n_count, af_str);
        (maf, info)
    };

    (info_str, maf, non_n_count)
}

fn write_vcf<W: Write>(
    writer: &mut W,
    ref_seq: &fasta::Record,
    vcf_data: &HashMap<usize, VcfRecord>,
    sample_ids: &[String],
    min_maf: f64,
    min_non_n: f64,
) -> Result<()> {
    let sample_count = sample_ids.len();
    if sample_count == 0 {
        eprintln!("No samples to write to VCF.");
        return Ok(());
    }

    // Write VCF header
    writeln!(writer, "##fileformat=VCFv4.2")?;
    writeln!(writer, "##contig=<ID={}>", ref_seq.id())?;
    writeln!(
        writer,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
    )?;

    let mut header = vec![
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
    ];
    header.extend(sample_ids.iter().map(|s| s.as_str()));
    writeln!(writer, "{}", header.join("\t"))?;

    // Write VCF records
    let mut positions: Vec<_> = vcf_data.keys().cloned().collect();
    positions.sort_unstable();

    let ref_len = ref_seq.seq().len();

    for position in positions {
        let record = &vcf_data[&position];

        if position - 1 >= ref_len {
            eprintln!(
                "Error: position {} exceeds reference length {}. \
                 That means the input MAPLE files are not consistent with the reference FASTA.",
                position, ref_len
            );
            std::process::exit(1);
        }

        let ref_base = ref_seq.seq()[position - 1] as char;
        let ref_base_upper = ref_base.to_ascii_uppercase();

        let alt_bases = if record.alt.is_empty() {
            "*".to_string()
        } else {
            record.alt.iter().fold(String::new(), |mut acc, &c| {
                if !acc.is_empty() {
                    acc.push(',');
                }
                acc.push(c);
                acc
            })
        };

        let (info_str, maf, non_n_count) = make_info_field(record, sample_count);
        let non_n_proportion = non_n_count as f64 / sample_count as f64;

        if maf < min_maf || non_n_proportion < min_non_n {
            continue;
        }

        let name = format!("{}{}{}", ref_base_upper, position, alt_bases);

        let mut vcf_line = vec![
            ref_seq.id().to_string(),
            position.to_string(),
            name,
            ref_base_upper.to_string(),
            alt_bases,
            ".".to_string(),
            "PASS".to_string(),
            info_str,
            "GT".to_string(),
        ];
        vcf_line.extend(record.genotypes.clone());

        writeln!(writer, "{}", vcf_line.join("\t"))?;
    }

    Ok(())
}

fn convert_maple_to_vcf(
    maple_directory: &Path,
    vcf_file: &Path,
    ref_seq: &fasta::Record,
    exclude_ids: &HashSet<String>,
    min_maf: f64,
    min_non_n: f64,
) -> Result<()> {
    let maple_data = parse_maple_files(maple_directory, exclude_ids)?;
    let mut sample_ids: Vec<String> = maple_data.keys().cloned().collect();
    sample_ids.sort();

    let vcf_records = build_vcf_records(&maple_data, &sample_ids);

    if vcf_file.extension().and_then(|s| s.to_str()) == Some("gz")
        || vcf_file
            .to_str()
            .map(|s| s.ends_with(".vcf.gz"))
            .unwrap_or(false)
    {
        let file = File::create(vcf_file)?;
        let mut encoder = GzEncoder::new(file, Compression::default());
        write_vcf(
            &mut encoder,
            ref_seq,
            &vcf_records,
            &sample_ids,
            min_maf,
            min_non_n,
        )?;
        encoder.finish()?;
    } else {
        let mut file = File::create(vcf_file)?;
        write_vcf(
            &mut file,
            ref_seq,
            &vcf_records,
            &sample_ids,
            min_maf,
            min_non_n,
        )?;
    }

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    if !args.input_directory.exists() {
        anyhow::bail!(
            "Input directory {} does not exist.",
            args.input_directory.display()
        );
    }

    let exclude_ids = if let Some(exclude_file) = &args.exclude_file {
        let file = File::open(exclude_file)?;
        let reader = BufReader::new(file);
        reader
            .lines()
            .map(|line| line.map(|s| s.trim().to_string()))
            .collect::<Result<HashSet<_>, _>>()?
    } else {
        HashSet::new()
    };

    let reader = fasta::Reader::from_file(&args.ref_fasta)?;
    let records: Vec<_> = reader.records().collect::<Result<_, _>>()?;

    if records.len() != 1 {
        anyhow::bail!("Reference FASTA must contain exactly one sequence.");
    }

    let ref_seq = &records[0];

    if args.min_maf < 0.0 || args.min_maf > 0.5 {
        anyhow::bail!("min-maf must be between 0.0 and 0.5");
    }

    if args.min_non_n < 0.0 || args.min_non_n > 1.0 {
        anyhow::bail!("min-non-N must be between 0.0 and 1.0");
    }

    convert_maple_to_vcf(
        &args.input_directory,
        &args.output,
        ref_seq,
        &exclude_ids,
        args.min_maf,
        args.min_non_n,
    )?;

    Ok(())
}
