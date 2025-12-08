#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from typing import TextIO
from collections import namedtuple
import gzip
import os
import sys

Record = namedtuple('Record', ['alt', 'genotypes'])


def parse_maple_file(file_path: str, exclude_ids: set):
    """Read in MAPLE file; if it looks like it's not valid MAPLE, or its ID is in exclude_ids,
    then return False, otherwise return True along with the ID and list of SNVs."""
    passes_filter = False
    sample_id = ""
    snvs = None
    with open(file_path, 'r') as mf:
        first_line = mf.readline()
        if first_line.startswith(">"):
            # ID is what follows ">":
            sample_id = first_line[1:].strip()
            if sample_id in exclude_ids:
                return False, sample_id, None
            snvs = []
            invalid = False
            for line in mf:
                words = line.strip().split("\t")
                if (words[0] == "N" or words[0] == "n") and len(words) == 3:
                    start = int(words[1])
                    length = int(words[2])
                    for pos in range(start, start + length):
                        snvs.append((pos, 'N'))
                elif words[0] == '-' and len(words) == 3:
                    # Deletion, skip
                    continue
                elif words[0] in ["A", "C", "G", "T", "a", "c", "g", "t"] and len(words) == 2:
                    position = int(words[1])
                    alt_base = words[0].upper()
                    snvs.append((position, alt_base))
                else:
                    # Not a valid maple file, skip this file
                    print(f"Invalid MAPLE file format in {file_path}, skipping.", file=sys.stderr)
                    invalid = True
                    break
            passes_filter = not invalid
        else:
            print(f"File {file_path} does not start with '>', skipping.", file=sys.stderr)
    return passes_filter, sample_id, snvs


def parse_maple_files(maple_directory: str, exclude_ids: set):
    """Attempt to parse all MAPLE files in the given directory, returning a dictionary
    mapping sample IDs to lists of SNVs."""
    maple_data = {}
    for maple_file in os.listdir(maple_directory):
        file_path = os.path.join(maple_directory, maple_file)
        passes, id, snvs = parse_maple_file(file_path, exclude_ids)
        if passes:
            maple_data[id] = snvs
    return maple_data


def build_vcf_records(maple_data: dict, sample_ids: list):
    """Build up a dict mapping positions to VCF alt bases and genotypes."""
    sample_count = len(sample_ids)
    vcf = {}
    for sample_id, snvs in maple_data.items():
        for position, alt_base in snvs:
            # For each SNV (or N), update the VCF alt and genotypes at that position
            ix = sample_ids.index(sample_id)
            if position not in vcf:
                vcf[position] = Record(alt=[], genotypes=["0"] * sample_count)
            if alt_base == 'N':
                # Don't add N as an alt base, just set genotype to "."
                vcf[position].genotypes[ix] = "."
            else:
                # Add alt base if not already present and set genotype to 1 + index of alt_base in alt list
                if alt_base not in vcf[position].alt:
                    vcf[position].alt.append(alt_base)
                vcf[position].genotypes[ix] = str(1 + vcf[position].alt.index(alt_base))
    return vcf


def make_info_field(record: Record, sample_count: int):
    """Calculate allele frequencies and build INFO field.
    Return INFO string, MAF, and non-N count for filtering."""
    N_count = record.genotypes.count(".")
    alt_counts = [0] * len(record.alt)
    for gt in record.genotypes:
        if gt != "." and gt != "0":
            alt_counts[int(gt) - 1] += 1
    non_N_count = sample_count - N_count
    alt_freqs = [0] * len(record.alt)
    for i in range(len(alt_counts)):
        alt_freqs[i] = alt_counts[i] / non_N_count if non_N_count > 0 else 0
    if len(alt_counts) == 0:
        maf = 0.0
        info_fields = [f"AC=0", f"AN={non_N_count}", f"AF=0" ]
    else:
        minor_alt_freqs = [x if x <= 0.5 else 1 - x for x in alt_freqs]
        maf = max(minor_alt_freqs)
        info_fields = [
            f"AC={','.join(map(str, alt_counts))}",
            f"AN={non_N_count}",
            f"AF={','.join(map(lambda x: f'{x:.6f}', alt_freqs))}"
        ]
    info_str = ";".join(info_fields)
    return info_str, maf, non_N_count


def write_vcf(vcf_file: TextIO, ref_seq: SeqIO.SeqRecord, vcf_data: dict, sample_ids: list, min_maf: float,
              min_non_N: float):
    """Write a VCF header and records to the given file, filtering by MAF and non-N proportion."""
    # Write VCF header
    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write(f"##contig=<ID={ref_seq.id}>\n")
    vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcf_file.write("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">\n")
    vcf_file.write("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n")
    vcf_file.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + sample_ids
    vcf_file.write("\t".join(header) + "\n")
    sample_count = len(sample_ids)
    if sample_count == 0:
        print("No samples to write to VCF.", file=sys.stderr)
        return
    # Write VCF records
    for position in sorted(vcf_data.keys()):
        record = vcf_data[position]
        if position - 1 >= len(ref_seq.seq):
            print(f"Error: position {position} exceeds reference length {len(ref_seq.seq)}. " +
                  "That means the input MAPLE files are not consistent with the reference FASTA.",
                  file=sys.stderr)
            sys.exit(1)
        ref_base = ref_seq.seq[position - 1].upper()
        alt_bases = ",".join(record.alt)
        if alt_bases == "":
            alt_bases = "*"
        info_str, maf, non_N_count = make_info_field(record, sample_count)
        non_N_proportion = non_N_count / sample_count
        if maf < min_maf or non_N_proportion < min_non_N:
            continue
        name = ref_base + str(position) + alt_bases
        vcf_line = [ref_seq.id, str(position), name, ref_base, alt_bases, ".", "PASS", info_str, "GT"] + record.genotypes
        vcf_file.write("\t".join(vcf_line) + "\n")


def convert_maple_to_vcf(maple_directory: str, vcf_file: str, ref_seq: SeqIO.SeqRecord, exclude_ids: set,
                         min_maf: float, min_non_N: float):
    # Build up a representation of a multi-sample VCF
    maple_data = parse_maple_files(maple_directory, exclude_ids)
    sample_ids = sorted(maple_data.keys())
    # Add SNVs to internal representation of multi-sample VCF
    vcf_records = build_vcf_records(maple_data, sample_ids)
    with gzip.open(vcf_file, 'wt') if vcf_file.endswith('.vcf.gz') else open(vcf_file, 'w') as vf:
        write_vcf(vf, ref_seq, vcf_records, sample_ids, min_maf, min_non_N)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert MAPLE files to multi-sample VCF.")
    parser.add_argument("-i", "--input-directory", required=True,
                        help="Directory containing input MAPLE files")
    parser.add_argument("-e", "--exclude-file",
                        help="File containing IDs to exclude from conversion")
    parser.add_argument("-m", "--min-maf", type=float, default=0.0,
                        help="Minimum minor allele frequency to include a variant (default: 0.0)")
    parser.add_argument("-n", "--min-non-N", type=float, default=0.5,
                        help="Minimum proportion of non-N genotypes to include a variant (default: 0.5)")
    parser.add_argument("-r", "--ref-fasta", required=True,
                        help="Reference FASTA file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output multi-sample VCF file")
    args = parser.parse_args()

    if not os.path.exists(args.input_directory):
        print(f"Input directory {args.input_directory} does not exist.")
        exit(1)

    if args.exclude_file:
        with open(args.exclude_file, 'r') as ef:
            exclude_ids = set(line.strip() for line in ef)
    else:
        exclude_ids = set()

    ref_seq_dict = SeqIO.to_dict(SeqIO.parse(args.ref_fasta, "fasta"))
    if len(ref_seq_dict) != 1:
        print("Reference FASTA must contain exactly one sequence.")
        exit(1)
    ref_seq = next(iter(ref_seq_dict.values()))

    if args.min_maf < 0.0 or args.min_maf > 0.5:
        print("min-maf must be between 0.0 and 0.5", file=sys.stderr)
        exit(1)

    if args.min_non_N < 0.0 or args.min_non_N > 1.0:
        print("min-non-N must be between 0.0 and 1.0", file=sys.stderr)
        exit(1)

    convert_maple_to_vcf(args.input_directory, args.output, ref_seq, exclude_ids,
                         args.min_maf, args.min_non_N)
