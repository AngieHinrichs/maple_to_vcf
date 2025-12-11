# maple_to_vcf

maple_to_vcf converts a given directory containing MAPLE/"diff" files and a reference FASTA file into a multi-sample VCF with some simple filtering options.

## Installation

For now, go to https://github.com/AngieHinrichs/maple_to_vcf/actions/runs/20116785061 (I haven't done a proper release yet) and click to download the package for your OS and architecture.
For example, if you are on a Linux machine and `uname -m` says `x86_64`, then download maple_to_vcf-x86_64-unknown-linux-gnu.
That will download the file maple_to_vcf-x86_64-unknown-linux-gnu.zip to your machine; unzip that to get the maple_to_vcf executable.

If you would prefer to install through conda, docker or some other means, then please [file an issue](https://github.com/AngieHinrichs/maple_to_vcf/issues/new/choose) to let me know!

## Usage

Required inputs:
* an input directory that contains MAPLE files that all share the same reference; each MAPLE file must describe only one sequence.
* the fasta file for the same reference sequence used to generate all MAPLE files in the input directory.

Run `maple_to_vcf --help` to see descriptions of command line options:

```console
Usage: maple-to-vcf [OPTIONS] --input-directory <INPUT_DIRECTORY> --ref-fasta <REF_FASTA> --output <OUTPUT>

Options:
  -i, --input-directory <INPUT_DIRECTORY>
          Directory containing input MAPLE files
  -e, --exclude-file <EXCLUDE_FILE>
          File containing IDs to exclude from conversion
  -m, --min-maf <MIN_MAF>
          Minimum minor allele frequency to include a variant (default: 0.0) [default: 0.0]
  -n, --min-non-N <MIN_NON_N>
          Minimum proportion of non-N genotypes to include a variant (default: 0.5) [default: 0.5]
  -r, --ref-fasta <REF_FASTA>
          Reference FASTA file
  -o, --output <OUTPUT>
          Output multi-sample VCF file
  -h, --help
          Print help
  -V, --version
          Print version
```


## Python version: maple_to_vcf.py

This script was the initial implementation but it was too slow and consumed too much memory for real usage.  I'm mostly just leaving it here for the record and may remove it if I end up developing the rust implementation further and don't want to maintain the python too.

It requires the Biopython library.  Install like this:
```
pip install biopython
git clone https://github.com/AngieHinrichs/maple_to_vcf
maple_to_vcf/maple_to_vcf.py --help
```

Example usages can be found in tests/run_tests.sh, with example input files in tests/inputs/ .

Usage/options:
```
usage: maple_to_vcf.py [-h] -i INPUT_DIRECTORY [-e EXCLUDE_FILE] [-m MIN_MAF] [-n MIN_NON_N] -r REF_FASTA -o OUTPUT

Convert MAPLE files to multi-sample VCF.

options:
  -h, --help            show this help message and exit
  -i INPUT_DIRECTORY, --input-directory INPUT_DIRECTORY
                        Directory containing input MAPLE files
  -e EXCLUDE_FILE, --exclude-file EXCLUDE_FILE
                        File containing IDs to exclude from conversion
  -m MIN_MAF, --min-maf MIN_MAF
                        Minimum minor allele frequency to include a variant (default: 0.0)
  -n MIN_NON_N, --min-non-N MIN_NON_N
                        Minimum proportion of non-N genotypes to include a variant (default: 0.5)
  -r REF_FASTA, --ref-fasta REF_FASTA
                        Reference FASTA file
  -o OUTPUT, --output OUTPUT
                        Output multi-sample VCF file
```

