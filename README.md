# maple_to_vcf.py

This script converts a given directory containing MAPLE/"diff" files and a reference FASTA file into a multi-sample VCF with some simple filtering options.

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

