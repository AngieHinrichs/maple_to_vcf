#!/bin/bash
set -euo pipefail
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$TEST_DIR"

rm -rf output
mkdir -p output

# Basic test: no filtering options
suffix=plain
../maple_to_vcf.py \
    --input-directory ./input \
    --ref-fasta ./input/ref.fasta \
    --output output/test_$suffix.vcf \
    2> output/test_$suffix.stderr
diff output/test_$suffix.vcf expected/test_$suffix.vcf
diff output/test_$suffix.stderr expected/test_all.py.stderr

# Exclude sample but no other filters
suffix=exclude
../maple_to_vcf.py \
    --input-directory ./input \
    --ref-fasta ./input/ref.fasta \
    --exclude-file ./input/exclude.txt \
    --output output/test_$suffix.vcf \
    2> output/test_$suffix.stderr
diff output/test_$suffix.vcf expected/test_$suffix.vcf
diff output/test_$suffix.stderr expected/test_all.py.stderr

# Min minor allele frequency filter, no exclude
suffix=min_maf
../maple_to_vcf.py \
    --input-directory ./input \
    --ref-fasta ./input/ref.fasta \
    --min-maf 0.1 \
    --output output/test_$suffix.vcf \
    2> output/test_$suffix.stderr
diff output/test_$suffix.vcf expected/test_$suffix.vcf
diff output/test_$suffix.stderr expected/test_all.py.stderr

# Min minor allele frequency filter with exclude
suffix=min_maf_exclude
../maple_to_vcf.py \
    --input-directory ./input \
    --ref-fasta ./input/ref.fasta \
    --exclude-file ./input/exclude.txt \
    --min-maf 0.1 \
    --output output/test_$suffix.vcf \
    2> output/test_$suffix.stderr
diff output/test_$suffix.vcf expected/test_$suffix.vcf
diff output/test_$suffix.stderr expected/test_all.py.stderr

# Higher minor allele frequency filter, no exclude
suffix=high_min_maf
../maple_to_vcf.py \
    --input-directory ./input \
    --ref-fasta ./input/ref.fasta \
    --min-maf 0.35 \
    --output output/test_$suffix.vcf \
    2> output/test_$suffix.stderr
diff output/test_$suffix.vcf expected/test_$suffix.vcf
diff output/test_$suffix.stderr expected/test_all.py.stderr

# Higher minor allele frequency filter with exclude
suffix=high_min_maf_exclude
../maple_to_vcf.py \
    --input-directory ./input \
    --ref-fasta ./input/ref.fasta \
    --min-maf 0.35 \
    --exclude-file ./input/exclude.txt \
    --output output/test_$suffix.vcf \
    2> output/test_$suffix.stderr
diff output/test_$suffix.vcf expected/test_$suffix.vcf
diff output/test_$suffix.stderr expected/test_all.py.stderr

# Min non-N proportion filter, no exclude
suffix=min_non_N
../maple_to_vcf.py \
    --input-directory ./input \
    --ref-fasta ./input/ref.fasta \
    --min-non-N 0.75 \
    --output output/test_$suffix.vcf \
    2> output/test_$suffix.stderr
diff output/test_$suffix.vcf expected/test_$suffix.vcf
diff output/test_$suffix.stderr expected/test_all.py.stderr

# Min non-N proportion filter with exclude
suffix=min_non_N_exclude
../maple_to_vcf.py \
    --input-directory ./input \
    --ref-fasta ./input/ref.fasta \
    --exclude-file ./input/exclude.txt \
    --min-non-N 0.75 \
    --output output/test_$suffix.vcf \
    2> output/test_$suffix.stderr
diff output/test_$suffix.vcf expected/test_$suffix.vcf
diff output/test_$suffix.stderr expected/test_all.py.stderr

echo "All tests passed."

# Clean up
rm -rf output
