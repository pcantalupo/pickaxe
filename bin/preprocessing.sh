#!/bin/bash
# Preprocess reads with
# 1) cutadapt for adapter removal
# 2) prinseq to remove low quality sequences

# Careful - do not use your original fastq file with this script since
# script will delete it at the end.  Rather, use this script by running the
# preprocessing_parallel.sh wrapper script.

set -euo pipefail

fastqfile=$1
quality_cutoff=$2
entropy_cutoff=$3
length_cutoff=${4:-0}  # default 0
trim_left=${5:-0}  # default 0 
trim_right=${6:-0} # default 0
adapteroptions=${7:-} # default "" # adapteroptions is a string of cutadapt options (ex. -g GTTCAGAGTTCTACAGTCCGACGATC -a TCGTATGCCGTCTTCTGCTTG)

if [[ ! -f $fastqfile ]]; then
  echo "$fastqfile was not found"
  exit 1
fi

echo "Starting preprocessing"

base=$(basename "$fastqfile")
base_prefix=${base%.fastq}

#-n COUNT, --times=COUNT    Remove up to COUNT adapters from each read. Default: 1
#-O MINLENGTH, --overlap=MINLENGTH    Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
echo "Running cutadapt: trimming TruSeq + others"
cutadapt -g GTTTCCCAGTCACGATA    -a TATCGTGACTGGGAAAC \
         -g GACCATCTAGCGACCTCCAC -a GTGGAGGTCGCTAGATGGTC \
         -g GTTTCCCACTGGAGGATA   -a TATCCTCCAGTGGGAAAC \
         -g TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT                         -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
         -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -a CTGTCTCTTATACACATCT \
         $adapteroptions \
         -n 15 -O 5 -o "${fastqfile%.*}".cutadapt.fastq $fastqfile

# $quality_cutoff = 15  (for average quality cutoff and for trimming ends)
# $entropy_cutoff = 60
# $length_cutoff = 30
echo "Running prinseq-lite"
prinseq-lite.pl -fastq "$base_prefix".cutadapt.fastq \
                -lc_method entropy -lc_threshold "$entropy_cutoff" \
                -min_qual_mean "$quality_cutoff" -ns_max_p 5 \
                -trim_qual_right "$quality_cutoff" -trim_qual_left "$quality_cutoff" \
                -trim_left "$trim_left" -trim_right "$trim_right" \
                -min_len "$length_cutoff" -no_qual_header \
                -out_good "$base_prefix".cutadapt.prinseq -out_bad null
mv -f "$base_prefix".cutadapt.prinseq.fastq "$base_prefix".preprocessed.fastq

# Post cleanup
rm -vf $fastqfile                     # remove unneeded input fastq file
rm -vf "$base_prefix".cutadapt.fastq  # remove unneeded cutadapt.fastq file

echo "Done preprocessing"

exit 0
