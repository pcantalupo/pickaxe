#!/bin/bash
# Requires wc, cut, split

set -euo pipefail

fastqfile=$1
adapterset=$2        # not implemented yet in "preprocessing.sh"
quality_cutoff=$3
entropy_cutoff=$4
length_cutoff=${5:-0} # default 0
threads=$6
trim_left=${7:-0}     # default 0
trim_right=${8:-0}    # default 0
adapteroptions=${9:-} # default "" # adapteroptions is a string of cutadapt options (ex. -g GTTCAGAGTTCTACAGTCCGACGATC -a TCGTATGCCGTCTTCTGCTTG)

if [ ! -f $fastqfile ]; then
  echo "$fastqfile was not found"
  exit 1
fi

totalL=$(wc -l $fastqfile | cut -f 1 -d ' ')
fastqreads=$((totalL/4))
echo -e "$fastqfile contains $fastqreads FASTQ reads"

# Split fastq file into smaller chunks to take advantage of parallel processing
#
tmp=$((totalL/(threads*4)))
lines=$((tmp*4))
split -l $lines -a4 $fastqfile

for x in xa*; do    # xaaaa, xaaab
  mv $x $x.fastq    # xaaaa.fastq (cutadapt requires suffix on filename)
  preprocessing.sh $x.fastq $quality_cutoff $entropy_cutoff $length_cutoff $trim_left $trim_right "$adapteroptions" \
                   > $x.preprocess.log 2>&1 &  # xaaaa.preprocess.log
done
wait


# Collate log files and preprocessed fastq files
#
base_fastq=$(basename "$fastqfile")   # temp.fastq
base_prefix_fastq=${base_fastq%.fastq}      # temp

# remove old files so they are not appended to
rm -vf $base_prefix_fastq.preprocess.log
rm -vf $base_prefix_fastq.preprocessed.fastq

cat xa*fastq > $base_prefix_fastq.preprocessed.fastq
rm -f xa*fastq
cat xa*preprocess.log
rm -f xa*preprocess.log

exit 0

