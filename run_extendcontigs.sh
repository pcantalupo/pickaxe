#!/bin/bash
set -euo pipefail

ID=$1
BLASTTSV=$2
ASSEMBLYFA=$3

CONTIGS=contigs.txt

grep ^"$ID"$'\t' "$BLASTTSV" | grep $'\t'virus$'\t' | cut -f 2 > $CONTIGS

# Getting viral contigs
echo "Getting viral contigs"
while read c; do echo $c; getfastaseqs.pl -f -s "$ASSEMBLYFA" $c > $c.fa & done < $CONTIGS
wait
echo

# 1st extension
echo "1st extension"
while read c; do
  echo $c;
  extendsequence.pl --target "$ASSEMBLYFA" --query $c.fa --blastout $c.ext1.blast.tsv > $c.ext1.fa;
done < $CONTIGS
echo

# 2nd extension
echo "2nd extension" 
while read c; do
  file="$c.ext1.fa"
  if [[ -s $file ]]; then
    echo $c
    extendsequence.pl --target "$ASSEMBLYFA" --query $file --blastout $c.ext2.blast.tsv > $c.ext2.fa;
  else
    echo "$c is empty...skipping"
  fi
done < $CONTIGS
echo

# Check seq lengths
echo "Getting sequence lengths"
for ext in *ext*.fa; do
  if [[ -s $ext ]]; then
    echo $ext;
    seqlength.pl $ext;
  fi
done > seqlength.tsv
echo



