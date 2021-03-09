# pickaxe
Virus discovery pipeline

# Overview
Pickaxe identifies potential infectious agents in NGS sequencing reads. It reports on the number of alignments to each virus and metrics such as sequence coverage and average depth. Reads are assembled into contigs for novel virus discovery using search algorithms such as BLAST.

# Documentation
This README.md and for more details run `pickaxe.pl -h`

# Installation
+ Install [annotater](https://github.com/pcantalupo/annotater) and test that it is working.
+ Add the `pickaxe` folder to PATH and `pickaxe/lib` to PERL5LIB.
+ Install the following software. The listed version may not be absolutely required.  Other versions were not tested
  + samtools 1.2
  + bowtie2 2.3.3
  + cutadapt 1.18
  + PRINSEQ-lite 0.20.4
  + RepeatMasker open-4.0.7
  + clc_assembler 5.1.1.184548-201811011136 (or megahit 1.2.9)
  + blast+ 2.7.1+
  + rapsearch v2.24 (optional)
  + sra-toolkit 2.9.2 (optional)

# Test Dataset

This example uses 5000 reads from `SRR11074364` which is RNA-Seq of Human
RPTE cells infected by the virus BKV.  Switch to the `t` directory and run
the following.  Depending on your computer, this will take about 1 hour. 
The test commands run pickaxe with many default parameters and requires the
following 1) viral.1.1.genomic bowtie2 index and the viral.1.1.genomic.fna
file used to create the index, 2) viral.1.1.genomic blast database, and 3)
the NCBI taxonomy database.  Note: add --shell parameter to run without
submitting a SLURM sbatch job. If you do not have these items installed, see the Docker section below.

```
pickaxe.pl --rm_default_index SRR11074364_
pickaxe.pl --collate both SRR11074364_
pickaxe.pl --extendcontigs SRR11074364_
pickaxe.pl --stats SRR11074364_
```

This will generate the following output files:

```
SRR11074364_.pickaxe/
├── annotator/
├── assembly.ra.tsv
├── assembly.RM.fa
├── assembly.RM.fa.fai
├── extendedcontigs/
├── kv.vrs.hg19.bam
├── okfiles/
├── pickaxe.job
├── pickaxe.job.SRR11074364_.out
├── pickaxe.OK
├── report_BLAST.tsv
├── report_ViralRefSeq.tsv
├── software_versions.tsv
├── SRR11074364_1_5k.fastq -> /path/to/t/SRR11074364_1_5k.fastq
├── SRR11074364_2_5k.fastq -> /path/to/t/SRR11074364_2_5k.fastq
└── unmapped.fq.gz
```

The description of the files:
+ `unmapped.fq.gz` - 'non-host' reads that survived QC filtering and subtraction against host genomes.
+ `assembly.RM.fa` - raw contigs
+ `kv.vrs.hg19.bam` - bowtie2 alignments to the Viral RefSeq database
+ `annotator/` - folder containing raw results of [annotater](https://github.com/pcantalupo/annotater)
+ `report_ViralRefSeq.tsv` - summarization of the alignments to each virus reference sequence
+ `report_BLAST.tsv` - contig annotation results of the BLAST pipeline
found in the BAM file
+ `extendedcontigs` - folder containing extended contigs

## Expected results

The `report_ViralRefSeq.tsv` file should show the following: number of alignments to BKV ~1500 (column 'aligns'), sequence coverage ~87% (column 'seqcov') and average depth ~22 (column 'avgdepth').

The longest contig annotated by the BLAST pipeline as BKV should be ~1200bp (see `report_BLAST.tsv`).


# Docker

You can test pickaxe using the Docker image if you do not have any of the required databases or dependencies installed. First install Docker on your system and build the pickaxe image:

`docker build -t pickaxe .`

Then enter the `t` directory and following the instructions in the `README.txt`



