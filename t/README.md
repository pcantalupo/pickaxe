# Testing Pickaxe

Testing pickaxe using default parameters with the commands below requires the following
- 1) ref_viruses_rep_genomes bowtie2 index and the ref_viruses_rep_genomes.fa file used to create the index. These files need to be in the folder specified by the environment variable 'BOWTIE2_INDEXES'
- 2) ref_viruses_rep_genomes blast database
- 3) the NCBI taxonomy database built by the `taxonomizr` R package
Note: add --shell parameter to run without submitting a SLURM sbatch job

```
pickaxe.pl --rm_default_index SRR11074364_
pickaxe.pl --collate both SRR11074364_
pickaxe.pl --extendcontigs SRR11074364_
pickaxe.pl --stats SRR11074364_
```

# Running Pickaxe with Docker

Build the docker image by running `docker build -t pickaxe .` in the top-level directory of this repository. Then in this directory, run the following commands to download the NCBI refseq virus Blast database and to build a small Bowtie2 index of 3 virus sequences:

```
mkdir -p refs && cd refs
wget https://ftp.ncbi.nlm.nih.gov/blast/db/ref_viruses_rep_genomes.tar.gz
tar xvzf ref_viruses_rep_genomes.tar.gz 
blastdbcmd -entry NC_001538.1,NC_001669.1,NC_001699.1 -db ref_viruses_rep_genomes > ref_viruses_rep_genomes.fa
bowtie2-build -f ref_viruses_rep_genomes.fa ref_viruses_rep_genomes
export BLASTDB=$(pwd -P)
export BOWTIE2_INDEXES=$(pwd -P)
cd ..
```

Afterwards, run the docker command below to align the reads against the 3 virus sequences, assemble contigs, and Blast the contigs against the NCBI viral refseq Blast database. If you do not have a local NCBI taxonomy database, add the `--remotetax` parameter which will obtain taxonomy info from NCBI.

```
docker run -w $(pwd) \
-v $(pwd):$(pwd) \
-v /PATH/TO/TAXONOMYDIR:/PATH/TO/TAXONOMYDIR \
-e TAXASQL=/PATH/TO/TAXONOMYDIR/accessionTaxa.sql \
-e NAMESDMP=/PATH/TO/TAXONOMYDIR/names.dmp \
-e NODESDMP=/PATH/TO/TAXONOMYDIR/nodes.dmp \
-e BLASTDB="$BLASTDB" \
-e BOWTIE2_INDEXES="$BOWTIE2_INDEXES" \
pickaxe pickaxe.pl --shell --assembler megahit --rm_default_index SRR11074364_
```

Then run the following docker command to generate the BLAST and ViralRefSeq reports:

`docker run -w $(pwd) -v $(pwd):$(pwd) -e BOWTIE2_INDEXES="$BOWTIE2_INDEXES" pickaxe pickaxe.pl --collate both SRR11074364_`

And finally run these two commands to extend the contigs and obtain stats about the analysis

```
docker run -v $(pwd):$(pwd) -w $(pwd) pickaxe pickaxe.pl --extendcontigs SRR11074364_
docker run -v $(pwd):$(pwd) -w $(pwd) pickaxe pickaxe.pl --stats SRR11074364_
```

