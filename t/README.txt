# Testing pickaxe using default parameters with the commands below requires the following
# 1) ref_viruses_rep_genomes bowtie2 index and the ref_viruses_rep_genomes.fa file used to create the index. These files need to be in the folder specified by the environment variable 'BOWTIE2_INDEXES'
# 2) ref_viruses_rep_genomes blast database
# 3) the NCBI taxonomy database
# Note: add --shell parameter to run without submitting a SLURM sbatch job
pickaxe.pl --rm_default_index SRR11074364_
pickaxe.pl --collate both SRR11074364_
pickaxe.pl --extendcontigs SRR11074364_
pickaxe.pl --stats SRR11074364_

# If testing pickaxe in the Docker container, run the following commands
# This aligns reads against the BKV genome and searches contigs against the
# BKV genome BLAST database.  It obtains taxonomy information for accession
# NC_001538 remotely from NCBI (precluding need to configure local taxonomy
# database). It skips the repeatmasker step since RM is difficult to
# install; one reason is that it must be configured interactively
docker run -v $(pwd):$(pwd) -w $(pwd) pickaxe pickaxe.pl --annotconf annotater.bkv.config --shell --assembler megahit --skip_repeatmasker --rm_default_index --virusindex /opt/bkv/bowtie2/NC_001538.1.fa --remotetax SRR11074364_
docker run -v $(pwd):$(pwd) -w $(pwd) pickaxe pickaxe.pl --annotconf annotater.bkv.config --collate both --virusfasta /opt/bkv/NC_001538.1.fa SRR11074364_
docker run -v $(pwd):$(pwd) -w $(pwd) pickaxe pickaxe.pl --extendcontigs SRR11074364_
docker run -v $(pwd):$(pwd) -w $(pwd) pickaxe pickaxe.pl --stats SRR11074364_

