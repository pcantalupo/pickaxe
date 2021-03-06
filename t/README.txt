# add --shell to the following command to run without SLURM
pickaxe.pl --rm_default_index SRR11074364_
pickaxe.pl --collate both SRR11074364_ 
pickaxe.pl --extendcontigs SRR11074364_
pickaxe.pl --stats SRR11074364_


# If testing pickaxe in the Docker container, run the following commands
# This aligns reads against the BKV genome and searches contigs against the
# BKV BLAST database.  It obtains taxonomy informatoin for accession
# NC_001538 remotely from NCBI (precluding need to configure local taxonomy
# database)
pickaxe.pl --shell --assembler megahit --skip_repeatmasker --rm_default_index --virusindex /opt/bkv/bowtie2/NC_001538.1.fa --remotetax SRR11074364_
pickaxe.pl --collate both --virusfasta /opt/bkv/NC_001538.1.fa SRR11074364_
pickaxe.pl --extendcontigs SRR11074364_
pickaxe.pl --stats SRR11074364_
