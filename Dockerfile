FROM virushunter/annotater

# Install Conda - https://fabiorosado.dev/blog/install-conda-in-docker/
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=/opt/conda/bin:$PATH

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict
RUN conda install samtools=1.14 bowtie2=2.4.5  cutadapt=4.1  prinseq=0.20.4  repeatmasker=4.1.2  megahit=1.2.9  diamond=2.0.15
RUN conda install perl-libwww-perl perl-xml-simple

# Install pickaxe from github
WORKDIR /opt
ADD https://api.github.com/repos/pcantalupo/pickaxe/git/refs/heads/master version.pickaxe.json
RUN git clone https://github.com/pcantalupo/pickaxe
ENV PATH=$PATH:/opt/pickaxe
ENV PERL5LIB "$PERL5LIB:/opt/pickaxe/lib"

RUN echo 'alias l="ls -l"' >> ~/.bashrc
RUN echo 'alias la="ls -la"' >> ~/.bashrc

