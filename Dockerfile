FROM virushunter/annotater

# Install samtools
RUN apt-get update && apt-get install --yes gcc libz-dev libncurses5-dev libncursesw5-dev
WORKDIR /
RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
RUN tar xjvf samtools-1.2.tar.bz2
WORKDIR /samtools-1.2
RUN make
WORKDIR /

# Install Text::Soundex for RepeatMasker
RUN cpanm Text::Soundex

# Install Conda
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=/miniconda/condabin:$PATH
RUN apt-get update --fix-missing
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh
RUN bash /miniconda.sh -b -p miniconda

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH=$PATH:/miniconda/envs/pickaxe-1.0/bin:/samtools-1.2

# Dump the details of the installed packages to a file for posterity
#RUN conda env export --name nf-core-nanoseq-1.1.0 > nf-core-nanoseq-1.1.0.yml
