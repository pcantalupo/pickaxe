FROM virushunter/annotater

WORKDIR /

# Install Conda
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=/miniconda/condabin:$PATH
RUN apt-get update --fix-missing
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh
RUN bash /miniconda.sh -b -p miniconda

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name pickaxe > pickaxe.yml

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH=$PATH:/miniconda/envs/pickaxe/bin:/opt/pickaxe

# Install pickaxe from github
WORKDIR /opt
ADD https://api.github.com/repos/pcantalupo/pickaxe/git/refs/heads/master version.pickaxe.json
RUN git clone https://github.com/pcantalupo/pickaxe
ENV PERL5LIB "$PERL5LIB:/opt/pickaxe/lib"

RUN echo 'alias l="ls -l"' >> ~/.bashrc
RUN echo 'alias la="ls -la"' >> ~/.bashrc

