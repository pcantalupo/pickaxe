FROM virushunter/annotater:v1.0.1

# Install Conda - https://fabiorosado.dev/blog/install-conda-in-docker/
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=/opt/conda/bin:$PATH

# Configure Conda channels to match environment.yml
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

COPY environment.yml /opt/environment.yml
RUN conda env create -f /opt/environment.yml
# Activate the environment by default by setting PATH directly
ENV PATH=/opt/conda/envs/pickaxe/bin:$PATH

# Install pickaxe from github
WORKDIR /opt
ADD https://api.github.com/repos/pcantalupo/pickaxe/git/refs/heads/master version.pickaxe.json
RUN git clone https://github.com/pcantalupo/pickaxe
ENV PATH=$PATH:/opt/pickaxe
ENV PERL5LIB="$PERL5LIB:/opt/pickaxe/lib"

RUN echo 'alias l="ls -l"' >> ~/.bashrc
RUN echo 'alias la="ls -la"' >> ~/.bashrc

