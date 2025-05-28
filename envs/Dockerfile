FROM condaforge/mambaforge:4.10.1-0
LABEL org.opencontainers.image.source=https://github.com/NBISweden/metagenomics-depth-app \
    org.opencontainers.image.description="Dockerfile for the metagenomics depth app" \
    org.opencontainers.image.licenses=MIT \
    author="Lokeshwaran Manoharan" \
    email="lokeshwaran.manoharan@nbis.se"
ARG QUARTO_VERSION="1.7.31"
ENV TZ=Europe/Copenhagen
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set working directory
WORKDIR /project

# Install system dependencies
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb -O quarto.deb \
    && dpkg -i quarto.deb \
    && rm quarto.deb \
    && apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libgsl-dev \
    libglpk-dev \
    libhdf5-dev \
    libfftw3-dev \
    libopenblas-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy and install the Conda environment
COPY environment.yml ./
RUN conda config --set channel_priority strict \
    && mamba env update --name base --file environment.yml \
    && mamba clean --all --force-pkgs-dirs --yes

# Install R packages directly
RUN Rscript -e 'devtools::install_github("jeromefroe/circlepackeR")'
RUN pip install argparse
RUN /bin/bash -c "echo \"PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\W\[\033[00m\]\$ '\" >> /root/.bashrc "

RUN apt-get update && apt-get install less

# Start Bash shell by default
CMD /bin/bash
