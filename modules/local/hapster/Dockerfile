FROM --platform=linux/amd64 condaforge/mambaforge:23.3.1-1

LABEL maintainer="orgeraj"
LABEL description="Hapster (formerly Polytect) - HLA somatic mutation calling pipeline"
LABEL version="1.0.1"

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    make \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    wget \
    curl \
    git \
    unzip \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment with Python and bioinformatics tools
# Using mamba for fast dependency resolution
RUN mamba create -n hapster python=3.7 -y && \
    mamba install -n hapster -c conda-forge -c bioconda \
    openjdk=8 \
    snakemake-minimal \
    samtools \
    bcftools \
    picard \
    minimap2 \
    gffread \
    pysam \
    biopython \
    tabix \
    htslib \
    -y && \
    mamba clean -afy

# Install R and required packages
RUN mamba install -n hapster -c conda-forge -c bioconda \
    r-base=4.2 \
    r-foreach \
    r-readr \
    r-dplyr \
    r-tibble \
    r-stringr \
    bioconductor-rtracklayer \
    bioconductor-biostrings \
    bioconductor-genomicfeatures \
    bioconductor-bsgenome \
    -y && \
    mamba clean -afy

# Activate the conda environment by default
ENV PATH=/opt/conda/envs/hapster/bin:$PATH
SHELL ["conda", "run", "-n", "hapster", "/bin/bash", "-c"]

# Clone Hapster repository and make CLI available
RUN git clone https://github.com/mctp/hapster.git /opt/hapster && \
    cd /opt/hapster && \
    if [ -f setup.py ] || [ -f pyproject.toml ]; then \
        pip install -e . ; \
    fi
ENV PATH=/opt/hapster:/opt/hapster/bin:$PATH

# Download tools tarball (bwa.kit 0.7.15, GATK 4.1.2.0)
# hapster setup_tools expects the tarball in resources/ and handles extraction
RUN mkdir -p /opt/hapster/resources && \
    wget -q https://storage.googleapis.com/mctp-open-share/hapster/tools-1.0.1.tar.gz \
        -O /opt/hapster/resources/tools-1.0.1.tar.gz

# Extract and build tools manually (setup_tools tries conda env create at the end
# which fails because our env already exists, so we run the steps ourselves)
RUN cd /opt/hapster && \
    tar -xf resources/tools-1.0.1.tar.gz -C build --skip-old-files && \
    bash -c '\
    set -e && \
    export HAPSTER_DIR=/opt/hapster && \
    export NCORES=$(nproc --all) && \
    cd $HAPSTER_DIR/build/tools-1.0.1 && \
    tar xf biobambam2-2.0.146-release-20191030105216-x86_64-linux-gnu.tar.xz && \
    cp -r biobambam2 $HAPSTER_DIR/bin && \
    ln -sfn biobambam2/x86_64-linux-gnu/2.0.146/bin/bamtofastq $HAPSTER_DIR/bin/bamtofastq && \
    tar xf hpseq-1.0.0.tar.gz && cd hpseq && make -j$NCORES deps && make -j$NCORES && \
    cp hpscan_cw $HAPSTER_DIR/bin && \
    cd $HAPSTER_DIR/build/tools-1.0.1 && \
    tar xf bwakit-0.7.15_x64-linux.tar.bz2 && \
    cp bwa.kit/k8 $HAPSTER_DIR/bin && cp bwa.kit/bwa-postalt.js $HAPSTER_DIR/bin && \
    unzip -q gatk-4.1.2.0.zip && \
    mv gatk-4.1.2.0/*sh $HAPSTER_DIR/bin && \
    mv gatk-4.1.2.0/*jar $HAPSTER_DIR/bin && \
    mv gatk-4.1.2.0/gatk $HAPSTER_DIR/bin' && \
    rm -rf /opt/hapster/build/tools-1.0.1 /opt/hapster/resources/tools-1.0.1.tar.gz

# Set working directory
WORKDIR /data

# Default entrypoint uses the conda environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "hapster"]
CMD ["hapster", "--help"]
