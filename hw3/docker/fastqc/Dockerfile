FROM --platform=linux/amd64 continuumio/miniconda3:24.7.1-0
RUN apt-get update \
    && apt-get install -y --no-install-recommends fontconfig fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/* \
    && conda install -y -c conda-forge -c bioconda fastqc=0.12.1 \
    && conda clean -ya
