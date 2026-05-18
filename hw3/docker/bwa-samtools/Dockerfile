FROM --platform=linux/amd64 continuumio/miniconda3:24.7.1-0
RUN conda install -y -c conda-forge -c bioconda bwa=0.7.17 samtools=1.19 \
    && conda clean -ya
