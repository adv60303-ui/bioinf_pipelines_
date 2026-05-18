FROM --platform=linux/amd64 continuumio/miniconda3:24.7.1-0
RUN conda install -y -c conda-forge -c bioconda samtools=1.19 matplotlib=3.8 \
    && conda clean -ya
