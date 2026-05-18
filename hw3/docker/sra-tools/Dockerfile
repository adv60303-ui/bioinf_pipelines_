FROM --platform=linux/amd64 continuumio/miniconda3:24.7.1-0
RUN conda install -y -c conda-forge -c bioconda sra-tools=3.1.0 \
    && conda clean -ya
