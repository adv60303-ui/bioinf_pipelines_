FROM --platform=linux/amd64 continuumio/miniconda3:24.7.1-0
RUN conda install -y -c conda-forge -c bioconda megahit=1.2.9 \
    && conda clean -ya
