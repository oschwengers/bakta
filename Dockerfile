FROM ubuntu:20.04

LABEL maintainer=TODO
LABEL description=TODO
LABEL version=TODO

ENV PATH=$PATH:/opt/conda/bin
RUN apt update && apt install -y wget \
    && wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda -u \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y mamba

# copy the conda env first, so that it can be setup and the same layer can used even when the other source files change
COPY environment.yml /tmp/source/
RUN mamba env create --file  /tmp/source/environment.yml
COPY . /tmp/source/
RUN python3 -m pip install /tmp/source/

ENTRYPOINT ["/opt/conda/bin/bakta"]
CMD ["/opt/conda/bin/bakta", "--help"]