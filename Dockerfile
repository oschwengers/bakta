FROM ubuntu:20.04

LABEL maintainer=TODO
LABEL description=TODO
LABEL version=TODO

# Move bashrc aside as the default bashrc won't be loaded in non-interactive mode
# but conda requires the file to be loaded. We create a emtpy one for this purpose.
# At the end we will replace the temporary bashrc with the original one and reregister
# conda to it
RUN mv /root/.bashrc /root/.bashrc.interative && touch /root/.bashrc && \
    apt update && apt install -y wget \
    && wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda -u \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && /opt/conda/bin/conda config --add channels defaults \
    && /opt/conda/bin/conda config --add channels bioconda \
    && /opt/conda/bin/conda config --add channels conda-forge \
    && /opt/conda/bin/conda install -y mamba \
    && /opt/conda/bin/conda init bash

# copy the conda env first, so that it can be setup and the same layer can used even when the other source files change
SHELL ["/bin/bash", "-l", "-c"]
COPY environment.yml /tmp/source/
RUN mamba env create --file  /tmp/source/environment.yml

COPY . /tmp/source/
RUN source /root/.bashrc && conda activate bakta && python3 -m pip install /tmp/source/

COPY entrypoint.sh /
ENTRYPOINT ["/entrypoint.sh"]
