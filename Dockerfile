FROM mambaorg/micromamba:2-debian13-slim AS build

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/source/
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && /opt/conda/bin/python -m pip install --no-cache-dir --no-deps /tmp/source/ \
    && micromamba clean --all --force-pkgs-dirs --yes \
    && find /opt/conda \( -name '*.a' -o -name '*.a.*' \) -delete \
    && find /opt/conda/lib/python3*/site-packages -type d -name tests -prune -exec rm -rf {} + \
    && rm -rf /opt/conda/pkgs /opt/conda/include /opt/conda/man /opt/conda/share/doc /opt/conda/share/man /opt/conda/share/info /opt/conda/share/locale /opt/conda/share/easel /opt/conda/etc/conda/test-files


FROM debian:trixie-slim

LABEL org.opencontainers.image.authors="oliver.schwengers@computational.bio.uni-giessen.de,lukas.jelonek@computational.bio.uni-giessen.de"
LABEL org.opencontainers.image.url='https://github.com/oschwengers/bakta'
LABEL org.opencontainers.image.documentation='https://github.com/oschwengers/bakta/readme.md'
LABEL org.opencontainers.image.title='Bakta'
LABEL org.opencontainers.image.description='Rapid & standardized annotation of bacterial genomes, MAGs & plasmids'

COPY --from=build --chown=root:root /opt/conda /opt/conda
COPY --from=build /usr/bin/micromamba /usr/bin/micromamba
# keep the env on the PATH in login shells, where /etc/profile resets it
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh

# the conda env needs no activation scripts, so putting it on the PATH enables it
# for every process, including --entrypoint overrides and workflow managers like nextflow
ENV PATH=/opt/conda/bin:$PATH \
    CONDA_PREFIX=/opt/conda \
    MAMBA_ROOT_PREFIX=/opt/conda \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    MPLCONFIGDIR=/tmp/matplotlib

COPY --chmod=755 docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
ENTRYPOINT ["docker-entrypoint.sh"]
