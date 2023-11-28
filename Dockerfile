FROM alpine:3.12

LABEL org.opencontainers.image.authors="oliver.schwengers@computational.bio.uni-giessen.de,lukas.jelonek@computational.bio.uni-giessen.de"
LABEL org.opencontainers.image.url='https://github.com/oschwengers/bakta'
LABEL org.opencontainers.image.documentation='https://github.com/oschwengers/bakta/readme.md'
LABEL org.opencontainers.image.title='Bakta'
LABEL org.opencontainers.image.description='Rapid & standardized annotation of bacterial genomes, MAGs & plasmids'

RUN apk update && apk add wget tar bash \
    && wget -q -O /etc/apk/keys/sgerrand.rsa.pub https://alpine-pkgs.sgerrand.com/sgerrand.rsa.pub \
    && wget https://github.com/sgerrand/alpine-pkg-glibc/releases/download/2.32-r0/glibc-2.32-r0.apk \
    && apk add glibc-2.32-r0.apk \
    && rm glibc-2.32-r0.apk \
    && wget --no-iri -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba \
    && touch /root/.bashrc \
    && ./bin/micromamba shell init -s bash -p /opt/conda  \
    && cp /root/.bashrc /opt/conda/bashrc

COPY environment.yml /tmp/
RUN echo -e '\n  - deepsig>=1.2.5' >> /tmp/environment.yml

SHELL ["bash", "-l" ,"-c"]

RUN source /opt/conda/bashrc && micromamba activate \
    && micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

COPY . /tmp/source/

RUN source /opt/conda/bashrc && micromamba activate \
    && python3 -m pip install --no-cache /tmp/source/ \
    && echo '#!/bin/bash' > /entrypoint.sh \
    && echo 'bakta "$@"' >> /entrypoint.sh \
    && chmod +x /entrypoint.sh \
    \
    # replace bash with a wrapper that initializes the micromamba env
    # everytime it is not initialized.
    # with this approach we can initialize the env automatically, even
    # in non-interactive, no-login sessions, e.g. when nextflow uses
    # containers
    && mv /bin/bash /bin/bash.orig \
    && echo '#!/bin/bash.orig' >> /bin/bash \
    && echo 'if [[ -z $MAMBA_INITIALIZED ]]' >> /bin/bash \
    && echo 'then' >> /bin/bash \
    && echo 'source /opt/conda/bashrc' >> /bin/bash \
    && echo 'micromamba activate' >> /bin/bash \
    && echo 'export MAMBA_INITIALIZED=1' >> /bin/bash \
    && echo 'fi' >> /bin/bash \
    && echo '/bin/bash.orig "$@"' >> /bin/bash \
    && chmod +x /bin/bash

ENTRYPOINT ["/entrypoint.sh"]
