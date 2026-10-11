#!/bin/sh
# Run bakta by default: docker run <image> --db <db> genome.fna
# Run other commands as given, e.g. when nextflow runs: docker run <image> /bin/bash -c ...
case "$1" in
    /bin/*|/usr/bin/*|/usr/local/bin/*|/opt/conda/bin/*)
        exec "$@" ;;
    ''|-*|*/*)
        ;;
    *)
        if command -v "$1" >/dev/null 2>&1; then
            exec "$@"
        fi ;;
esac
exec bakta "$@"
