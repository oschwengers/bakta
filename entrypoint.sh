#!/bin/bash
source /opt/conda/bashrc
if [[ $# -lt 1 ]]; then
    bakta --help
else
    bakta "$@"
fi
