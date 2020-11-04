#!/bin/bash
source /root/.bashrc
conda activate bakta
if [[ $# -lt 1 ]]; then
    bakta --help
else
    bakta "$@"
fi
