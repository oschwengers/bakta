#!/usr/bin/env bash

set -euo pipefail

ORI_FINDER_ARCH_URL='https://tubic.tju.edu.cn/Ori-Finder-Arch/download/Ori-Finder-Arch.tar.gz'
ORI_FINDER_ARCH_SHA256='f6a29d446b82df93b28144c02d404eeee6a3ba83fd9e0a0f31f63ba6b606b703'
target_dir="${1:-${CONDA_PREFIX:?Provide an installation directory or activate the Conda environment}/bin}"
tmp_dir=$(mktemp -d)
trap 'rm -rf "${tmp_dir}"' EXIT

mkdir -p "${target_dir}"
wget -q -O "${tmp_dir}/Ori-Finder-Arch.tar.gz" "${ORI_FINDER_ARCH_URL}"
echo "${ORI_FINDER_ARCH_SHA256}  ${tmp_dir}/Ori-Finder-Arch.tar.gz" | sha256sum --check --status
tar -xzf "${tmp_dir}/Ori-Finder-Arch.tar.gz" -C "${tmp_dir}"
install -m 0755 "${tmp_dir}/OriFinderArch" "${target_dir}/OriFinderArch"
