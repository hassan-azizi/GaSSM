#!/usr/bin/env bash
set -euo pipefail

CUDA_HOME="${CUDA_HOME:-/usr/local/cuda}"
NVCC="${CUDA_HOME}/bin/nvcc"

if [[ ! -x "${NVCC}" ]]; then
    echo "nvcc not found at: ${NVCC}"
    exit 1
fi

"${NVCC}" --version

"${NVCC}" \
    -O3 \
    -arch=sm_89 \
    -I./source_code \
    ./source_code/RC_main.cu \
    -o GPU_string_polyatmoic \
    -lm

echo "Build successful: GPU_string_polyatmoic"