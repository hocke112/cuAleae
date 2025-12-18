#!/usr/bin/env bash
set -euo pipefail

# Make sure output directory exists
mkdir -p results

# Loop over seeds 1 through 20
for s in {1..20}; do
  echo "Seed $s (CUDA)…"
  ./cuAleae input_files/test.in input_files/test.r 5000 -1 0 -1 --seed "$s" --csv "results/gpu_seed${s}.csv"

  echo "Seed $s (CPU)…"
  ./old_src_code/aleae input_files/test.in input_files/test.r 5000 -1 0 -1 --seed "$s" --csv "results/cpu_seed${s}.csv"
done