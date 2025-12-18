# Verification Branch

This branch exists to verify the CUDA implementation (cuAleae) against the original CPU implementation on the same chemical reaction network inputs. It adds deterministic seeding, CSV output, and helper scripts to run paired tests and analyze statistical equivalence.

## What it does
- Builds two binaries:
  - CUDA: `cuAleae`
  - CPU: `old_src_code/aleae`
- Runs both for a set of seeds with identical inputs/params, emits per-chemical stats (mean, variance, threshold hits) to CSVs.
- Provides an analysis script (`analyze.py`) to merge CPU/GPU results, compute t-tests, and generate a scatter plot.

## How to run
1) Build:
   - CUDA: `make clean && make build`
   - CPU:  `make -C old_src_code clean && make -C old_src_code`

2) Generate results (20 seeds, 5,000 trials on test inputs):

`./run_seeds.sh`

This writes `results/gpu_seed*.csv` and `results/cpu_seed*.csv`.

If you want a different seed range or trials, edit `run_seeds.sh`.

3) Analyze:

`python3 analyze.py`

Prints a summary table and saves `results/mean_scatter.png`.

## Notes
- Inputs live in `input_files/` (e.g., `test.in`, `test.r`).
- CLI flags `--seed` and `--csv` are available on both binaries for ad hoc runs.
- Deterministic seeding ensures paired CPU/GPU runs are comparable.