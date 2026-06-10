#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_repo_root="$(cd "${script_dir}/../.." && pwd)"
submit_repo_root="${SLURM_SUBMIT_DIR:-}"

if [[ -n "${BLUP_OUTCOME_REPO_ROOT:-}" ]]; then
  repo_root="${BLUP_OUTCOME_REPO_ROOT}"
elif [[ -n "${submit_repo_root}" && -f "${submit_repo_root}/blup_outcome/mlm_random_slope_blup_outcome_sim.R" ]]; then
  repo_root="${submit_repo_root}"
else
  repo_root="${script_repo_root}"
fi

job_script="${repo_root}/blup_outcome/slurm/blup_outcome_condition_array.sbatch"
rscript="${BLUP_OUTCOME_RSCRIPT:-Rscript}"
grid_mode="${BLUP_OUTCOME_GRID_MODE:-base}"
chunk_size="${BLUP_OUTCOME_CHUNK_SIZE:-5}"
max_conditions="${BLUP_OUTCOME_MAX_CONDITIONS:-NA}"

if [[ ! -f "${job_script}" ]]; then
  echo "Could not find BLUP-outcome SLURM job script: ${job_script}" >&2
  exit 2
fi

array_end="$(
  BLUP_OUTCOME_REPO_ROOT="${repo_root}" \
  BLUP_OUTCOME_GRID_MODE="${grid_mode}" \
  BLUP_OUTCOME_CHUNK_SIZE="${chunk_size}" \
  BLUP_OUTCOME_MAX_CONDITIONS="${max_conditions}" \
  "${rscript}" --vanilla - <<'RS'
suppressPackageStartupMessages({
  library(tibble)
  library(tidyr)
  library(dplyr)
})

repo_root <- Sys.getenv("BLUP_OUTCOME_REPO_ROOT")
source(file.path(repo_root, "blup_outcome", "designs.R"), local = TRUE)

grid_mode <- Sys.getenv("BLUP_OUTCOME_GRID_MODE", "base")
chunk_size <- Sys.getenv("BLUP_OUTCOME_CHUNK_SIZE", "5")
max_conditions <- Sys.getenv("BLUP_OUTCOME_MAX_CONDITIONS", "NA")

cat(blup_outcome_chunk_count(grid_mode = grid_mode, chunk_size = as.integer(chunk_size), max_conditions = max_conditions))
RS
)"

echo "Submitting BLUP-outcome array with ${array_end} tasks for grid_mode=${grid_mode}, chunk_size=${chunk_size}, max_conditions=${max_conditions}."

exec sbatch \
  --export=ALL \
  --array=1-"${array_end}" \
  "$@" \
  "${job_script}"