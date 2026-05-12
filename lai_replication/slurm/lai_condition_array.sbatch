#!/usr/bin/env bash
#SBATCH --job-name=lai-condition-array
#SBATCH --array=1-117
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=lai_condition_array_%A_%a.out
#SBATCH --error=lai_condition_array_%A_%a.err

set -euo pipefail

# Submit from the repository root, or set LAI_REPO_ROOT explicitly.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="${LAI_REPO_ROOT:-$(cd "${script_dir}/../.." && pwd)}"
runner="${repo_root}/lai_replication/mlm_random_slope_lai_apples_to_apples_sim.R"
rscript="${LAI_RSCRIPT:-Rscript}"

if [[ ! -f "${runner}" ]]; then
  echo "Could not find Lai runner: ${runner}" >&2
  exit 2
fi

n_sim="${LAI_N_SIM:-1000}"
study="${LAI_STUDY:-all}"
out_dir="${LAI_OUT_DIR:-${repo_root}/outputs/lai_apples_to_apples_slurm}"
n_cores="${LAI_N_CORES:-${SLURM_CPUS_PER_TASK:-1}}"
max_conditions="${LAI_MAX_CONDITIONS:-NA}"
chunk_index="${LAI_CHUNK_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}"
chunk_size="${LAI_CHUNK_SIZE:-5}"
resume_existing="${LAI_RESUME_EXISTING:-1}"
include_tempered_eiv="${LAI_INCLUDE_TEMPERED_EIV:-0}"

mkdir -p "${out_dir}/logs"
log_path="${out_dir}/logs/lai_condition_array_${SLURM_JOB_ID:-manual}_${chunk_index}.log"

{
  echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "Host: $(hostname)"
  echo "Repository: ${repo_root}"
  echo "Rscript: ${rscript}"
  echo "Output directory: ${out_dir}"
  echo "Study selector: ${study}"
  echo "Replications per condition: ${n_sim}"
  echo "Chunk index: ${chunk_index}"
  echo "Chunk size: ${chunk_size}"
  echo "Cores: ${n_cores}"
  echo "Max conditions: ${max_conditions}"
  echo "Resume existing: ${resume_existing}"
  echo "Include tempered EIV: ${include_tempered_eiv}"
} | tee "${log_path}"

"${rscript}" "${runner}" \
  "${n_sim}" \
  "${study}" \
  "${out_dir}" \
  "${n_cores}" \
  "${max_conditions}" \
  "${chunk_index}" \
  "${chunk_size}" \
  "${resume_existing}" \
  "${include_tempered_eiv}" 2>&1 | tee -a "${log_path}"

echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${log_path}"
