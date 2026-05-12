#!/usr/bin/env bash
#SBATCH --job-name=lai-condition-array
#SBATCH --array=1-117
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Submit from the repository root, or set LAI_REPO_ROOT explicitly.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_repo_root="$(cd "${script_dir}/../.." && pwd)"
submit_repo_root="${SLURM_SUBMIT_DIR:-}"

if [[ -n "${LAI_REPO_ROOT:-}" ]]; then
  repo_root="${LAI_REPO_ROOT}"
elif [[ -n "${submit_repo_root}" && -f "${submit_repo_root}/lai_replication/mlm_random_slope_lai_apples_to_apples_sim.R" ]]; then
  repo_root="${submit_repo_root}"
else
  repo_root="${script_repo_root}"
fi

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
log_prefix="${out_dir}/logs/lai_condition_array_${SLURM_JOB_ID:-manual}_${chunk_index}"
stdout_path="${log_prefix}.out"
stderr_path="${log_prefix}.err"
log_path="${log_prefix}.log"

# SLURM parses #SBATCH --output/--error before shell variables such as
# LAI_OUT_DIR exist. Redirect after computing out_dir so each run keeps its
# array stdout/stderr beside the combined wrapper log.
exec > >(tee -a "${stdout_path}" "${log_path}") 2> >(tee -a "${stderr_path}" "${log_path}" >&2)

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
echo "Stdout log: ${stdout_path}"
echo "Stderr log: ${stderr_path}"
echo "Combined log: ${log_path}"

"${rscript}" "${runner}" \
  "${n_sim}" \
  "${study}" \
  "${out_dir}" \
  "${n_cores}" \
  "${max_conditions}" \
  "${chunk_index}" \
  "${chunk_size}" \
  "${resume_existing}" \
  "${include_tempered_eiv}"

echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
