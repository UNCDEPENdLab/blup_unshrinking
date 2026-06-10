#!/usr/bin/env bash
#SBATCH --job-name=blup-outcome-array
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# The array range depends on the selected design grid and chunk size, so it is
# supplied at submission time by the companion launcher or via `sbatch --array`.

# Submit from the repository root, or set BLUP_OUTCOME_REPO_ROOT explicitly.
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

runner="${repo_root}/blup_outcome/mlm_random_slope_blup_outcome_sim.R"
rscript="${BLUP_OUTCOME_RSCRIPT:-Rscript}"

if [[ ! -f "${runner}" ]]; then
  echo "Could not find BLUP-outcome runner: ${runner}" >&2
  exit 2
fi

n_sim="${BLUP_OUTCOME_N_SIM:-100}"
out_dir="${BLUP_OUTCOME_OUT_DIR:-${repo_root}/outputs/blup_outcome_slurm}"
n_cores="${BLUP_OUTCOME_N_CORES:-${SLURM_CPUS_PER_TASK:-1}}"
derivative_method="${BLUP_OUTCOME_DERIVATIVE_METHOD:-handcoded}"
grid_mode="${BLUP_OUTCOME_GRID_MODE:-base}"
analysis_mode="${BLUP_OUTCOME_ANALYSIS_MODE:-full}"
max_conditions="${BLUP_OUTCOME_MAX_CONDITIONS:-NA}"
chunk_index="${BLUP_OUTCOME_CHUNK_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}"
chunk_size="${BLUP_OUTCOME_CHUNK_SIZE:-5}"
resume_existing="${BLUP_OUTCOME_RESUME_EXISTING:-1}"
execution_mode="${BLUP_OUTCOME_EXECUTION_MODE:-run}"

mkdir -p "${out_dir}/logs"
log_prefix="${out_dir}/logs/blup_outcome_condition_array_${SLURM_JOB_ID:-manual}_${chunk_index}"
stdout_path="${log_prefix}.out"
stderr_path="${log_prefix}.err"
log_path="${log_prefix}.log"

# SLURM parses #SBATCH --output/--error before shell variables exist. Redirect
# after computing out_dir so each run keeps its logs beside the combined log.
exec > >(tee -a "${stdout_path}" "${log_path}") 2> >(tee -a "${stderr_path}" "${log_path}" >&2)

echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Host: $(hostname)"
echo "Repository: ${repo_root}"
echo "Rscript: ${rscript}"
echo "Output directory: ${out_dir}"
echo "Grid mode: ${grid_mode}"
echo "Analysis mode: ${analysis_mode}"
echo "Derivative backend: ${derivative_method}"
echo "Replications per condition: ${n_sim}"
echo "Chunk index: ${chunk_index}"
echo "Chunk size: ${chunk_size}"
echo "Cores: ${n_cores}"
echo "Max conditions: ${max_conditions}"
echo "Resume existing: ${resume_existing}"
echo "Execution mode: ${execution_mode}"
echo "Stdout log: ${stdout_path}"
echo "Stderr log: ${stderr_path}"
echo "Combined log: ${log_path}"

"${rscript}" "${runner}" \
  "${n_sim}" \
  "${out_dir}" \
  "${n_cores}" \
  "${derivative_method}" \
  "${grid_mode}" \
  "${analysis_mode}" \
  "${chunk_index}" \
  "${chunk_size}" \
  "${resume_existing}" \
  "${execution_mode}" \
  "${max_conditions}"

echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
