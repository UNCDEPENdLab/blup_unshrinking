#!/usr/bin/env bash
#SBATCH --job-name=score-correction-condition-array
#SBATCH --array=1-117
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=20:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Submit from the repository root, or set SCORE_CORRECTION_REPO_ROOT explicitly.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_repo_root="$(cd "${script_dir}/../.." && pwd)"
submit_repo_root="${SLURM_SUBMIT_DIR:-}"

if [[ -n "${SCORE_CORRECTION_REPO_ROOT:-}" ]]; then
  repo_root="${SCORE_CORRECTION_REPO_ROOT}"
elif [[ -n "${submit_repo_root}" && -f "${submit_repo_root}/score_correction_methods/score_correction_comparison.R" ]]; then
  repo_root="${submit_repo_root}"
else
  repo_root="${script_repo_root}"
fi

runner="${repo_root}/score_correction_methods/score_correction_comparison.R"
rscript="${SCORE_CORRECTION_RSCRIPT:-Rscript}"

if [[ ! -f "${runner}" ]]; then
  echo "Could not find score correction runner: ${runner}" >&2
  exit 2
fi

n_sim="${SCORE_CORRECTION_N_SIM:-1000}"
out_dir="${SCORE_CORRECTION_OUT_DIR:-${repo_root}/outputs/score_correction_comparison_slurm}"
n_cores="${SCORE_CORRECTION_N_CORES:-${SLURM_CPUS_PER_TASK:-1}}"
max_conditions="${SCORE_CORRECTION_MAX_CONDITIONS:-NA}"
chunk_index="${SCORE_CORRECTION_CHUNK_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}"
chunk_size="${SCORE_CORRECTION_CHUNK_SIZE:-5}"
resume_existing="${SCORE_CORRECTION_RESUME_EXISTING:-1}"

mkdir -p "${out_dir}/logs"
log_prefix="${out_dir}/logs/score_correction_condition_array_${SLURM_JOB_ID:-manual}_${chunk_index}"
stdout_path="${log_prefix}.out"
stderr_path="${log_prefix}.err"
log_path="${log_prefix}.log"

# SLURM parses #SBATCH --output/--error before shell variables such as
# SCORE_CORRECTION_OUT_DIR exist. Redirect after computing out_dir so each run keeps its
# array stdout/stderr beside the combined wrapper log.
exec > >(tee -a "${stdout_path}" "${log_path}") 2> >(tee -a "${stderr_path}" "${log_path}" >&2)

echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Host: $(hostname)"
echo "Repository: ${repo_root}"
echo "Rscript: ${rscript}"
echo "Output directory: ${out_dir}"
echo "Replications per condition: ${n_sim}"
echo "Chunk index: ${chunk_index}"
echo "Chunk size: ${chunk_size}"
echo "Cores: ${n_cores}"
echo "Max conditions: ${max_conditions}"
echo "Resume existing: ${resume_existing}"
echo "Stdout log: ${stdout_path}"
echo "Stderr log: ${stderr_path}"
echo "Combined log: ${log_path}"

"${rscript}" "${runner}" \
  "${n_sim}" \
  "${out_dir}" \
  "${n_cores}" \
  "${max_conditions}" \
  "${chunk_index}" \
  "${chunk_size}" \
  "${resume_existing}"

echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
