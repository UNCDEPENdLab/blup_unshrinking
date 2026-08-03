#!/usr/bin/env bash
#SBATCH --job-name=lai-study1-vh-figures
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_repo_root="$(cd "${script_dir}/../.." && pwd)"
submit_repo_root="${SLURM_SUBMIT_DIR:-}"
if [[ -n "${LAI_VH_REPO_ROOT:-}" ]]; then
  repo_root="${LAI_VH_REPO_ROOT}"
elif [[ -n "${submit_repo_root}" && -f "${submit_repo_root}/lai_study1_replication_vh/make_figure2_analogue.R" ]]; then
  repo_root="${submit_repo_root}"
else
  repo_root="${script_repo_root}"
fi

rscript="${LAI_VH_RSCRIPT:-Rscript}"
results_dir="${LAI_VH_OUT_DIR:-${repo_root}/outputs/lai_study1_vh_production}"
analysis_dir="${LAI_VH_ANALYSIS_DIR:-${results_dir}/figure2_analogue}"
force="${LAI_VH_FIGURE_FORCE:-0}"

mkdir -p "${analysis_dir}/logs"
log_prefix="${analysis_dir}/logs/lai_study1_vh_figures_${SLURM_JOB_ID:-manual}"
exec > >(tee -a "${log_prefix}.out" "${log_prefix}.log") \
     2> >(tee -a "${log_prefix}.err" "${log_prefix}.log" >&2)

echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results: ${results_dir}"
echo "Analysis: ${analysis_dir}"
echo "Force historical cache rebuild: ${force}"

"${rscript}" "${repo_root}/lai_study1_replication_vh/make_figure2_analogue.R" \
  "${results_dir}" "${analysis_dir}" "${force}"

echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
