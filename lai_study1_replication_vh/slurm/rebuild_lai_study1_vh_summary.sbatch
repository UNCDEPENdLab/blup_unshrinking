#!/usr/bin/env bash
#SBATCH --job-name=lai-study1-vh-aggregate
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_repo_root="$(cd "${script_dir}/../.." && pwd)"
submit_repo_root="${SLURM_SUBMIT_DIR:-}"
if [[ -n "${LAI_VH_REPO_ROOT:-}" ]]; then
  repo_root="${LAI_VH_REPO_ROOT}"
elif [[ -n "${submit_repo_root}" && -f "${submit_repo_root}/lai_study1_replication_vh/run.R" ]]; then
  repo_root="${submit_repo_root}"
else
  repo_root="${script_repo_root}"
fi
rscript="${LAI_VH_RSCRIPT:-Rscript}"
runner="${repo_root}/lai_study1_replication_vh/run.R"
n_sim="${LAI_VH_N_SIM:-5000}"
out_dir="${LAI_VH_OUT_DIR:-${repo_root}/outputs/lai_study1_vh_production}"
methods="${LAI_VH_METHODS:-oracle_dual,naive_dual_blup,closed_form_dual,fuller_closed_form,fuller_alpha_stepdown_closed_form,lai_2spa,msem}"

"${rscript}" "${runner}" \
  "${n_sim}" "${out_dir}" 1 NA 1 "${methods}" NA 1 1
