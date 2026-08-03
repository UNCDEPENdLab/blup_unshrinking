#!/usr/bin/env bash
#SBATCH --job-name=lai-study1-vh
#SBATCH --array=1-486%100
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# The seven-method bundle includes one Mplus fit per replication.  Keep the
# within-task core count at one; use the array throttle to respect the Mplus
# license and cluster capacity.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_repo_root="$(cd "${script_dir}/../.." && pwd)"
submit_repo_root="${SLURM_SUBMIT_DIR:-}"
if [[ -n "${LAI_VH_REPO_ROOT:-}" ]]; then
  repo_root="${LAI_VH_REPO_ROOT}"
elif [[ -n "${submit_repo_root}" && -f "${submit_repo_root}/lai_study1_replication_vh/run.R" ]]; then
  # Slurm may spool a copy of this script under /var/spool, so BASH_SOURCE is
  # not reliable after submission.  The submission directory remains stable.
  repo_root="${submit_repo_root}"
else
  repo_root="${script_repo_root}"
fi
runner="${repo_root}/lai_study1_replication_vh/run.R"
rscript="${LAI_VH_RSCRIPT:-Rscript}"

if [[ ! -f "${runner}" ]]; then
  echo "Could not find Lai Study 1 VH runner: ${runner}" >&2
  exit 2
fi

n_sim="${LAI_VH_N_SIM:-5000}"
out_dir="${LAI_VH_OUT_DIR:-${repo_root}/outputs/lai_study1_vh_production}"
methods="${LAI_VH_METHODS:-oracle_dual,naive_dual_blup,closed_form_dual,fuller_closed_form,fuller_alpha_stepdown_closed_form,lai_2spa,msem}"
condition_id="${LAI_VH_CONDITION_ID:-${SLURM_ARRAY_TASK_ID:?SLURM array task ID is required}}"
resume_existing="${LAI_VH_RESUME_EXISTING:-1}"

mkdir -p "${out_dir}/logs"
log_prefix="${out_dir}/logs/lai_study1_vh_${SLURM_JOB_ID:-manual}_${condition_id}"
exec > >(tee -a "${log_prefix}.out" "${log_prefix}.log") \
     2> >(tee -a "${log_prefix}.err" "${log_prefix}.log" >&2)

echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Condition: ${condition_id}; replications: ${n_sim}; cores: 1"
echo "Methods: ${methods}"
echo "Output: ${out_dir}"

# Arguments 7--9: condition ID, aggregate-only flag, write-aggregate flag.
# Array tasks never write the shared aggregate summary file.
"${rscript}" "${runner}" \
  "${n_sim}" "${out_dir}" 1 NA "${resume_existing}" "${methods}" \
  "${condition_id}" 0 0

echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
