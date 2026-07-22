#!/bin/bash
set -euo pipefail

cd ${MAIZE_PROJECT_ROOT:?Set MAIZE_PROJECT_ROOT}

mkdir -p new_ko_genes_scripts/double_ko/logs

prep_job=$(sbatch --parsable new_ko_genes_scripts/double_ko/prepare_double_ko_inputs.sbatch)
echo "Submitted prepare job: ${prep_job}"

for group in B73xA554 B73xOh40B; do
  job=$(sbatch --parsable --dependency=afterok:${prep_job} \
    --job-name="double_ko_${group}" \
    new_ko_genes_scripts/double_ko/run_double_ko.sbatch "${group}")
  echo "Submitted ${group}: ${job}"
done
