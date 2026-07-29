#!/bin/bash

set -euo pipefail

cd ${MAIZE_PROJECT_ROOT:?Set MAIZE_PROJECT_ROOT}

JOB_TABLE="new_ko_genes_scripts/sctenifold_inputs/job_table.tsv"

if [ ! -f "$JOB_TABLE" ]; then
  echo "Cannot find job table: $JOB_TABLE" >&2
  exit 1
fi

mkdir -p new_ko_genes_scripts/sctenifold_logs new_ko_genes_scripts/sctenifold_results

tail -n +2 "$JOB_TABLE" | while IFS=$'\t' read -r group gene present runner; do
  if [ "$present" != "TRUE" ]; then
    echo "SKIP missing gene in matrix: $group $gene"
    continue
  fi

  if [ "$group" = "B73xA554" ]; then
    sbatch --job-name="safe_${group}_${gene}" new_ko_genes_scripts/run_sctenifold_knk_safe.sbatch "$group" "$gene"
  else
    sbatch --job-name="${group}_${gene}" new_ko_genes_scripts/run_sctenifold_knk.sbatch "$group" "$gene"
  fi
done
