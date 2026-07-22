#!/bin/bash
# ============================================================
# submit_jobs.sh
#
# 作用：
#   一次性提交新基因的全部 scTenifoldKnk 虚拟敲除任务。
#
# 任务数量：
#   job_table.tsv 中每一行就是一个任务。
#   当前应为 22 genes x 2 groups = 44 jobs。
#
# 提交规则：
#   B73xA554  -> run_sctenifold_knk_safe.sbatch
#   B73xOh40B -> run_sctenifold_knk.sbatch
#
# 使用方式：
#   cd ${MAIZE_PROJECT_ROOT:?Set MAIZE_PROJECT_ROOT}
#   bash new_ko_genes_scripts/submit_jobs.sh
# ============================================================

set -euo pipefail

cd ${MAIZE_PROJECT_ROOT:?Set MAIZE_PROJECT_ROOT}

# job table 由 prepare_inputs_annotated.R 生成。
# 列为：group, gene, present_in_matrix, suggested_runner
JOB_TABLE="new_ko_genes_scripts/sctenifold_inputs/job_table.tsv"

if [ ! -f "$JOB_TABLE" ]; then
  echo "Cannot find job table: $JOB_TABLE" >&2
  exit 1
fi

mkdir -p new_ko_genes_scripts/sctenifold_logs new_ko_genes_scripts/sctenifold_results

# 跳过表头，从第二行开始逐行读取任务。
tail -n +2 "$JOB_TABLE" | while IFS=$'\t' read -r group gene present runner; do
  # 如果某个基因没有进入输入矩阵，就不提交。
  if [ "$present" != "TRUE" ]; then
    echo "SKIP missing gene in matrix: $group $gene"
    continue
  fi

  # B73xA554 之前标准 manifoldAlignment 会报维度退化错误，使用 safe runner。
  if [ "$group" = "B73xA554" ]; then
    sbatch --job-name="safe_${group}_${gene}" new_ko_genes_scripts/run_sctenifold_knk_safe.sbatch "$group" "$gene"
  else
    sbatch --job-name="${group}_${gene}" new_ko_genes_scripts/run_sctenifold_knk.sbatch "$group" "$gene"
  fi
done
