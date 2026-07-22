#!/usr/bin/env bash

input_dir="./data_output/cortex_h5ad"
output_root="./data_output/mmd_result"

mkdir -p "$output_root"

for file in "$input_dir"/*.hvg2000.pca10.h5ad
do
    name=$(basename "$file" .hvg2000.pca10.h5ad)

    out_dir="$output_root/$name"
    mkdir -p "$out_dir"

    echo "======================================"
    echo "Running scmmd for $name"
    echo "Input: $file"
    echo "Output: $out_dir"
    echo "======================================"

    scmmd \
        --data "$file" \
        --out_path "$out_dir" \
        --groupby group_mmd \
        --contrast celltype \
        --representation X_pca \
        --sample_size 100\
        --n_iter 100

done

echo "All jobs finished."
