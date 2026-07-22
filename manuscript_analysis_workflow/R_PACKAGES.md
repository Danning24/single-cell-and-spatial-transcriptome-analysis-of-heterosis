# R package notes

The Conda environment covers the broadly distributed dependencies. Install the following project-specific packages from their official repositories when they are not available through Conda/Bioconductor:

- monocle3 and SeuratWrappers for pseudotime conversion and trajectory analysis;
- UCell for rank-based spatial signature scoring;
- PlantCellChat for plant ligand-receptor inference;
- scTenifoldKnk and scTenifoldNet for virtual knockout and network construction;
- zellkonverter for H5AD export.

Package versions used for the final ROS integration are recorded in the root README. Earlier exploratory stages were recovered from their original scripts; use `SOURCE_PROVENANCE.tsv` to trace each public file to its source snapshot.
