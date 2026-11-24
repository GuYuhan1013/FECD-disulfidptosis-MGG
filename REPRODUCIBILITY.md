# Reproducibility check

This file records the clean-room rerun performed before the MGG submission. It is
not a results file; it documents what was executed and what was compared.

## Environment

- R 4.6.1 (2026-06-24), Windows 11 x64
- Python 3.13.18 (statistics and plotting of the optional S07 module)
- Project root contained only `02_原始数据/`; `S00_setup_directories.R` created
  the output tree before step 1.
- Random seed `20261121` is set in S03, S05_02, S06_03 and S06_08.

## What was run

`S00` through `S06_08b` (18 R scripts) were run in order on a copy of the raw
data. The optional S07 statistics and plotting scripts were run on the
checked-in 19-gene matrix. `S07_prepare_GSE74123.py` was not rerun because the
seven refine.bio PCL files and the annotation are not redistributed; its output
is checked in as `GSE74123_19gene_matrix.csv`.

## Key values reproduced

| Step | Value |
|---|---|
| S01 | 50,259 genes × 25 samples; 471,969,630 total counts; 9/4/6/6 by group |
| S03 | 16,754 filtered genes; 7,630 at adjusted P < 0.05; 3,964 DEGs (2,040 up, 1,924 down); 9/19 disulfidptosis genes significant; 0/19 with \|log2FC\| > 1 |
| S04 | 6,720 most variable genes; power 12; 16 modules merged to 6; turquoise r = 0.984, blue r = −0.963 |
| S05 | 3,543 GSEA sets tested; 909 at adjusted P < 0.05 (773 positive, 136 negative) |
| S06 | 9 candidates; GSE112201 direction concordance 5/9 (core 2/4); validation DEGs 687; Pre_S/REP/NR vs control 211/4,653/3,048 DEGs; energy axis Kruskal–Wallis P = 0.0004 |
| S07 | core 4/4, nine candidates 6/9, all 19 genes 11/19; r = 0.555, P = 0.0137 |

## Table comparison

- 30 of the 35 CSV tables in `04_分析结果/01_结果表格` were byte-identical to the
  archived versions.
- T09 and T10 (GO tables) and T22 (gene–pathway correlations) contained the same
  values but a different row order; the numeric columns matched exactly.
- T23 (single-gene GSEA) had identical enrichment scores (ES) and gene-set sizes;
  NES and Monte Carlo P values differed slightly before the seed was added to
  `S06_08`. With the seed, two consecutive runs produced byte-identical T23 files.
- `TableS2_GSE74123_validation.csv` was byte-identical to the archived file
  (SHA-256 `B7EE82E643ECF7225E9D4C7DA4D6178B27314D3EA9D4AA4FF2D2BC7D8187E32A`).

## Caveats

- `S05_01_moduleGO_KEGG.R` queries KEGG online. In the clean-room run the KEGG
  request timed out and the script skipped KEGG; the GO biological-process
  results used in the manuscript were unaffected.
- WGCNA and glmnet were built under R 4.6.2 and load with a version warning under
  R 4.6.1; this does not change the reported values.
- The optional S07 preparation step needs `pyreadr` and the raw refine.bio files.
  Steps 19–20 are reproducible from the checked-in matrix alone.
