# S07 independent validation (GSE74123)

This optional module reproduces the independent direction check added to the
MGG manuscript. It is not part of the 17-step R pipeline.

## Data

- GEO series matrix: `GSE74123_series_matrix.txt.gz`
  `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74123/matrix/GSE74123_series_matrix.txt.gz`
- Seven processed samples from refine.bio (SCAN log2): `GSM1910841`–`GSM1910847`
  Experiment endpoint: `https://api.refine.bio/v1/experiments/GSE74123`
- `GSM1910840` (FECD A13) was not processed by refine.bio. It is recovered from
  the series matrix by probe correlation, as described below.
- Annotation: `AnnoProbe::humanGTF` (`anno/AnnoProbe/R/sysdata.rda`), or any
  ENSEMBL-to-SYMBOL table with the same columns.

## Samples

| Group | Samples |
|---|---|
| Control (4) | GSM1910841 (A17), GSM1910843 (A25), GSM1910845 (A28), GSM1910847 (A31) |
| FECD (4) | GSM1910840 (A13), GSM1910842 (A22), GSM1910844 (A39), GSM1910846 (A41) |

## Steps

All three Python scripts resolve paths relative to this folder, so they can be
run from any directory with `python <path-to-script>`. Step 1 requires
`pyreadr`; steps 2–3 do not.

1. Place the series matrix, the seven `GSM*.PCL` files and the annotation in
   this folder.
2. Run the preparation script:

   `python S07_prepare_GSE74123.py`

   This writes `GSE74123_19gene_matrix.csv`. For each gene, the probe with the
   highest absolute correlation to the refine.bio gene vector across the seven
   processed samples is selected (`|r| > 0.98`); if r is negative, the probe row
   is multiplied by −1 so that its direction matches the processed vector, and
   used across all eight samples.
3. Run the statistics script:

   `python S07_GSE74123_validation.py`

   This writes `TableS2_GSE74123_validation.csv`.
4. Run the plotting script:

   `python S07_plot_FigS4.py`

   This writes `FigS4_GSE74123_validation.pdf/.png` from the checked-in matrix
   and TableS2.

The checked-in `GSE74123_19gene_matrix.csv` and
`TableS2_GSE74123_validation.csv` let step 3 be run without downloading the
raw files.

## Scope

Four samples per group. The analysis is a direction check, not a replication of
differential-expression counts. Nominal Welch and Mann–Whitney P values are
reported without multiplicity adjustment.
