# -*- coding: utf-8 -*-
"""S07 preparation: build a 4 FECD / 4 control gene-level matrix from GSE74123.

Required files (place in this folder):
  GSE74123_series_matrix.txt.gz
  GSM1910841_hyb8262.PCL ... GSM1910847_hyb8268.PCL  (refine.bio SCAN log2)
  anno/AnnoProbe/R/sysdata.rda                       (humanGTF)

The GEO series matrix probe IDs do not match the GPL6244 platform table.
We therefore orient each gene by the probe that correlates most strongly with
the refine.bio gene vector across the seven processed samples (|r| > 0.98). If
r is negative, the entire probe row is multiplied by -1 so that its direction
matches the processed vector; the resulting direction-oriented values are not
negative RMA. This recovers GSM1910840 (FECD A13), which refine.bio did not
process.
"""
import glob
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr

HERE = Path(__file__).resolve().parent

PROCESSED = {
    "GSM1910841": "control_A17", "GSM1910843": "control_A25",
    "GSM1910845": "control_A28", "GSM1910847": "control_A31",
    "GSM1910842": "FECD_A22", "GSM1910844": "FECD_A39",
    "GSM1910846": "FECD_A41",
}

# Ensembl -> symbol map from AnnoProbe::humanGTF.
h = pyreadr.read_r(HERE / "anno/AnnoProbe/R/sysdata.rda")["humanGTF"]
h = h.dropna(subset=["ENSEMBL", "SYMBOL"])[["ENSEMBL", "SYMBOL"]].drop_duplicates()
h["ENSEMBL"] = h["ENSEMBL"].str.replace(r"\..*$", "", regex=True)
ens2sym = dict(zip(h["ENSEMBL"], h["SYMBOL"]))

# Read the seven processed PCL files; collapse multiple Ensembl IDs per symbol
# by the median.
pcl = {}
for gsm in PROCESSED:
    f = glob.glob(str(HERE / f"{gsm}_*.PCL"))[0]
    d = pd.read_csv(f, sep="\t", header=None, skiprows=1, names=["ens", "value"])
    d["ens"] = d["ens"].str.replace("_at$", "", regex=True)
    d["symbol"] = d["ens"].map(ens2sym)
    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    d = d.dropna(subset=["symbol", "value"])
    pcl[gsm] = d.groupby("symbol")["value"].median()
pcl = pd.DataFrame(pcl)

# Read the GEO series matrix.
rows, samples = [], None
with gzip.open(HERE / "GSE74123_series_matrix.txt.gz", "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            samples = [x.strip().strip('"') for x in next(fh).rstrip("\n").split("\t")[1:]]
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if line.startswith("!"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 9:
            rows.append([parts[0].strip('"')] + [float(x) for x in parts[1:9]])
mat = pd.DataFrame(rows, columns=["probe"] + samples).set_index("probe")

# Probe orientation and extraction.
overlap = list(PROCESSED)
all_samples = ["GSM1910840"] + list(PROCESSED)
vals = {}
for gene in pcl.index:
    y = pcl.loc[gene, overlap].astype(float).values
    y = (y - y.mean()) / y.std(ddof=0)
    x = mat[overlap].astype(float)
    x = (x - x.mean(axis=1).values[:, None]) / x.std(axis=1, ddof=0).values[:, None]
    r = (x.values * y).mean(axis=1)
    probe = mat.index[int(np.nanargmax(np.abs(r)))]
    if r[np.nanargmax(np.abs(r))] < 0:
        mat.loc[probe] = -mat.loc[probe]
    vals[gene] = mat.loc[probe, all_samples].astype(float)

out = pd.DataFrame(vals).T
out.index.name = "gene"
out.to_csv(HERE / "GSE74123_19gene_matrix.csv")
print("saved GSE74123_19gene_matrix.csv", out.shape)
