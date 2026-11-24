# -*- coding: utf-8 -*-
"""S07: independent direction check in GSE74123 (4 FECD vs 4 control).

Reads GSE74123_19gene_matrix.csv and writes TableS2_GSE74123_validation.csv.
Figure S4 is drawn separately by S07_plot_FigS4.py.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

HERE = Path(__file__).resolve().parent

DISCOVERY = {
    "OXSM": -0.762295699686256, "LRPPRC": -0.514086753831831,
    "IQGAP1": 0.689336167560357, "GYS1": -0.694114928734513,
    "FLNA": 0.932221804113881, "SLC7A11": 0.838902233047518,
    "FLNB": 0.716513389684057, "NCKAP1": 0.520498337691829,
    "RPN1": -0.317932209699466, "SLC3A2": -0.311520445752705,
    "ACTB": 0.445075108726993, "NUBPL": 0.311485245341571,
    "MYH10": 0.14840266408768, "ACTN4": 0.462479993546308,
    "CD2AP": -0.0903040119066736, "CAPZB": -0.00981607247612093,
    "MYH9": 0.229699173184198, "DSTN": 0.0939376494667901,
    "PDLIM1": -0.315655602061351,
}
CORE = ["OXSM", "GYS1", "FLNB", "IQGAP1"]
NINE = ["OXSM", "GYS1", "LRPPRC", "RPN1", "FLNA", "SLC7A11", "FLNB", "IQGAP1", "NCKAP1"]
CONTROL = ["GSM1910841", "GSM1910843", "GSM1910845", "GSM1910847"]
FECD = ["GSM1910840", "GSM1910842", "GSM1910844", "GSM1910846"]

mat = pd.read_csv(HERE / "GSE74123_19gene_matrix.csv", index_col=0)
rows = []
for gene in mat.columns:
    xc = mat.loc[CONTROL, gene].astype(float)
    xf = mat.loc[FECD, gene].astype(float)
    lfc = xf.mean() - xc.mean()
    welch = stats.ttest_ind(xf, xc, equal_var=False)
    mwu = stats.mannwhitneyu(xf, xc, alternative="two-sided")
    disc = DISCOVERY.get(gene, np.nan)
    rows.append([gene, disc, lfc, welch.pvalue, mwu.pvalue,
                 "concordant" if np.sign(lfc) == np.sign(disc) else "discordant"])
res = pd.DataFrame(rows, columns=["gene", "discovery_log2FC", "validation_log2FC",
                                  "welch_P", "mannwhitney_P", "direction"]).set_index("gene")
res.to_csv(HERE / "TableS2_GSE74123_validation.csv", float_format="%.4f")

for label, genes in (("core four", CORE), ("nine candidates", NINE), ("all 19", list(res.index))):
    print(f"{label}: {(res.loc[genes,'direction']=='concordant').sum()}/{len(genes)}")
r, p = stats.pearsonr(res["discovery_log2FC"], res["validation_log2FC"])
print(f"r = {r:.3f}, P = {p:.4f}")

print("saved TableS2_GSE74123_validation.csv")
