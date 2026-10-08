# CBIO055 results (results/results.json, src/run.py, results/run.log). Poore et al. TCGA Kraken Voom-SNM, 17,625 samples, 1,287 genera; 5,492 early (stage I/II) primary tumors, 1,206 solid-tissue normals. Grouped (case) 5-fold CV, L2 logistic regression.
- D2 pan-cancer early tumor vs solid normal (types that have normals; 5,245 tumors, 1,206 normals): microbial AUROC 0.863; metadata baseline (disease+center+strategy) 0.696; center+strategy only 0.612. Gate (>=0.80 and >= baseline+0.05) MET.
- D3 within-type CV, 14 cancer types: microbial AUROC 0.845-0.992 (see json); center+strategy-only baseline 0.29-0.70 (near or below chance).
- D4 leave-one-cancer-type-out mean AUROC 0.79 (range in json; lowest thyroid 0.637, highest rectum 0.92). Gate (>=0.70) MET.
- D5 train on late-stage tumor + half of normal cases, test on early tumor + other normal cases: 0.857 (metadata 0.627).
- D6 21-class cancer type among early tumors: macro-F1 0.583 vs center+strategy-only 0.093 (majority class 16.8%).
- D7 label shuffle control AUROC 0.498 (no pipeline leakage).
- D9 top-20 genera overlap across folds, mean Jaccard 0.44; 9 genera in all folds (listed in json, mostly odd taxa such as Anaplasma, Gemmata, Legionella; no biological claims made).
- D10 vs blood-derived normal: microbial 0.906 but metadata alone 0.968 - confounded by sample type, not evidence of detection.
Limits: this is TCGA tissue, not plasma or healthy donors. Tumor and normal differ by tissue handling; Voom-SNM already removed plate/center effects, and the metadata baselines do not capture all batch structure. The Poore dataset is contested for contamination (Gihawi 2023), and the distinguishing genera look contamination-like. Gates being met does not show a clinical early-detection signal. Contamination filtering not run. Replication/audit, no novelty claim.

## Audit (second pass) - D2 gate DOWNGRADED to "met against a weak baseline, not robust"
Prereg AUDIT_PREREG.md committed before running (one column-name fix noted there; Q1/Q2 unaffected). src/audit.py, results/audit.json, results/audit.log. Data SHA256 re-verified; D2 reproduced exactly (0.8631).
- Q1 leave-one-center-out (7 centers with enough samples): Baylor 0.733, Broad 0.672, BC Genome Sciences 0.888, Harvard 0.556, MD Anderson 0.557, UNC 0.603 (4,734 samples, 73% of the data), WashU 0.669. Mean 0.668, below the pre-specified 0.70 bar: MISSED. The signal transfers poorly to centers not seen in training (UNC, the largest, 0.60).
- Q2 within (disease x center x strategy) strata (22 strata with >=15 per class): case-weighted mean AUROC 0.892. Bar (>=0.70) MET: signal is not only between-batch.
- Q3 extended metadata baseline (adds tissue source site, FFPE, analyte, platform): 0.809 (original baseline 0.696). Microbial 0.863 vs 0.809: the gate (>= baseline + 0.05, i.e. 0.859) is met by a margin of 0.004: fragile, not robust.
Revised reading: within-batch separation is real (Q2) but a richer batch-aware baseline recovers most of the AUROC (Q3) and cross-center transfer is weak (Q1). D2's "gate MET" now reads "met against a weak baseline, not robust to center hold-out"; D4 and others unchanged. Contamination remains untested. Still not a clinical early-detection signal.
