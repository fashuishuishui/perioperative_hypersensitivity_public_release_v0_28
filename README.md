# Perioperative anaphylaxis: background-dependent reporting contrasts

**Current version: v1.0 — analytical lock 18 October 2026.**

This is the public code and aggregate-data companion to *Background dependence of within-class disproportionality in perioperative anaphylaxis reports: a three-database methodological case study*.
The repository URL retains its historical `_v0_28` suffix. The current default branch and `v1.0` tag contain the analysis accompanying the current four-drug manuscript, not the older v0.28 numerical snapshot. Earlier code and results remain retrievable through Git history and the pre-existing historical tag; they must not be mixed with this version.

## Use this snapshot

- `Table_S*.csv`: the complete, full-precision S1–S17 table groups, matching Online Resource 1 of the local submission candidate finalized on 18 October 2026.
- `code/` and `config/`: source-qualified analysis scripts, bounded drug/event dictionaries, public source inventories and a configurable input-path template.
- `figures/`: the three current manuscript figures (PDF/SVG and PNG previews).
- `README.txt`: input setup, analysis order, table conventions and interpretation limits.
- `ANALYSIS_MANIFEST.json`: source windows and the analysis software versions.
- `FILE_INDEX.csv`: inventory of the complete public snapshot.
- `REPRODUCIBILITY.md`: verification and figure-generation commands.

No raw regulatory reports, patient/report-level caches, licensed MedDRA dictionaries or complete SMQ membership, private input paths, manuscript/cover-letter files, or internal author-review archives are distributed.

## Study scope and main interpretation

Rocuronium, succinylcholine, cefazolin and teicoplanin are exploratory case studies of same-class reporting contrasts. FAERS uses 2004 Q1–2025 Q4 releases; Canada Vigilance covers through 31 May 2026; JADER uses the September 2026 release. Systems are analyzed separately, without pooling or cross-system patient linkage.

Medication-context FAERS RORs (95% CIs) are 1.38 (1.27–1.50), 0.89 (0.80–1.00), 3.92 (3.53–4.35) and 1.32 (1.00–1.73), respectively. RORs are reporting associations, not incidence or clinical culprit-risk rankings.

The monotonic rocuronium medication-burden gradient in FAERS was not reproduced in JADER. Teicoplanin fingerprint analyses cover only 36 of 66 target-core reports; 30 have incomplete required fields. Strict and PT-relaxed merging are scenario sensitivities, not bounds on unique patients or true duplicates. `Table_S4a_primary_four_cells.csv`, S5 and S15 provide the full counts and uncertainty behind these statements.

## Check the released aggregates

Use Python 3.12 and the dependencies listed in `requirements.txt` for analysis. The aggregate verifier uses only the Python standard library and does not require raw data, credentials or a MedDRA licence:

```console
python code/verify_aggregate_results.py .
```

Expected result: 753 four-cell rows recomputed, status `passed`. This checks RORs, log-ROR confidence intervals, Newcombe share-difference intervals, sparse-cell handling and Bonferroni intervals. It is an arithmetic check, not clinical adjudication or proof that no bias remains.

Full source reconstruction requires the named regulator releases and, for the official FAERS SMQ sensitivity, licensed MedDRA v28.0 files. See `README.txt` and `REPRODUCIBILITY.md`. Local input paths belong in the ignored `config/paths.json`, never in tracked files.

## Version and citation

Cite the `v1.0` tag or its commit for this analysis; `CITATION.cff` describes this snapshot. The DOI formerly present in this repository's metadata refers to an earlier archived snapshot and is not the DOI for v1.0. No new Zenodo deposit or DOI is asserted here. See `CHANGELOG.md` for the transition from the superseded snapshot.

Original code retains the existing MIT licence. Regulatory data and licensed terminology retain their own access and redistribution conditions; see `LICENSE_BOUNDARIES.md`.
