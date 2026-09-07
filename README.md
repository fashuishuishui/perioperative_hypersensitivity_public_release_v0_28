# Perioperative anaphylaxis pharmacovigilance audit

This repository contains the reproducibility materials for the manuscript *When suspected-drug attribution changes with the reference set: a multi-database pharmacovigilance audit of perioperative anaphylaxis reports*.

The study asks what an agent-level suspected-drug contrast in a spontaneous-reporting system describes when its reporting reference setting is made explicit. It is not a clinical risk, incidence, or culprit-adjudication analysis. NAP6 is used as clinical context, not as a validation dataset. The repository is an empirical pharmacovigilance audit package, not a clinical risk calculator or a general-purpose signal-detection framework.

## Current state

- FAERS 2004 Q1-2025 Q4 was rebuilt from 88 official quarterly ASCII archives. ZIP/CRC/SHA-256 checks, historical-schema mapping, DELETE handling, case-version retention, and source-period-safe joins are documented in `work/faers_raw_2004q1_2025q4_schema_fixed/`.
- Canada Vigilance and JADER were rerun from frozen local extracts after identifier and linkage audits. Their source-provenance boundary is explicitly recorded in the manuscript and supplement.
- The canonical aggregate output lock is `aggregate_outputs/schema_fixed/canonical_artifact_lock_v0_26.json` and reads `"result": "passed"` for the released analysis branch.
- This repository contains only code, bounded project mappings, aggregate outputs, figures, and manifests. It does not contain raw reports, SQLite databases, or licensed MedDRA/SMQ content.
- The archived release DOI will be recorded here after the public GitHub release is connected to Zenodo.

## Repository layout

- `code/`: data-build, analysis, validation, figure, and supplementary-package scripts.
- `raw_rebuild/`: FAERS archive retrieval and raw-source reconstruction scripts.
- `aggregate_outputs/schema_fixed/`: validated aggregate result tables.
- `figures/v0_26/`: manuscript figures in PDF, SVG, and PNG-preview formats.
- `supplementary/v0_26/`: aggregate-only supplementary tables and manifests.
- `run_v0_26_rebuild.ps1`: parameterized end-to-end runner.

## Reproduce the analysis

Obtain permitted source data before running. The runner requires local Canada Vigilance and JADER source directories. It also requires a licensed local MedDRA distribution for the optional official-SMQ sensitivity. The public repository intentionally does not include any of those source files.

```powershell
.\run_v0_26_rebuild.ps1 `
  -CanadaDir "<path-to-CANADA_20250731>" `
  -JaderDir "<path-to-JADER_202511>" `
  -PythonExe "<path-to-python>" `
  -RscriptExe "<path-to-Rscript>" `
  -FaersZipDir "<path-to-faers_ascii_2004q1_2025q4>" `
  -ReuseValidatedFaersRawIntegrity `
  -MedDraAscii "<path-to-licensed-MedDRA-MedAscii>"
```

Use `-ReuseValidatedFaersRawIntegrity` only when the named raw SQLite branch has already passed the full source-build integrity check; the rebuilt manifest records this audit path. Without this switch, reuse mode performs another full SQLite integrity scan before rebuilding the retained-report cache. Use `-SkipLicensedSmq` when a licensed MedDRA distribution is unavailable. This skips only the post hoc official-SMQ sensitivity; it does not alter the primary three-preferred-term analysis. Use `-RebuildFaersRaw` only when the local raw FAERS SQLite source branch should be rebuilt from the 88 audited ZIP archives.

See `REPRODUCIBILITY.md` and `LICENSE_BOUNDARIES.md` before release or rerun.
