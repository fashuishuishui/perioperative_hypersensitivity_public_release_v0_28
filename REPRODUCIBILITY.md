# Reproducibility guide: v0.26 raw rebuild

## Analysis unit and source boundary

The analysis unit is a database-specific retained reporting unit. FAERS uses one retained report per case family after DELETE handling and latest-version selection. Canada Vigilance uses the supplied report identifier after a duplicate and linkage audit. JADER uses the supplied report identifier with a latest-round rule; the audited November 2025 local extract contained no repeated identifier.

The primary outputs use official re-downloaded FAERS archives, a frozen local Canada Vigilance extract dated 31 July 2025, and a frozen local JADER extract dated November 2025. The Canada and JADER local files are not redistributed. An official Open Canada archive through 30 November 2024 was separately downloaded and audited but was not substituted into the study because it would change the planned date window. PMDA terms must be accepted by the user before any official JADER archive can be downloaded.

## Required software

- Python 3.12.13
- SQLite 3.50.4
- R 4.5.1
- R package `DescTools` 0.99.60
- A licensed MedDRA v28.0 distribution only for the optional official *Anaphylactic reaction* SMQ sensitivity

## Required local inputs

1. FAERS official quarterly ASCII ZIP archives in `raw_inputs/faers_ascii_2004q1_2025q4/` or another path passed through the raw-rebuild script.
2. A Canada Vigilance source directory containing `reports.txt*`, `reactions.txt*`, and `report_drug.txt*`.
3. A JADER source directory containing `demo202511.csv*`, `reac202511.csv*`, and `drug202511.csv*`.
4. Licensed MedDRA files when the official-SMQ sensitivity is required.

The optional `*` suffix supports a local downloader suffix without copying, renaming, or redistributing source material.

## End-to-end runner

Run the following from the repository root after adapting paths:

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

Use `-ReuseValidatedFaersRawIntegrity` only after the named raw SQLite branch has passed the full source-build integrity check. The source manifest then records `raw_integrity_check_method` as `reused_from_source_build`; cache integrity is still checked afresh. Omit the switch when an independent full SQLite integrity rescan is required. Use `-SkipLicensedSmq` if licensed MedDRA files are unavailable. To reconstruct the FAERS raw SQLite branch from the source archives, add `-RebuildFaersRaw`; this is time- and storage-intensive and rewrites only the explicitly named source-build outputs.

The runner performs these stages in order:

1. audit or rebuild the FAERS raw source and retained-case cache;
2. rebuild FAERS, Canada Vigilance, and JADER report-level flags;
3. calculate the primary reference-setting, country, role, burden, broad-PT, teicoplanin-time, and JADER-structure diagnostics;
4. run the licensed official-SMQ sensitivity when requested;
5. build and validate the canonical artifacts;
6. independently reproduce Breslow-Day statistics in R;
7. regenerate figures and the aggregate-only supplementary package.

## Validation gate

The final check is `aggregate_outputs/schema_fixed/canonical_artifact_lock_v0_26.json`. It validates:

- report-level sum-backs and 2 x 2 tables;
- displayed shares and non-sparse RORs;
- country-stratum sum-backs;
- medication-burden and role-profile sum-backs;
- JADER structural-audit consistency;
- official-SMQ aggregate output consistency; and
- independent R/DescTools Breslow-Day calculations.

Do not copy a number into the manuscript, table, or figure unless it derives from a passing canonical lock. The public release must preserve the source manifests and aggregate result hashes, but must not include raw reports, report-level SQLite databases, or licensed MedDRA/SMQ terms.
