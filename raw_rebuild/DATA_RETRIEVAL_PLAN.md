# v0.26 raw-source retrieval plan

## Freeze target

The first reproducibility target is the v0.25 study window, not the newest
available calendar quarter:

| Source | Intended source window | Required unit | Reproducibility target |
|---|---|---|---|
| FAERS | 2004Q1-2025Q4 | latest retained `PRIMARYID` per `CASEID` | Official quarterly ASCII ZIPs plus DELETE ledger |
| Canada Vigilance | historical extract dated 2025-07-31 | Canada report identifier | Official historical extract, or a documented same-schema contemporaneous replacement |
| JADER | monthly extract 2025-11 | latest report round per report identifier | Official PMDA 2025-11 extract |

The rebuilt data files stay outside the eventual public repository.  The
public release will contain only download manifests, checksums, source URLs,
build code, aggregate outputs, and permitted dictionaries.

## Retrieval gates

1. Capture the official landing page, source URL, access time, file name,
   byte count, and SHA-256 before parsing.
2. Validate archive integrity before extraction.  A source file with a partial
   suffix, failed checksum, or malformed required table is not accepted.
3. Preserve the original compressed file and a read-only extracted copy under
   `raw_inputs/`.  Processing writes only to `work/`.
4. Do not substitute an unofficial mirror, a third-party compiled SQLite file,
   or a renamed partial download for an official source file.
5. Stop and document the gap if an exact historical Canada or JADER release is
   unavailable from the agency.  A new release is a protocol amendment, not a
   silent replacement.

## FAERS retention rule

For every official ASCII quarter, parse `DEMO`, `DRUG`, `REAC`, and any
available `DELETE` member.  Delete-file identifiers are audited against both
`CASEID` and `PRIMARYID`; the identifier type is established empirically from
the source files before exclusion.  For nondeleted case families, retain one
report version by descending valid `FDA_DT`, then valid numeric `CASEVERSION`,
then numeric `PRIMARYID`.  A blank `CASEID` is retained under its `PRIMARYID`
and reported separately.  Only DRUG and REAC rows attached to retained IDs are
eligible for downstream flags.

## Completion criteria

The raw build can advance only when each source has a complete manifest,
source-specific retention audit, deterministic report-level flags, and a
numeric comparison against the frozen v0.25 branch.  The manuscript will be
revised after this comparison, not before.
