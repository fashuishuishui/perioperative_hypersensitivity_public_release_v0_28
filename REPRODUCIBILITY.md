# Reproduce this snapshot

## Without raw source data

From the repository root:

```console
python code/verify_aggregate_results.py .
python code/generate_figures.py
```

The verifier recomputes all 753 exported four-cell rows with standard-library formulas. Figure generation uses pandas and matplotlib and reads only the released tables. The generated figures describe the same reporting contrasts as the manuscript; native CSV fields retain full precision.

## Full source reconstruction

Use `README.txt` for the exact preparation steps and analysis order. Copy `config/paths.example.json` to the ignored `config/paths.json`, edit locations for your own permitted inputs, and start with an empty `work/` for a different source/dictionary lock. Do not reuse report-level caches after changing sources, aliases, role definitions or outcomes.

The FAERS initializer is `code/prepare_faers_sources.py`. It imports the 88 listed regulator archives, applies historical schema mappings, case-version/deletion handling and source-period-qualified joins. The named complete Canada/JADER extracts are rebuilt with their native identifiers/rounds; they are not recovered from the earlier 2025 caches. The official FAERS SMQ sensitivity requires local licensed MedDRA v28.0 membership and LLT/PT files. They are not supplied. JADER native-term diagnostics are not an official cross-language SMQ reconstruction.

All source-based analyses use the paths template; no report-level records are distributed. Analyses generate intermediate and aggregate output directories that are intentionally ignored by Git. The numbered-table exporter maps generated output names to S1–S17.

## Acceptance evidence

The released aggregate tables and approved manuscript figures are taken directly from the local submission candidate finalized on 18 October 2026. All 753 four-cell rows were independently recalculated; JADER burden strata sum to the corresponding primary tables, and fingerprint/burden checks are recorded in the local author audit. The public FAERS diagnostic runner was executed against the same source branch and matched 16 aggregate outputs; unchanged source extractions were reused. The raw-source import was not rerun solely to publish this code snapshot.

These checks concern numerical consistency and execution, not verification of individual clinical episodes. The authors remain responsible for final review, interpretation and manuscript content; publication authorization is not a substitute for that review.
