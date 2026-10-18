Online Resource 1 — aggregate tables and reproducible analysis code
Analytical lock: 18 October 2026

Tables are continuously grouped S1–S17. FILE_INDEX.csv identifies each file. Readable supporting methods and selected results are Online Resource 2 (ESM_2.pdf). Tables S4a and main Table 1 share one canonical numeric basis. Report counts are not persons, drug administrations, clinical incidence or adjudicated culprits.

Cell convention: a=agent-present core; b=agent-absent core; c=agent-present non-core; d=agent-absent non-core. ROR=ad/(bc). Both outcome arms use the same role/context/exclusion rules. Any zero cell or min(a,c)<5 withholds ROR without continuity correction. Share differences use Newcombe method 10. S17 uses 98.75% individual intervals for a four-agent Bonferroni sensitivity; other ratio intervals are descriptive 95% intervals. Fixed-country-mixture ORs are descriptive without CIs. Burden heterogeneity is explicitly Woolf, not Breslow–Day; MH is descriptive. Fingerprint matches are possible duplicate clusters, not verified unique patients. Missing/incomplete strict keys remain singletons. All expanded diagnostics are exploratory.

Input setup
Copy config/paths.example.json to config/paths.json and edit paths; paths can be relative to the configuration file or absolute on your own machine. Alternatively set PV_INPUT_CONFIG to a private configuration file. Use a separate output directory for a different source/dictionary lock rather than reusing cached work files. Prepare folders code/, config/, work/, aggregate_data/, metadata/ under one analysis root. No local paths, source report rows or complete licensed MedDRA membership are distributed.

FAERS preparation: obtain the 88 regulator ZIPs listed in config/faers_quarters.csv. prepare_faers_sources.py imports DEMO/DRUG/REAC with pre/post-2012 schema aliases and builds a source-qualified latest-case cache. Invoke it with --manifest config/faers_quarters.csv --zip-dir <your input ZIP directory> --output-dir <your prepared source directory>. Point raw_sqlite/latest_sqlite in paths.json to those outputs. ZIP CRC and required-member checks are performed. The existing source window used here ends 2025 Q4. config/faers_DEMO_inventory.csv lists the source DEMO members for demographic fields.

Canada/JADER preparation: supply complete reports.txt/reactions.txt/report_drug.txt for Canada (coverage through 2026-05-31) and demo202609.csv/drug202609.csv/reac202609.csv for JADER. JADER is CP932 CSV, Canada UTF-8 dollar-delimited. The source reader checks identifiers, latest rounds and linkage; it reconstructs receipt/report-quarter coverage. Different database source windows and native suspect constructs are not treated as equivalent.

Execution order (Python; dependencies and exact versions in ANALYSIS_MANIFEST.json)
1. code/build_drug_event_flags.py
2. code/run_targeted_diagnostics.py
3. code/rebuild_current_Canada_JADER.py
4. code/supporting_diagnostics.py
5. code/fingerprint_and_burden.py
6. code/fingerprint_followup.py
7. code/context_dictionary_sensitivity.py
8. code/burden_PS_feasibility.py
9. code/additional_composition_diagnostics.py
10. code/export_numbered_tables.py (copies analysis outputs to the manuscript table filenames).

Official FAERS SMQ sensitivity requires licensed MedDRA v28.0 pt.asc, llt.asc, smq_list.asc and smq_content.asc. Set their directories in the private input config. The supplied parser reads those local files; complete membership is not redistributed. Official SMQ sensitivity was not constructed in Japanese. The project core and native-term pharmacologic families are documented separately and are not called official SMQs.

The alias dictionary is bounded, source-specific lexical matching, not exhaustive WHO-DD/RxNorm normalization. S2c is a co-medication name audit, not an unmapped-exposure percentage. S8's expanded-marker populations are separate from the primary five-group definition. Adjacent anaphylaxis manifestations are not excluded. Pharmacologic exclusions apply only to NMBAs.

To verify tables, run code/verify_aggregate_results.py against the distributed table directory. Reproduction from public sources requires the named releases and a MedDRA licence; the supplied aggregate tables allow four-cell arithmetic checks without source report data. The private work directory must never be publicly uploaded.

Additional composition diagnostics
S5f/S5g contain JADER burden strata and descriptive MH/Woolf summaries. Count distinct native drug-sequence values only at the retained report round; this is not a standardized active-ingredient count shared with FAERS. All strata sum to the JADER primary cells. A monotonic FAERS rocuronium pattern does not reproduce in JADER.
S15g/S15h contain PT-relaxed FAERS merging and completeness by outcome arm. Keep the same complete demographic/date/country fields and core status, but omit PT identities in BOTH arms. Incomplete fingerprints remain singletons. As with strict merging, retain actual labels from the latest FDA_DT/numeric PRIMARYID representative. These rules are scenarios, not bounds on true duplicates or patients.
S8 reports all four agents under two expanded dictionary definitions. The primary dictionary is not preregistered or validated as identifying surgery. Expanded brand coverage and broader pharmacologic coverage are distinct sensitivities.
