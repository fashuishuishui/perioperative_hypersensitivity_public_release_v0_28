# License and redistribution boundaries

This file is a scope notice, not a grant of rights for third-party source material.

## Included

- Original analysis code written for this project.
- Original build notes, manifests, figure-generation code, and aggregate derived outputs.
- Bounded drug aliases and project event-group mappings, subject to the source terms and any required author review before public release.

## Excluded

- Raw FAERS, Canada Vigilance, and JADER extracts, including report-level records and SQLite databases.
- Licensed MedDRA and MedDRA/J dictionaries, source tables, PT/LLT files, and SMQ membership files.
- Any local cache or intermediate database containing report-level records.
- Local absolute paths, credentials, tokens, private correspondence, and unpublished raw narratives.

## Official SMQ sensitivity

The repository may retain the SMQ identifier, version, algorithm, source-file hashes, aggregate membership counts, and aggregate estimates. It must not redistribute the licensed SMQ membership rows or MedDRA source files. The JADER branch must not be described as an official SMQ analysis without a licensed MedDRA/J PT-code crosswalk.

## License selection

Original code in this repository is released under the MIT License in `LICENSE`. Documentation and aggregate derived outputs are provided under the terms described in the repository documentation; this does not grant rights to redistribute regulatory source data or licensed MedDRA/MedDRA/J/SMQ content. Bounded drug aliases and project event-group mappings remain subject to the applicable source terms and should not be interpreted as a redistribution of the underlying proprietary terminology.
