# The Golden Target — reconciliation pipeline

Builds one golden record per drug-discovery target by reconciling five
overlapping source extracts (ChEMBL, UniProt, BindingDB, an internal target
registry, and a literature-mention table) against the EBI Proteins API
(https://www.ebi.ac.uk/proteins/api) as the authority, per the challenge
brief.

Supports Python 3.10, 3.11, and 3.12.

## Structure

```
.
├── pipeline.py            # entry point - the whole tool, stdlib only
├── goldentarget_config.json  # grader manifest: run_command, runtime_version, requirements_file
├── requirements.txt       # no third-party deps needed; Artifactory-ready if that changes
└── README.md
```

## Contract

```
python3 pipeline.py --data <pack_dir> --out submission.csv
```
`<pack_dir>` is expected to contain `source_chembl.csv`, `source_uniprot.csv`,
`source_bindingdb.csv`, `source_internal.csv`, and `source_publications.csv`.

Writes a single CSV to `--out` (default `submission.csv`) with one row per
golden record and one row per finding, distinguished by `record_type`:

| column | golden_record | finding |
|---|---|---|
| `record_type` | `golden_record` | `finding` |
| `gene` | ✓ | ✓ |
| `primary_accession` | ✓ | |
| `sources` | ✓ (`;`-joined) | |
| `observed` | | ✓ |
| `correct` | | ✓ |
| `retrieved_evidence` | | ✓ |
| `evidence_source` | | ✓ |
| `severity` | | ✓ |
| `classification` | | ✓ |

Unused columns are left blank on each row.

## Testing locally

```
cd r3-claude-olympics-golden-data
pip install -r requirements.txt
python3 pipeline.py --data ../path/to/data --out submission.csv
cat submission.csv
```

No `pip install` step is actually required — `pipeline.py` uses only the
Python standard library, so there's nothing to fetch from Artifactory or
PyPI. `requirements.txt` is kept for parity with the standard project layout.

The grader drives the tool via `goldentarget_config.json`'s `run_command`
(`python3 pipeline.py --data data/ --out submission.csv`), which expects the
hidden pack placed at `data/` in the repo root before running. To reproduce
that exact invocation locally, symlink or copy your pack to `data/` (already
gitignored) instead of passing an arbitrary `--data` path.

## Notes

- Every accession/gene/organism claim is verified live against the EBI
  Proteins API — nothing is hardcoded to this specific exam pack, so it
  should hold up against the hidden dataset with the same schema.
- `source_publications.csv`'s `pmid` values are never fetched (per the
  brief); literature-mention checks are resolved from the shipped
  `context_sentence` text only.
- A handful of "needs_review" findings (classification
  `ambiguous_mention_candidate`) are intentionally left as unconfirmed
  candidates rather than auto-corrected — those need a human/Claude to
  actually read the sentence before they're reported as a confirmed finding.
- Cold run against the real exam pack (~640 unique accessions): ~19s.
  Warm (cached) run: ~3s. Well inside the 5-minute budget.

