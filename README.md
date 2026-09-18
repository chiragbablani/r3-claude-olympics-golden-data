# The Golden Target — reconciliation tool

Builds one golden record per drug-discovery target by reconciling five
overlapping source extracts (ChEMBL, UniProt, BindingDB, an internal target
registry, and a literature-mention table) against the EBI Proteins API
(https://www.ebi.ac.uk/proteins/api) as the authority, per the challenge
brief.

## Structure

```
.
├── reconcile.py       # entry point - the whole tool, stdlib only
├── requirements.txt   # no third-party deps needed; Artifactory-ready if that changes
└── README.md
```

## Contract

```
python3 reconcile.py <pack_dir>
```
prints exactly one JSON object to stdout:
```json
{
  "unique_target_count": <int>,
  "golden_records": [{"gene": "...", "primary_accession": "...", "sources": ["..."]}],
  "findings": [
    {"gene": "...", "observed": "...", "correct": "...",
     "retrieved_evidence": "...", "evidence_source": "...",
     "severity": "...", "classification": "..."}
  ]
}
```
`<pack_dir>` is expected to contain `source_chembl.csv`, `source_uniprot.csv`,
`source_bindingdb.csv`, `source_internal.csv`, and `source_publications.csv`.

## Testing locally

```
python3 reconcile.py exam/ > output.json
cat output.json
```

No `pip install` step is required — `reconcile.py` uses only the Python
standard library, so there's nothing to fetch from Artifactory or PyPI.

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
