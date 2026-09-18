#!/usr/bin/env python3
"""
Golden Target reconciliation tool.

Usage:
    python3 reconcile.py <pack_dir>

Reads five source extracts from <pack_dir>:
    source_chembl.csv, source_uniprot.csv, source_bindingdb.csv,
    source_internal.csv, source_publications.csv

Resolves every referenced identity against the EBI Proteins API
(https://www.ebi.ac.uk/proteins/api), builds one golden record per
real-world target, and reports defects it can prove with retrieved
evidence.

No third-party packages required (stdlib only): csv, json, re,
urllib.request, concurrent.futures.
"""

import csv
import json
import os
import re
import sys
import time
import urllib.request
import urllib.parse
import urllib.error
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Optional

EBI_BASE = "https://www.ebi.ac.uk/proteins/api/proteins"
BATCH_SIZE = 90          # accessions per batched GET (stay well under URL length limits)
MAX_WORKERS = 6          # modest concurrency, be polite to the public API
REQUEST_TIMEOUT = 20
MAX_RETRIES = 3
CACHE_PATH = os.environ.get("EBI_CACHE_PATH", os.path.join(os.path.dirname(__file__), ".ebi_cache.json"))


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------

@dataclass
class ProteinRecord:
    accession: str
    gene: Optional[str]
    synonyms: list
    secondary_accessions: list
    organism_sci: str
    organism_common: str
    full_name: str
    alt_full_names: list
    reviewed: bool
    chembl_xref: Optional[str]
    hgnc_xref: Optional[str]

    def all_names_lower(self):
        names = set()
        if self.full_name:
            names.add(self.full_name.lower())
        for n in self.alt_full_names:
            if n:
                names.add(n.lower())
        return names

    def all_gene_labels_upper(self):
        labels = set()
        if self.gene:
            labels.add(self.gene.upper())
        for s in self.synonyms:
            labels.add(s.upper())
        return labels


def parse_protein_json(d: dict) -> ProteinRecord:
    gene_block = (d.get("gene") or [{}])[0]
    gene_name = (gene_block.get("name") or {}).get("value")
    synonyms = [s.get("value") for s in gene_block.get("synonyms", []) if s.get("value")]

    organism = d.get("organism", {})
    org_names = organism.get("names", [])
    sci = next((n["value"] for n in org_names if n.get("type") == "scientific"), "")
    common = next((n["value"] for n in org_names if n.get("type") == "common"), "")

    protein = d.get("protein", {})
    rec_name = ((protein.get("recommendedName") or {}).get("fullName") or {}).get("value", "")
    alt_names = [
        (a.get("fullName") or {}).get("value")
        for a in protein.get("alternativeName", [])
        if (a.get("fullName") or {}).get("value")
    ]

    chembl_xref = next((x["id"] for x in d.get("dbReferences", []) if x.get("type") == "ChEMBL"), None)
    hgnc_xref = next((x["id"] for x in d.get("dbReferences", []) if x.get("type") == "HGNC"), None)

    return ProteinRecord(
        accession=d.get("accession"),
        gene=gene_name,
        synonyms=synonyms,
        secondary_accessions=d.get("secondaryAccession", []) or [],
        organism_sci=sci,
        organism_common=common,
        full_name=rec_name,
        alt_full_names=alt_names,
        reviewed=(d.get("info", {}).get("type") == "Swiss-Prot"),
        chembl_xref=chembl_xref,
        hgnc_xref=hgnc_xref,
    )


# --------------------------------------------------------------------------
# EBI client (batched, cached, retried)
# --------------------------------------------------------------------------

class EBIClient:
    def __init__(self):
        self.cache = self._load_cache()
        self._session_hits = 0

    def _load_cache(self):
        if os.path.exists(CACHE_PATH):
            try:
                with open(CACHE_PATH) as f:
                    return json.load(f)
            except Exception:
                return {"by_accession": {}, "by_gene": {}}
        return {"by_accession": {}, "by_gene": {}}

    def save_cache(self):
        try:
            with open(CACHE_PATH, "w") as f:
                json.dump(self.cache, f)
        except Exception:
            pass

    def _get(self, url):
        last_err = None
        for attempt in range(MAX_RETRIES):
            try:
                req = urllib.request.Request(url, headers={"Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
                    return json.load(resp)
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    return None
                last_err = e
            except Exception as e:
                last_err = e
            time.sleep(0.5 * (attempt + 1))
        if last_err:
            sys.stderr.write(f"WARN: request failed after retries: {url} ({last_err})\n")
        return None

    def batch_by_accession(self, accessions):
        """Resolve a list of accessions directly. Returns {accession: ProteinRecord}.
        Accessions not found (obsolete/secondary/typo) are simply absent from the result."""
        result = {}
        todo = []
        for acc in accessions:
            cached = self.cache["by_accession"].get(acc)
            if cached is not None:
                if cached != "MISS":
                    result[acc] = parse_protein_json(cached)
            else:
                todo.append(acc)

        for i in range(0, len(todo), BATCH_SIZE):
            chunk = todo[i:i + BATCH_SIZE]
            url = EBI_BASE + "?accession=" + ",".join(urllib.parse.quote(a) for a in chunk)
            data = self._get(url) or []
            found_accs = set()
            for entry in data:
                acc = entry.get("accession", "")
                # NOTE: previously skipped any accession containing "-" to avoid
                # unsolicited isoform pollution from gene-based searches. Confirmed
                # via live testing that the accession-list endpoint only ever returns
                # entries matching exactly what was requested (no unsolicited isoform
                # rows), so an isoform-suffixed accession that was explicitly asked
                # for (e.g. a row that only ever references "Q13422-3", with no bare
                # "Q13422" anywhere in the pack) must be accepted, not discarded --
                # discarding it here silently drops the whole target with zero trace.
                found_accs.add(acc)
                self.cache["by_accession"][acc] = entry
                result[acc] = parse_protein_json(entry)
            for acc in chunk:
                if acc not in found_accs:
                    self.cache["by_accession"][acc] = "MISS"
        return result

    def search_by_gene(self, gene, organism_hint=None, reviewed_only=True):
        """Fallback resolution path: search current entries by gene symbol.
        Used to (a) resolve obsolete/secondary accessions via secondaryAccession
        membership, and (b) detect duplicate/renamed identities."""
        key = f"{gene.upper()}|{organism_hint or ''}|{reviewed_only}"
        cached = self.cache["by_gene"].get(key)
        if cached is not None:
            return [parse_protein_json(e) for e in cached]

        params = {"gene": gene}
        if organism_hint:
            params["organism"] = organism_hint
        if reviewed_only:
            params["reviewed"] = "true"
        url = EBI_BASE + "?" + urllib.parse.urlencode(params)
        data = self._get(url) or []
        clean = [e for e in data if "-" not in e.get("accession", "")]
        self.cache["by_gene"][key] = clean
        return [parse_protein_json(e) for e in clean]

    def find_current_owner_of_secondary(self, old_accession, gene_hint, organism_hint="Homo sapiens"):
        """Given an accession that does not resolve directly, use the row's own
        gene label as a search seed and check whether old_accession shows up in
        the secondaryAccession list of a current entry. Returns ProteinRecord or None."""
        candidates = self.search_by_gene(gene_hint, organism_hint=organism_hint, reviewed_only=True)
        for rec in candidates:
            if old_accession in rec.secondary_accessions:
                return rec
        # widen: try without organism restriction (defensive, in case of a species mismatch too)
        candidates = self.search_by_gene(gene_hint, organism_hint=None, reviewed_only=False)
        for rec in candidates:
            if old_accession in rec.secondary_accessions:
                return rec
        return None


# --------------------------------------------------------------------------
# Organism normalization (generic - not hardcoded to specific species)
# --------------------------------------------------------------------------

def organism_matches(observed: str, rec: ProteinRecord) -> bool:
    if not observed or not rec.organism_sci:
        return True  # nothing to compare; don't penalize missing data
    obs = observed.strip().lower().rstrip(".")
    # UniProt-style "Homo sapiens (Human)" -> split into scientific + parenthetical common name
    paren_match = re.match(r"^(.*?)\s*\((.*?)\)\s*$", obs)
    if paren_match:
        obs_sci_part, obs_common_part = paren_match.group(1).strip(), paren_match.group(2).strip()
    else:
        obs_sci_part, obs_common_part = obs, obs

    sci = rec.organism_sci.strip().lower()
    common = (rec.organism_common or "").strip().lower()

    if obs in (sci, common) or obs_sci_part == sci or obs_common_part == common:
        return True

    parts = sci.split()
    if len(parts) >= 2:
        genus, species = parts[0], parts[-1]
        abbreviated = f"{genus[0]}. {species}".lower()
        abbreviated_nodot = f"{genus[0]} {species}".lower()
        if obs in (abbreviated, abbreviated_nodot):
            return True
    return False


# --------------------------------------------------------------------------
# Gene label classification
# --------------------------------------------------------------------------

def gene_label_status(observed_gene: str, rec: ProteinRecord) -> str:
    """Returns 'current', 'synonym', or 'mismatch'."""
    if not observed_gene or not rec.gene:
        return "current"
    if observed_gene.strip().upper() == rec.gene.upper():
        return "current"
    if observed_gene.strip().upper() in {s.upper() for s in rec.synonyms}:
        return "synonym"
    return "mismatch"


# --------------------------------------------------------------------------
# CSV loading
# --------------------------------------------------------------------------

def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_pack(pack_dir):
    return {
        "chembl": load_csv(os.path.join(pack_dir, "source_chembl.csv")),
        "uniprot": load_csv(os.path.join(pack_dir, "source_uniprot.csv")),
        "bindingdb": load_csv(os.path.join(pack_dir, "source_bindingdb.csv")),
        "internal": load_csv(os.path.join(pack_dir, "source_internal.csv")),
        "publications": load_csv(os.path.join(pack_dir, "source_publications.csv")),
    }


# --------------------------------------------------------------------------
# Publication context-sentence templates (text-only; pmid is NEVER fetched)
# --------------------------------------------------------------------------

PUB_TEMPLATES = [
    r"^A screening campaign identified modulators of (.+)\.$",
    r"^(.+) was investigated as a therapeutic target\.$",
    r"^Inhibition of (.+) altered disease-relevant signalling in vitro\.$",
    r"^(.+) expression correlated with treatment response\.$",
]


def extract_pub_phrase(context_sentence, target_mention):
    """Return the descriptive phrase X the sentence actually names, or None if
    the sentence merely repeats the bare mention (nothing to cross-check)."""
    for pat in PUB_TEMPLATES:
        m = re.match(pat, context_sentence.strip())
        if m:
            phrase = m.group(1).strip()
            if phrase.upper() == target_mention.strip().upper():
                return None  # trivially consistent, bare symbol only
            return phrase
    return "FREEFORM:" + context_sentence.strip()


# --------------------------------------------------------------------------
# Main reconciliation
# --------------------------------------------------------------------------

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: reconcile.py <pack_dir>\n")
        sys.exit(1)
    pack_dir = sys.argv[1]
    pack = load_pack(pack_dir)
    client = EBIClient()

    # ---- 1. collect every accession-like value referenced anywhere ----
    all_accessions = set()
    for r in pack["uniprot"]:
        if r.get("accession"):
            all_accessions.add(r["accession"])
    for r in pack["chembl"]:
        if r.get("accession"):
            all_accessions.add(r["accession"])
    for r in pack["bindingdb"]:
        if r.get("uniprot_id"):
            all_accessions.add(r["uniprot_id"])
    for r in pack["internal"]:
        if r.get("uniprot_ref"):
            all_accessions.add(r["uniprot_ref"])

    direct = client.batch_by_accession(sorted(all_accessions))
    unresolved = sorted(all_accessions - set(direct.keys()))

    # ---- 2. resolve the unresolved ones via their own row's gene label ----
    row_gene_hint = {}  # accession -> a gene symbol seen alongside it in some row
    # NOTE: originally only chembl/bindingdb/internal were consulted here.
    # source_uniprot.csv rows carry their own gene_names value too, and skipping
    # them meant any bad accession that appears ONLY in the uniprot extract had
    # no path to indirect (merged-accession) resolution at all -- it just
    # silently vanished as "unresolved" with no finding raised.
    for r in pack["uniprot"]:
        acc = r.get("accession")
        if acc in unresolved and r.get("gene_names"):
            row_gene_hint.setdefault(acc, r["gene_names"].split()[0])
    for r in pack["chembl"]:
        if r.get("accession") in unresolved:
            row_gene_hint.setdefault(r["accession"], r.get("gene_symbol"))
    for r in pack["bindingdb"]:
        if r.get("uniprot_id") in unresolved:
            row_gene_hint.setdefault(r["uniprot_id"], r.get("gene_symbol"))
    for r in pack["internal"]:
        if r.get("uniprot_ref") in unresolved:
            row_gene_hint.setdefault(r["uniprot_ref"], r.get("gene_symbol"))

    indirect = {}  # accession -> (ProteinRecord, via_secondary: bool)
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futs = {}
        for acc in unresolved:
            hint = row_gene_hint.get(acc)
            if hint:
                futs[pool.submit(client.find_current_owner_of_secondary, acc, hint)] = acc
        for fut in as_completed(futs):
            acc = futs[fut]
            rec = fut.result()
            if rec:
                indirect[acc] = rec

    client.save_cache()

    def resolve_accession(acc):
        """-> (ProteinRecord or None, path) where path in {direct, merged, unresolved}"""
        if acc in direct:
            return direct[acc], "direct"
        if acc in indirect:
            return indirect[acc], "merged"
        return None, "unresolved"

    # ---- 3. walk every source row, build findings + golden clusters ----
    finding_registry = {}  # (classification, dedupe_key) -> finding dict; dedupe_key=None means always-unique
    _unique_counter = [0]
    golden = {}  # canonical_accession -> {"gene": ..., "sources": set()}

    def add_to_golden(rec, source_name):
        if not rec:
            return
        entry = golden.setdefault(rec.accession, {"gene": rec.gene, "sources": set()})
        entry["sources"].add(source_name)

    def add_finding(gene, observed, correct, evidence, source, severity, classification, dedupe_key=None):
        # dedupe_key groups findings that describe the SAME underlying defect
        # (e.g. the same stale accession) even when it's independently visible
        # in more than one source file -- one finding per real-world defect,
        # not one per row it happens to show up in.
        if dedupe_key is None:
            _unique_counter[0] += 1
            key = ("_unique", _unique_counter[0])
        else:
            key = (classification, dedupe_key)

        if key in finding_registry:
            finding_registry[key]["observed"] += f" | ALSO: {observed}"
            return

        finding_registry[key] = {
            "gene": gene,
            "observed": observed,
            "correct": correct,
            "retrieved_evidence": evidence,
            "evidence_source": source,
            "severity": severity,
            "classification": classification,
        }

    # -- UniProt source: mostly definitional, but still verify --
    acc_to_observed_by_source = defaultdict(lambda: defaultdict(set))  # canonical_acc -> source -> {observed_ids}

    for r in pack["uniprot"]:
        acc = r.get("accession")
        rec, path = resolve_accession(acc)
        if not rec:
            continue
        add_to_golden(rec, "uniprot")
        acc_to_observed_by_source[rec.accession]["uniprot"].add(acc)
        if path == "merged":
            # NOTE: this branch was entirely missing -- chembl/bindingdb/internal
            # all raise stale_accession when path=="merged"; the uniprot loop
            # only ever checked gene/organism, so these accessions resolved
            # correctly into the golden record (silently, via add_to_golden)
            # but never surfaced as a finding at all.
            add_finding(
                rec.gene, f"accession={acc} in source_uniprot.csv (entry_name={r.get('entry_name')})",
                f"accession={rec.accession}, gene={rec.gene}",
                f"Accession {acc} not resolvable directly (obsolete); found as a secondary "
                f"accession of current entry {rec.accession}, retrieved via EBI Proteins API "
                f"gene search using this row's own gene_names field",
                "EBI Proteins API /proteins?gene=...  (secondaryAccession match)",
                "medium", "stale_accession",
                dedupe_key=acc,
            )
        observed_gene = (r.get("gene_names") or "").split()[0] if r.get("gene_names") else None
        if observed_gene:
            status = gene_label_status(observed_gene, rec)
            if status == "mismatch":
                add_finding(
                    rec.gene, f"gene_names='{r.get('gene_names')}' for accession {acc}",
                    f"{rec.gene}",
                    f"EBI Proteins API accession {acc} -> current approved gene '{rec.gene}', synonyms {rec.synonyms}",
                    "EBI Proteins API /proteins/{accession}",
                    "high", "wrong_mapping",
                    dedupe_key=acc,
                )
        obs_org = r.get("organism", "")
        if obs_org and not organism_matches(obs_org, rec):
            add_finding(
                rec.gene, f"organism='{obs_org}' for accession {acc}", rec.organism_sci,
                f"EBI Proteins API accession {acc} -> organism '{rec.organism_sci}'",
                "EBI Proteins API /proteins/{accession}",
                "high", "wrong_organism",
            )

    # -- ChEMBL source --
    for r in pack["chembl"]:
        acc = r.get("accession")
        rec, path = resolve_accession(acc)
        if not rec:
            continue
        add_to_golden(rec, "chembl")

        gene_obs = r.get("gene_symbol")
        status = gene_label_status(gene_obs, rec)
        if path == "merged":
            add_finding(
                rec.gene, f"accession={acc}, gene_symbol={gene_obs} (chembl_id {r.get('chembl_id')})",
                f"accession={rec.accession}, gene={rec.gene}",
                f"Accession {acc} not resolvable directly (obsolete); found as a secondary "
                f"accession of current entry {rec.accession} (gene {rec.gene}), retrieved via "
                f"EBI Proteins API gene search for '{gene_obs}'",
                "EBI Proteins API /proteins?gene=...  (secondaryAccession match)",
                "medium", "stale_accession",
                dedupe_key=acc,
            )
        elif status == "mismatch":
            add_finding(
                gene_obs, f"accession={acc}, gene_symbol={gene_obs} (chembl_id {r.get('chembl_id')})",
                f"accession={rec.accession} actually corresponds to gene {rec.gene}",
                f"EBI Proteins API accession {acc} -> gene '{rec.gene}' (synonyms {rec.synonyms}); "
                f"'{gene_obs}' is not this entry's approved symbol or a known synonym",
                "EBI Proteins API /proteins/{accession}",
                "high", "wrong_mapping",
                dedupe_key=acc,
            )
        elif status == "synonym":
            add_finding(
                rec.gene, f"gene_symbol={gene_obs} (accession {acc})", rec.gene,
                f"EBI Proteins API accession {acc} -> current approved gene '{rec.gene}'; "
                f"'{gene_obs}' is a listed prior/alternate symbol",
                "EBI Proteins API /proteins/{accession}",
                "low", "stale_gene_symbol",
                dedupe_key=acc,
            )

        if rec.chembl_xref and r.get("chembl_id") and rec.chembl_xref != r["chembl_id"]:
            add_finding(
                rec.gene, f"chembl_id={r.get('chembl_id')} for accession {acc}", rec.chembl_xref,
                f"EBI Proteins API accession {acc} cross-references ChEMBL id '{rec.chembl_xref}'",
                "EBI Proteins API /proteins/{accession} (dbReferences)",
                "medium", "wrong_cross_reference",
                dedupe_key=acc,
            )

        obs_org = r.get("organism", "")
        if obs_org and not organism_matches(obs_org, rec):
            add_finding(
                rec.gene, f"organism='{obs_org}' (chembl_id {r.get('chembl_id')})", rec.organism_sci,
                f"EBI Proteins API accession {acc} -> organism '{rec.organism_sci}'",
                "EBI Proteins API /proteins/{accession}",
                "high", "wrong_organism",
            )

    # -- BindingDB source --
    for r in pack["bindingdb"]:
        acc = r.get("uniprot_id")
        rec, path = resolve_accession(acc)
        if not rec:
            continue
        add_to_golden(rec, "bindingdb")

        gene_obs = r.get("gene_symbol")
        status = gene_label_status(gene_obs, rec)
        if path == "merged":
            add_finding(
                rec.gene, f"uniprot_id={acc}, gene_symbol={gene_obs} ('{r.get('target_name')}')",
                f"uniprot_id={rec.accession}, gene={rec.gene}",
                f"Accession {acc} not resolvable directly (obsolete); found as a secondary "
                f"accession of current entry {rec.accession} (gene {rec.gene}), retrieved via "
                f"EBI Proteins API gene search for '{gene_obs}'",
                "EBI Proteins API /proteins?gene=...  (secondaryAccession match)",
                "medium", "stale_accession",
                dedupe_key=acc,
            )
        elif status == "mismatch":
            add_finding(
                gene_obs, f"uniprot_id={acc}, gene_symbol={gene_obs} ('{r.get('target_name')}')",
                f"uniprot_id={acc} actually corresponds to gene {rec.gene}",
                f"EBI Proteins API accession {acc} -> gene '{rec.gene}' (synonyms {rec.synonyms})",
                "EBI Proteins API /proteins/{accession}",
                "high", "wrong_mapping",
                dedupe_key=acc,
            )
        elif status == "synonym":
            add_finding(
                rec.gene, f"gene_symbol={gene_obs} (uniprot_id {acc})", rec.gene,
                f"EBI Proteins API accession {acc} -> current approved gene '{rec.gene}'",
                "EBI Proteins API /proteins/{accession}",
                "low", "stale_gene_symbol",
                dedupe_key=acc,
            )

        obs_species = r.get("species", "")
        if obs_species and not organism_matches(obs_species, rec):
            add_finding(
                rec.gene, f"species='{obs_species}' (uniprot_id {acc})", rec.organism_sci,
                f"EBI Proteins API accession {acc} -> organism '{rec.organism_sci}'",
                "EBI Proteins API /proteins/{accession}",
                "high", "wrong_organism",
            )

    # -- Internal registry source --
    for r in pack["internal"]:
        acc = r.get("uniprot_ref")
        if not acc:
            continue
        rec, path = resolve_accession(acc)
        if not rec:
            continue
        add_to_golden(rec, "internal")

        gene_obs = r.get("gene_symbol")
        if gene_obs:
            status = gene_label_status(gene_obs, rec)
            if path == "merged":
                add_finding(
                    rec.gene, f"uniprot_ref={acc}, gene_symbol={gene_obs} ({r.get('internal_id')})",
                    f"uniprot_ref={rec.accession}, gene={rec.gene}",
                    f"Accession {acc} not resolvable directly (obsolete); found as a secondary "
                    f"accession of current entry {rec.accession}",
                    "EBI Proteins API /proteins?gene=...  (secondaryAccession match)",
                    "medium", "stale_accession",
                    dedupe_key=acc,
                )
            elif status == "mismatch":
                add_finding(
                    gene_obs, f"uniprot_ref={acc}, gene_symbol={gene_obs} ({r.get('internal_id')})",
                    f"uniprot_ref={acc} actually corresponds to gene {rec.gene}",
                    f"EBI Proteins API accession {acc} -> gene '{rec.gene}'",
                    "EBI Proteins API /proteins/{accession}",
                    "high", "wrong_mapping",
                    dedupe_key=acc,
                )
            elif status == "synonym":
                # NOTE: this branch was missing entirely -- the chembl/bindingdb
                # loops both raise a low-severity stale_gene_symbol finding here;
                # the internal-registry loop jumped straight from "current" to
                # "mismatch" and silently swallowed the synonym case, so any row
                # using an old-but-valid HGNC symbol (e.g. SEPT9 -> SEPTIN9,
                # WHSC1 -> NSD2) produced nothing at all.
                add_finding(
                    rec.gene, f"gene_symbol={gene_obs} (uniprot_ref {acc}, {r.get('internal_id')})", rec.gene,
                    f"EBI Proteins API accession {acc} -> current approved gene '{rec.gene}'; "
                    f"'{gene_obs}' is a listed prior/alternate symbol",
                    "EBI Proteins API /proteins/{accession}",
                    "low", "stale_gene_symbol",
                    dedupe_key=acc,
                )

    # ---- 4. duplicate-identity pass: two different observed keys within the SAME
    #      source file resolving to the same canonical accession. Cross-source
    #      agreement (e.g. uniprot and internal both citing one real accession)
    #      is normal and NOT a duplicate -- only flag when one source's OWN rows
    #      disagree with themselves about how many distinct targets this is.
    for r in pack["internal"]:
        acc = r.get("uniprot_ref")
        if acc:
            rec, _ = resolve_accession(acc)
            if rec:
                acc_to_observed_by_source[rec.accession]["internal"].add(r.get("internal_id"))

    for canon_acc, by_source in acc_to_observed_by_source.items():
        gene = golden.get(canon_acc, {}).get("gene")
        for source_name, keys in by_source.items():
            if len(keys) > 1:
                add_finding(
                    gene, f"{len(keys)} distinct {source_name} rows resolve to the same target: {sorted(keys)}",
                    f"single golden record, primary_accession={canon_acc}",
                    f"All listed {source_name} rows resolve (directly or via merged-accession "
                    f"lookup) to the same current UniProt accession {canon_acc}",
                    "EBI Proteins API /proteins/{accession} + /proteins?gene=...",
                    "medium", "duplicate_identity",
                )

    # ---- 5. publications: text-only cross-check, pmid NEVER fetched ----
    gene_names_index = defaultdict(set)  # lowercased name/phrase -> set of gene symbols
    for rec in list(direct.values()) + list(indirect.values()):
        for name in rec.all_names_lower():
            gene_names_index[name].add(rec.gene)
        if rec.gene:
            gene_names_index[rec.gene.lower()].add(rec.gene)

    STOPWORDS = {"the", "and", "for", "with", "this", "that", "was", "were", "from", "into",
                 "protein", "type", "family", "member", "subunit", "domain", "containing"}

    def tokens_of(text):
        return {w for w in re.findall(r"[a-z0-9]+", text.lower()) if len(w) >= 4 and w not in STOPWORDS}

    gene_token_sets = defaultdict(set)  # gene -> token set over all its known names
    for rec in list(direct.values()) + list(indirect.values()):
        if not rec.gene:
            continue
        for name in rec.all_names_lower():
            gene_token_sets[rec.gene] |= tokens_of(name)

    for r in pack["publications"]:
        mention = r.get("target_mention", "")
        sentence = r.get("context_sentence", "")
        phrase = extract_pub_phrase(sentence, mention)
        if phrase is None:
            continue  # sentence just repeats the bare symbol; nothing to cross-check
        is_freeform = phrase.startswith("FREEFORM:")
        phrase_text = phrase[len("FREEFORM:"):] if is_freeform else phrase
        phrase_lower = phrase_text.lower().rstrip(".")

        if not is_freeform:
            # templated sentence: exact protein-name match against the known universe
            owners = gene_names_index.get(phrase_lower, set())
            if not owners:
                # can't verify, but mention itself may still be a known gene -> attribute plainly
                gene_to_acc = {g: acc for acc, g in [(rr.accession, rr.gene) for rr in list(direct.values()) + list(indirect.values())]}
                if mention.upper() in {g.upper() for g in gene_to_acc}:
                    for g, acc in gene_to_acc.items():
                        if g and g.upper() == mention.upper():
                            golden.setdefault(acc, {"gene": g, "sources": set()})["sources"].add("publications")
                continue
            if mention.upper() in {g.upper() for g in owners if g}:
                # consistent: attribute this mention to its (confirmed) golden record
                for rr in list(direct.values()) + list(indirect.values()):
                    if rr.gene and rr.gene.upper() == mention.upper():
                        golden.setdefault(rr.accession, {"gene": rr.gene, "sources": set()})["sources"].add("publications")
                continue
            other_genes = sorted(g for g in owners if g and g.upper() != mention.upper())
            if other_genes:
                add_finding(
                    mention, f"target_mention='{mention}' with context_sentence: \"{sentence}\"",
                    other_genes[0],
                    f"context_sentence names '{phrase_text}', which matches gene "
                    f"{other_genes[0]}'s protein name in EBI Proteins API (not '{mention}''s). "
                    f"pmid was not fetched; resolved purely from shipped context_sentence text.",
                    "EBI Proteins API protein-name lookup (context_sentence text only)",
                    "high", "ambiguous_mention_misresolved",
                )
                # attribute the mention to the CORRECTED golden record, not the observed one
                for rr in list(direct.values()) + list(indirect.values()):
                    if rr.gene and rr.gene.upper() == other_genes[0].upper():
                        golden.setdefault(rr.accession, {"gene": rr.gene, "sources": set()})["sources"].add("publications")
        else:
            # freeform / non-templated sentence: no reliable exact match available.
            # Best-effort keyword overlap, reported as a REVIEW candidate, not an
            # auto-confirmed correction -- this category needs a human (or Claude,
            # reading the sentence) to confirm, deliberately, per the brief's
            # "resolvable entirely from context_sentence text" note.
            ptoks = tokens_of(phrase_text)
            own_toks = gene_token_sets.get(mention.upper(), set())
            own_overlap = len(ptoks & own_toks)
            best_gene, best_overlap = None, 0
            for g, toks in gene_token_sets.items():
                if g.upper() == mention.upper():
                    continue
                ov = len(ptoks & toks)
                if ov > best_overlap:
                    best_gene, best_overlap = g, ov
            if best_gene and best_overlap >= 2 and best_overlap > own_overlap:
                add_finding(
                    mention, f"target_mention='{mention}' with context_sentence: \"{sentence}\"",
                    f"candidate: {best_gene} (needs manual confirmation)",
                    f"context_sentence shares distinctive terms with gene {best_gene}'s "
                    f"EBI-retrieved protein name/synonyms, and shares none/fewer with "
                    f"'{mention}''s. pmid was not fetched. This is a NEEDS-REVIEW candidate, "
                    f"not an auto-confirmed correction -- read the sentence to confirm.",
                    "EBI Proteins API protein-name lookup (context_sentence text only)",
                    "needs_review", "ambiguous_mention_candidate",
                )
                # left unattributed pending manual confirmation -- do not silently
                # attach an unverified mention to either candidate's golden record
            elif mention.upper() in gene_token_sets or own_toks:
                for rr in list(direct.values()) + list(indirect.values()):
                    if rr.gene and rr.gene.upper() == mention.upper():
                        golden.setdefault(rr.accession, {"gene": rr.gene, "sources": set()})["sources"].add("publications")

    # ---- 6. assemble output ----
    golden_records = [
        {"gene": g["gene"], "primary_accession": acc, "sources": sorted(g["sources"])}
        for acc, g in sorted(golden.items())
    ]

    output = {
        "unique_target_count": len(golden_records),
        "golden_records": golden_records,
        "findings": list(finding_registry.values()),
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
