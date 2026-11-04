#!/usr/bin/env python3
"""Check source identity, original-claim coverage, and manuscript references."""
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


def check(root):
    root = Path(root)
    provenance = json.loads((root / "source/provenance.json").read_text())
    digest = hashlib.sha256((root / "source/original.pdf").read_bytes()).hexdigest()
    assert digest == provenance["sha256"]
    source = "\n".join(p.read_text() for p in (root / "paper").rglob("*.tex"))
    labels = re.findall(r"\\label\{([^}]+)\}", source)
    assert all(n == 1 for n in Counter(labels).values()), "Duplicate label"
    refs = re.findall(r"\\(?:eq)?ref\{([^}]+)\}", source)
    assert not set(refs) - set(labels), sorted(set(refs) - set(labels))
    assert not re.search(r"(?<!\\)qquad", source), "Literal qquad"
    assert "\\tag{\\ref" not in source
    bib = (root / "paper/references.bib").read_text()
    keys = set(re.findall(r"@\w+\{([^,]+),", bib))
    cites = set()
    for group in re.findall(r"\\cite(?:[pt])?(?:\[[^]]*\])?\{([^}]+)\}", source):
        cites.update(group.split(","))
    assert not cites - keys, sorted(cites - keys)
    audit = "\n".join(p.read_text() for p in (root / "audit").glob("*.md"))
    covered = set(re.findall(r"(?<![\d.])\d+\.\d+(?![\d.])", audit))
    for m in re.finditer(r"(?<![\d.])(\d+)\.(\d+)\s*(?:--|–|—|至|到)\s*(?:(\d+)\.)?(\d+)", audit):
        section, start, end_section, end = int(m[1]), int(m[2]), int(m[3] or m[1]), int(m[4])
        if section == end_section and 0 <= end - start < 50:
            covered.update(f"{section}.{i}" for i in range(start, end + 1))
    with (root / "audit/original-claims.csv").open() as stream:
        claims = list(csv.DictReader(stream))
    missing = sorted({row["id"] for row in claims} - covered)
    assert not missing, missing
    return {"original_sha256": digest, "original_named_items": len(claims),
            "items_with_audit_mapping": len(claims), "labels": len(labels),
            "references": len(refs), "cited_bibliography_keys": len(cites),
            "bibliography_entries": len(keys),
            "scope": "Structural coverage and reference consistency; mathematical proofs are reviewed separately."}


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    result = check(root)
    target = root / "results/manuscript_checks.json"
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False))
