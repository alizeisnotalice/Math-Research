# math-g04-rubio-de-francia-square-function-scope-check portable audit access

The shared evidence corpus is one complete package root with sibling `workspace_revised/`, `audit_current/`, and `evidence/` directories. `EVIDENCE_ROOT` and `MATH64_SOURCE_PACKAGE_ROOT` mean that package root. A standalone installed Skill does not include the shared source PDFs or audit ledgers. Do not infer a path from the working directory or a WorkBuddy cache.

## Locate the absolute resolver through the user-local registry

The resolver lives at `<package_root>/audit_current/delivery/resolve_evidence.py`. To locate it from an installed Skill without knowing the bundle path first, read the unique `math64-20261007` entry from `${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json`, expand its `package_root`, then pass the same registry and root ID to the resolver:

```sh
set -e
REGISTRY="${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json"
PACKAGE_ROOT="$(python3 - "$REGISTRY" <<'PY'
import json, os, re, sys
from pathlib import Path
registry = Path(sys.argv[1]).expanduser()
try:
    doc = json.loads(registry.read_text(encoding="utf-8"))
    rows = doc.get("roots", [])
    matches = [r for r in rows if isinstance(r, dict) and r.get("id") == "math64-20261007"]
    if doc.get("schema") != "math-skill-evidence-root-registry-v1" or len(matches) != 1:
        raise ValueError("missing/ambiguous registry root math64-20261007")
    value = os.path.expandvars(os.path.expanduser(str(matches[0].get("package_root", ""))))
    if not value or re.search(r"\$[A-Za-z_{]", value):
        raise ValueError("package_root is empty or has an unset environment variable")
    print(Path(value).resolve())
except Exception as exc:
    raise SystemExit(f"configuration gap: {exc}")
PY
)" || exit 2
RESOLVER="$PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
SKILL_DIR="$PACKAGE_ROOT/workspace_revised/math-g04-rubio-de-francia-square-function-scope-check"
test -f "$RESOLVER" && test -f "$SKILL_DIR/references/handoff-evidence-index.csv" || { echo "configuration gap: bundle resolver or Skill index missing" >&2; exit 2; }
env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id math64-20261007 --registry "$REGISTRY" \
  --skill-dir "$SKILL_DIR" --evidence-id E0422
env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id math64-20261007 --registry "$REGISTRY" \
  --artifact-path audit_current/g/claims.json
env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id math64-20261007 --registry "$REGISTRY" \
  --artifact-path audit_current/g/cases.json
```

## Explicit package-root route

When a root is supplied explicitly, the resolver can be derived without using the registry. `EVIDENCE_ROOT` and `MATH64_SOURCE_PACKAGE_ROOT` both name the package root; the first set value below is used:

```sh
set -e
PACKAGE_ROOT="${EVIDENCE_ROOT:-${MATH64_SOURCE_PACKAGE_ROOT:-}}"
test -n "$PACKAGE_ROOT" || { echo "configuration gap: set EVIDENCE_ROOT/MATH64_SOURCE_PACKAGE_ROOT or configure the user-local registry" >&2; exit 2; }
RESOLVER="$PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
SKILL_DIR="$PACKAGE_ROOT/workspace_revised/math-g04-rubio-de-francia-square-function-scope-check"
python3 "$RESOLVER" --package-root "$PACKAGE_ROOT" \
  --skill-dir "$SKILL_DIR" --evidence-id E0422
python3 "$RESOLVER" --package-root "$PACKAGE_ROOT" \
  --artifact-path audit_current/g/claims.json
python3 "$RESOLVER" --package-root "$PACKAGE_ROOT" \
  --artifact-path audit_current/g/cases.json
```

If the registry, root ID, complete bundle, resolver, or index is missing, stop and report the configuration gap; do not invent a package path or registry entry. An unset `${MATH64_SOURCE_PACKAGE_ROOT}` placeholder in a registry is not a usable root. The indexed source lookup verifies the index row against the current PDF SHA. `--artifact-path` returns the artifact path and current SHA so the artifact can be read from that resolved path. `PASS` verifies path resolution and current bytes only; it does not certify full-text reading, a proved claim, or mathematical acceptance.
