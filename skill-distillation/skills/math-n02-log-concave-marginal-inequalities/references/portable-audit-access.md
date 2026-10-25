# Accessing shared evidence from an installed L–N Skill

A standalone Skill directory does not contain the shared audit ledgers or PDF corpus. Keep the complete unpacked bundle root available; the package-level resolver is stored once at `<bundle-root>/audit_current/delivery/resolve_evidence.py`, not copied into each Skill. These commands do not rely on the current working directory. When the user-local registry `${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json` exists, its unique `math64-20261007` record is authoritative and overrides any inherited `EVIDENCE_ROOT`; only when no registry exists may an explicitly set absolute `EVIDENCE_ROOT` select the bundle.

## Discover the bundle root

Set `MATH64_SKILL_DIR` if this Skill is installed outside the default Codex Skills directory. The bundle registry must contain one record whose `id` is `math64-20261007`; set `EVIDENCE_ROOT` directly when no registry is configured. The registry's `package_root` must be absolute and point to the directory containing `workspace_revised/`, `audit_current/`, and `evidence/`. The discovery and resolver blocks use `set -eu`; any missing file or failed resolver call stops the block.

```sh
set -eu
ROOT_REGISTRY="${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json"
SKILL_DIR="${MATH64_SKILL_DIR:-${CODEX_HOME:-$HOME/.codex}/skills/math-n02-log-concave-marginal-inequalities}"
ROOT_ID="math64-20261007"
if [ -f "$ROOT_REGISTRY" ]; then
  EVIDENCE_ROOT="$(python3 - "$ROOT_REGISTRY" "$ROOT_ID" <<'PYROOT'
import json, os, re, sys
from pathlib import Path
registry = Path(sys.argv[1]).expanduser()
root_id = sys.argv[2]
if not registry.is_file():
    raise SystemExit("Evidence-root registry disappeared during discovery.")
doc = json.loads(registry.read_text(encoding="utf-8"))
if doc.get("schema") != "math-skill-evidence-root-registry-v1":
    raise SystemExit("Unsupported evidence-root registry schema.")
rows = [row for row in doc.get("roots", []) if isinstance(row, dict) and row.get("id") == root_id]
if len(rows) != 1:
    raise SystemExit("Registry must resolve math64-20261007 exactly once.")
raw = os.path.expandvars(os.path.expanduser(str(rows[0].get("package_root", ""))))
if not raw or re.search(r"\$[A-Za-z_{]", raw):
    raise SystemExit("Registry package_root is empty or contains an unset variable.")
root = Path(raw).expanduser()
if not root.is_absolute():
    raise SystemExit("Registry package_root must be absolute; current-directory resolution is disabled.")
root = root.resolve()
missing = [name for name in ("workspace_revised", "audit_current", "evidence") if not (root / name).is_dir()]
if missing:
    raise SystemExit("Bundle root is missing: " + ", ".join(missing))
print(root)
PYROOT
  )"
elif [ -n "${EVIDENCE_ROOT:-}" ]; then
  EVIDENCE_ROOT="$(python3 - "$EVIDENCE_ROOT" <<'PYROOT'
import os, re, sys
from pathlib import Path
raw = os.path.expandvars(os.path.expanduser(sys.argv[1]))
if not raw or re.search(r"\$[A-Za-z_{]", raw):
    raise SystemExit("EVIDENCE_ROOT is empty or contains an unset variable.")
root = Path(raw).expanduser()
if not root.is_absolute():
    raise SystemExit("EVIDENCE_ROOT must be absolute; current-directory resolution is disabled.")
root = root.resolve()
missing = [name for name in ("workspace_revised", "audit_current", "evidence") if not (root / name).is_dir()]
if missing:
    raise SystemExit("Bundle root is missing: " + ", ".join(missing))
print(root)
PYROOT
  )"
else
  echo "No bundle registry found; set EVIDENCE_ROOT to the absolute unpacked bundle root." >&2
  exit 1
fi
export EVIDENCE_ROOT
```

If the registry does not exist, create it from `audit_current/delivery/source-root-registry.example.json`, replace the placeholder with the actual absolute bundle root, and save it as `${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json`. Otherwise set `EVIDENCE_ROOT` to the actual absolute root. Run the discovery block first in the same shell; a configured registry always takes precedence over an ambient `EVIDENCE_ROOT`.

## Resolve a shared ledger or source

The resolver accepts only package-relative paths, rejects absolute paths and `..` traversal, and returns the resolved path with its current SHA-256. Resolve shared L–N ledgers as follows:

```sh
set -eu
RESOLVER="$EVIDENCE_ROOT/audit_current/delivery/resolve_evidence.py"
if [ -f "$ROOT_REGISTRY" ]; then
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$ROOT_REGISTRY" --artifact-path audit_current/ln/claims.json
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$ROOT_REGISTRY" --artifact-path audit_current/ln/cases.json
else
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --package-root "$EVIDENCE_ROOT" --artifact-path audit_current/ln/claims.json
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --package-root "$EVIDENCE_ROOT" --artifact-path audit_current/ln/cases.json
fi
```

To locate a PDF, set `PAPER_ID` to an exact `P-…` ID listed in this Skill's `references/handoff-evidence-index.csv`. The resolver reads the installed Skill's index and checks the expected PDF SHA before returning its path:

```sh
set -eu
PAPER_ID="P-…"
RESOLVER="$EVIDENCE_ROOT/audit_current/delivery/resolve_evidence.py"
if [ -f "$ROOT_REGISTRY" ]; then
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$ROOT_REGISTRY" --skill-dir "$SKILL_DIR" --paper-id "$PAPER_ID" --verify-sha256
else
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --package-root "$EVIDENCE_ROOT" --skill-dir "$SKILL_DIR" --paper-id "$PAPER_ID" --verify-sha256
fi
```

Use `--evidence-id E…` when resolving an exact evidence row instead of a paper ID. To resolve another package-level artifact, pass its relative path to `--artifact-path`. The resolver provides path and byte-integrity checks; it does not establish a mathematical claim or make a Skill-only directory self-contained. Keep the bundle root and shared `evidence/` tree available once rather than duplicating them under every Skill.
