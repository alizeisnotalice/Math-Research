# Accessing shared F audit records from an installed Skill

This standalone Skill refers to shared package files such as `audit_current/f/production_audit_20261007.md`, `audit_current/f/claims.json`, `audit_current/f/cases.json`, and `audit_current/f/papers-read.json`. They live at the bundle root, not inside the installed Skill folder. Keep the unpacked full bundle available; its root contains `workspace_revised/`, `audit_current/`, and `evidence/`.

## Discover the bundle root without using the current directory

When `${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json` exists, its unique `math64-20261007` row is authoritative, even if an inherited `EVIDENCE_ROOT` points somewhere else. The registry's `package_root` must be one concrete absolute path. If no registry file exists, set `EVIDENCE_ROOT` or `MATH64_SOURCE_PACKAGE_ROOT` to the absolute bundle root. If the configured registry is malformed, missing the ID, ambiguous, unresolved, or points to an incomplete bundle, stop and report that configuration issue; do not silently fall back to another root.

Set `MATH64_SKILL_DIR` only when this Skill is installed somewhere other than the default Codex path. The discovery block works from any working directory. It uses `set -eu` so a failed lookup or check stops before later commands run.

```sh
set -eu
ROOT_REGISTRY="${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json"
ROOT_ID="math64-20261007"
MATH64_SKILL_DIR="${MATH64_SKILL_DIR:-${CODEX_HOME:-$HOME/.codex}/skills/math-f04-poisson-gamma-abel-square-function}"
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
  ROOT_MODE="registry"
elif [ -n "${EVIDENCE_ROOT:-}" ]; then
  ROOT_CANDIDATE="$EVIDENCE_ROOT"
  ROOT_MODE="environment"
elif [ -n "${MATH64_SOURCE_PACKAGE_ROOT:-}" ]; then
  ROOT_CANDIDATE="$MATH64_SOURCE_PACKAGE_ROOT"
  ROOT_MODE="environment"
else
  echo "No bundle registry found; set EVIDENCE_ROOT or MATH64_SOURCE_PACKAGE_ROOT to the absolute bundle root." >&2
  exit 1
fi
if [ "$ROOT_MODE" = environment ]; then
  EVIDENCE_ROOT="$(python3 - "$ROOT_CANDIDATE" <<'PYROOT'
import os, re, sys
from pathlib import Path
raw = os.path.expandvars(os.path.expanduser(sys.argv[1]))
if not raw or re.search(r"\$[A-Za-z_{]", raw):
    raise SystemExit("Bundle root is empty or contains an unset variable.")
root = Path(raw).expanduser()
if not root.is_absolute():
    raise SystemExit("Bundle root must be absolute; current-directory resolution is disabled.")
root = root.resolve()
missing = [name for name in ("workspace_revised", "audit_current", "evidence") if not (root / name).is_dir()]
if missing:
    raise SystemExit("Bundle root is missing: " + ", ".join(missing))
print(root)
PYROOT
  )"
fi
export EVIDENCE_ROOT ROOT_REGISTRY ROOT_ID ROOT_MODE MATH64_SKILL_DIR
RESOLVER="$EVIDENCE_ROOT/audit_current/delivery/resolve_evidence.py"
```

Run the discovery block first. If the registry exists, it overrides any ambient root. In later registry-mode commands, `env -u EVIDENCE_ROOT` makes the resolver use the same registry row selected above. When no registry exists, pass `--package-root "$EVIDENCE_ROOT"` explicitly. Do not run resolver commands from a Skill folder by relying on the current directory.

## Resolve current audit records

Each resolver call prints the located path and current SHA-256. The resolver rejects absolute artifact paths and `..` traversal. Keep `set -eu` in each block: a failed first lookup must stop the block, rather than be hidden by a successful later lookup.

```sh
set -eu
if [ "$ROOT_MODE" = registry ]; then
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$ROOT_REGISTRY" --artifact-path audit_current/f/production_audit_20261007.md
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$ROOT_REGISTRY" --artifact-path audit_current/f/claims.json
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$ROOT_REGISTRY" --artifact-path audit_current/f/cases.json
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$ROOT_REGISTRY" --artifact-path audit_current/f/papers-read.json
else
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --package-root "$EVIDENCE_ROOT" --artifact-path audit_current/f/production_audit_20261007.md
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --package-root "$EVIDENCE_ROOT" --artifact-path audit_current/f/claims.json
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --package-root "$EVIDENCE_ROOT" --artifact-path audit_current/f/cases.json
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --package-root "$EVIDENCE_ROOT" --artifact-path audit_current/f/papers-read.json
fi
```

Other F records, including `audit_current/f/step_coverage_20261007.md` and files under `audit_current/f/source-recovery/`, use the same `--artifact-path` form.

## Resolve an indexed source PDF

Set `EVIDENCE_ID` to an exact evidence row ID (for example `E0000`) in this Skill's `references/handoff-evidence-index.csv`. The resolver reads the absolute installed Skill directory and checks the indexed PDF SHA before returning its path:

```sh
set -eu
EVIDENCE_ID="E0000"
if [ "$ROOT_MODE" = registry ]; then
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$ROOT_REGISTRY" --skill-dir "$MATH64_SKILL_DIR" --evidence-id "$EVIDENCE_ID" --verify-sha256
else
  env -u EVIDENCE_ROOT python3 "$RESOLVER" --package-root "$EVIDENCE_ROOT" --skill-dir "$MATH64_SKILL_DIR" --evidence-id "$EVIDENCE_ID" --verify-sha256
fi
```

For a different source, set `EVIDENCE_ID` to the exact evidence-row ID in this Skill’s `references/handoff-evidence-index.csv` (not the paper ID); keep the resolver option `--evidence-id "$EVIDENCE_ID"`. To inspect another shared package file, use `--artifact-path` with its package-relative path. A resolver pass verifies path and current bytes only; it does not establish that a source was read, a theorem proved, or a Skill installed. Keep shared ledgers and PDFs at the bundle root rather than duplicating them inside each Skill.
