# Portable case and source evidence access

<a id="case-ledger"></a>

This Skill's source index is `references/handoff-evidence-index.csv`. The C/D claims and case ledgers are the bundle artifacts `audit_current/cd/claims.json` and `audit_current/cd/cases.json`; their portable artifact keys are independent of this Skill's install path. The shared resolver and PDFs live at the root of the complete evidence bundle, alongside `workspace_revised/`, `audit_current/`, and `evidence/`. A Skill-only installation does not include those shared files. `EVIDENCE_ROOT`, when used, means the complete bundle root, not its `evidence/` directory. Never resolve a bundle path from the current working directory. Both access routes below use the one shared `audit_current/delivery/resolve_evidence.py`; do not copy a resolver into the Skill.

The registry route below works from any working directory. It uses the user-local registry `${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json` and requires one concrete entry with ID `math64-20261007`; set `MATH64_REGISTRY` or `MATH64_ROOT_ID` only when your local registry uses another file or exact root ID. The shipped `source-root-registry.example.json` is a template, not a live registry. These commands resolve the current C/D claims and case ledgers and verify a real PDF SHA indexed by this Skill. Change `MATH64_EVIDENCE_ID` only to another exact ID from this Skill’s local index when needed.

```sh
set -eu
MATH64_CODEX_HOME="${CODEX_HOME:-$HOME/.codex}"
MATH64_REGISTRY="${MATH64_REGISTRY:-$MATH64_CODEX_HOME/math-skill-evidence-roots.json}"
MATH64_ROOT_ID="${MATH64_ROOT_ID:-math64-20261007}"
MATH64_SKILL_ID="math-d01-kantorovich-duality"
MATH64_EVIDENCE_ID="E1066"

if ! MATH64_PACKAGE_ROOT="$(python3 - "$MATH64_REGISTRY" "$MATH64_ROOT_ID" <<'PYROOT'
import json, os, sys
from pathlib import Path
registry = Path(sys.argv[1]).expanduser()
doc = json.loads(registry.read_text(encoding="utf-8"))
if doc.get("schema") != "math-skill-evidence-root-registry-v1":
    raise SystemExit("unsupported evidence-root registry schema")
roots = [r for r in doc.get("roots", []) if isinstance(r, dict) and r.get("id") == sys.argv[2]]
if len(roots) != 1:
    raise SystemExit("root ID is missing or ambiguous in the registry")
raw = os.path.expandvars(os.path.expanduser(str(roots[0].get("package_root", ""))))
if not raw or "$" in raw:
    raise SystemExit("package_root is empty or contains an unresolved variable")
root = Path(raw).resolve()
missing = [name for name in ("workspace_revised", "audit_current", "evidence") if not (root / name).is_dir()]
if missing:
    raise SystemExit("bundle root is missing: " + ", ".join(missing))
print(root)
PYROOT
)"; then
  echo "Could not resolve the selected evidence bundle from the registry." >&2
  exit 2
fi
MATH64_RESOLVER="$MATH64_PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
if [ ! -f "$MATH64_RESOLVER" ]; then
  echo "Evidence resolver is missing: $MATH64_RESOLVER" >&2
  exit 2
fi
MATH64_SKILL_DIR="${MATH64_SKILL_DIR:-}"
if [ -z "$MATH64_SKILL_DIR" ]; then
  for MATH64_CANDIDATE in "$MATH64_CODEX_HOME/skills/$MATH64_SKILL_ID" "$HOME/.agents/skills/$MATH64_SKILL_ID" "$MATH64_PACKAGE_ROOT/workspace_revised/$MATH64_SKILL_ID"; do
    if [ -f "$MATH64_CANDIDATE/references/handoff-evidence-index.csv" ]; then
      MATH64_SKILL_DIR="$MATH64_CANDIDATE"
      break
    fi
  done
fi
case "$MATH64_SKILL_DIR" in
  /*) ;;
  *) echo "Set MATH64_SKILL_DIR to this Skill's absolute directory." >&2; exit 2 ;;
esac
if [ "$(basename "$MATH64_SKILL_DIR")" != "$MATH64_SKILL_ID" ] || [ ! -f "$MATH64_SKILL_DIR/references/handoff-evidence-index.csv" ]; then
  echo "This Skill's evidence index was not found at $MATH64_SKILL_DIR" >&2
  exit 2
fi

# Every registry lookup explicitly ignores a stale/conflicting EVIDENCE_ROOT.
env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" --artifact-path audit_current/cd/claims.json
env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" --artifact-path audit_current/cd/cases.json
env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" --skill-dir "$MATH64_SKILL_DIR" --evidence-id "$MATH64_EVIDENCE_ID" --verify-sha256
```

## Explicit package-root route

When a registry is unavailable, supply the absolute bundle path directly. This route also works from any working directory, and `--package-root` takes precedence over a conflicting `EVIDENCE_ROOT`; `env -u` makes that choice explicit. Set `MATH64_SOURCE_PACKAGE_ROOT` to the directory that contains `workspace_revised/`, `audit_current/`, and `evidence/`.

```sh
set -eu
: "${MATH64_SOURCE_PACKAGE_ROOT:?Set an absolute complete evidence-bundle root}"
case "$MATH64_SOURCE_PACKAGE_ROOT" in
  /*) ;;
  *) echo "MATH64_SOURCE_PACKAGE_ROOT must be absolute." >&2; exit 2 ;;
esac
MATH64_ROOT="$(python3 - "$MATH64_SOURCE_PACKAGE_ROOT" <<'PYROOT'
from pathlib import Path
import sys
root = Path(sys.argv[1]).expanduser().resolve()
missing = [name for name in ("workspace_revised", "audit_current", "evidence") if not (root / name).is_dir()]
if missing:
    raise SystemExit("bundle root is missing: " + ", ".join(missing))
print(root)
PYROOT
)"
MATH64_RESOLVER="$MATH64_ROOT/audit_current/delivery/resolve_evidence.py"
MATH64_SKILL_ID="math-d01-kantorovich-duality"
MATH64_SKILL_DIR="${MATH64_SKILL_DIR:-$MATH64_ROOT/workspace_revised/$MATH64_SKILL_ID}"
MATH64_EVIDENCE_ID="E1066"
if [ ! -f "$MATH64_RESOLVER" ] || [ "$(basename "$MATH64_SKILL_DIR")" != "$MATH64_SKILL_ID" ] || [ ! -f "$MATH64_SKILL_DIR/references/handoff-evidence-index.csv" ]; then
  echo "The explicit bundle root lacks this Skill index or the resolver." >&2
  exit 2
fi
env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --package-root "$MATH64_ROOT" --artifact-path audit_current/cd/claims.json
env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --package-root "$MATH64_ROOT" --artifact-path audit_current/cd/cases.json
env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --package-root "$MATH64_ROOT" --skill-dir "$MATH64_SKILL_DIR" --evidence-id "$MATH64_EVIDENCE_ID" --verify-sha256
```

If either route cannot find a unique registry ID, a complete bundle root, the resolver, this Skill's index, or the indexed PDF bytes, stop and report the exact missing item. Do not guess paths or substitute a different PDF/version. The resolver checks current paths and source SHA only; it does not establish correct bibliographic identity, complete reading, theorem validity, or mathematical acceptance.
