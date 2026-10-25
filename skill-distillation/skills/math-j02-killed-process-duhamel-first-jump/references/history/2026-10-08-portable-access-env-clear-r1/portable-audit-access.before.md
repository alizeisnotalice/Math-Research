# Portable source and audit access for `math-j02-killed-process-duhamel-first-jump`

An installed Skill folder does not contain the shared PDF corpus or package audit records. Resolve both through the package's `audit_current/delivery/resolve_evidence.py`; do not assemble a path from this Skill's directory or depend on the current working directory. Keep one complete bundle root containing `workspace_revised/`, `audit_current/`, and `evidence/` available. The resolver checks indexed PDF bytes and reports the artifact SHA.

## Discover the bundle root from any directory

Set `MATH64_SKILL_DIR` when this Skill is installed somewhere other than the default directory. For a Skill used from the unpacked package, set it to `$EVIDENCE_ROOT/workspace_revised/math-j02-killed-process-duhamel-first-jump`. Use the user-local registry entry `math64-20261007` in `${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json` first. The registry must have the schema `math-skill-evidence-root-registry-v1` and exactly one matching ID with a concrete absolute `package_root`. If the registry file is absent, fall back to explicit `EVIDENCE_ROOT` or `MATH64_SOURCE_PACKAGE_ROOT`; if a readable registry is malformed or ambiguous, stop instead of overriding it with an environment value.

```sh
set -eu
MATH64_ROOT_ID="math64-20261007"
MATH64_CODEX_HOME="${CODEX_HOME:-$HOME/.codex}"
MATH64_REGISTRY="${MATH64_CODEX_HOME}/math-skill-evidence-roots.json"
MATH64_SKILL_NAME="math-j02-killed-process-duhamel-first-jump"
MATH64_SKILL_DIR="${MATH64_SKILL_DIR:-}"

if [ -e "$MATH64_REGISTRY" ]; then
  if [ ! -r "$MATH64_REGISTRY" ]; then
    printf '%s\n' "Evidence root registry exists but is unreadable: $MATH64_REGISTRY"
    exit 2
  fi
  MATH64_PACKAGE_ROOT="$(python3 - "$MATH64_REGISTRY" "$MATH64_ROOT_ID" <<'PYROOT'
import json, os, re, sys
from pathlib import Path
registry = Path(sys.argv[1]).expanduser()
with registry.open(encoding="utf-8") as stream:
    doc = json.load(stream)
rows = [row for row in doc.get("roots", []) if isinstance(row, dict) and row.get("id") == sys.argv[2]]
if doc.get("schema") != "math-skill-evidence-root-registry-v1" or len(rows) != 1:
    raise SystemExit("Registry schema or math64-20261007 entry is missing/ambiguous.")
raw = os.path.expandvars(os.path.expanduser(str(rows[0].get("package_root", ""))))
if not raw or re.search(r"\$[A-Za-z_{]", raw) or not Path(raw).is_absolute():
    raise SystemExit("Registry package_root must be a concrete absolute path.")
print(Path(raw).resolve())
PYROOT
  )" || exit 2
  MATH64_ROOT_MODE="registry"
elif [ -n "${EVIDENCE_ROOT:-}" ]; then
  MATH64_PACKAGE_ROOT="$EVIDENCE_ROOT"
  MATH64_ROOT_MODE="package-root"
elif [ -n "${MATH64_SOURCE_PACKAGE_ROOT:-}" ]; then
  MATH64_PACKAGE_ROOT="$MATH64_SOURCE_PACKAGE_ROOT"
  MATH64_ROOT_MODE="package-root"
else
  printf '%s\n' "Bundle root is not configured; set EVIDENCE_ROOT/MATH64_SOURCE_PACKAGE_ROOT or configure the registry."
  exit 2
fi

case "$MATH64_PACKAGE_ROOT" in
  /*) ;;
  *) printf '%s\n' "Bundle root must be absolute: $MATH64_PACKAGE_ROOT"; exit 2 ;;
esac
for MATH64_REQUIRED_DIR in workspace_revised audit_current evidence; do
  if [ ! -d "$MATH64_PACKAGE_ROOT/$MATH64_REQUIRED_DIR" ]; then
    printf '%s\n' "Bundle root is missing $MATH64_REQUIRED_DIR/: $MATH64_PACKAGE_ROOT"
    exit 2
  fi
done
if [ -z "$MATH64_SKILL_DIR" ]; then
  if [ -f "$MATH64_PACKAGE_ROOT/workspace_revised/$MATH64_SKILL_NAME/references/handoff-evidence-index.csv" ]; then
    MATH64_SKILL_DIR="$MATH64_PACKAGE_ROOT/workspace_revised/$MATH64_SKILL_NAME"
  else
    MATH64_SKILL_DIR="$MATH64_CODEX_HOME/skills/$MATH64_SKILL_NAME"
  fi
fi
MATH64_RESOLVER="$MATH64_PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
if [ ! -f "$MATH64_RESOLVER" ] || [ ! -f "$MATH64_SKILL_DIR/references/handoff-evidence-index.csv" ]; then
  printf '%s\n' "Resolver or this Skill's evidence index is unavailable; set MATH64_SKILL_DIR to the installed/package Skill directory."
  exit 2
fi
```

Run the bootstrap block and one resolver block in the same shell session so its variables are available. Each resolver block uses fail-fast shell options. For a registry-selected root, retain `--root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY"` in the resolver command. For either environment-selected root, pass `--package-root "$MATH64_PACKAGE_ROOT"`. A packaged Skill can use `$MATH64_PACKAGE_ROOT/workspace_revised/math-j02-killed-process-duhamel-first-jump` as `MATH64_SKILL_DIR`; an isolated installed Skill normally uses `$MATH64_CODEX_HOME/skills/math-j02-killed-process-duhamel-first-jump`.

## Resolve and verify an indexed PDF

Use the exact `paper_id` or `evidence_id` from this Skill's `references/handoff-evidence-index.csv`. `--verify-sha256` makes the resolver hash the local PDF and compare it with the indexed SHA before returning its path.

```sh
set -eu
PAPER_ID="P-…"
if [ "$MATH64_ROOT_MODE" = registry ]; then
  python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" \
    --skill-dir "$MATH64_SKILL_DIR" --paper-id "$PAPER_ID" --verify-sha256
else
  python3 "$MATH64_RESOLVER" --package-root "$MATH64_PACKAGE_ROOT" \
    --skill-dir "$MATH64_SKILL_DIR" --paper-id "$PAPER_ID" --verify-sha256
fi
```

Use `--evidence-id E…` instead when the evidence row is the stable locator. If a row has no local portable source or its SHA does not match, retain the unresolved/mismatch status; do not infer that the source was read.

## Resolve a shared audit artifact

Pass the package-relative path printed in the Skill's evidence notes. The resolver rejects absolute paths and `..`, then returns the path and current artifact SHA-256.

```sh
set -eu
ARTIFACT_PATH="audit_current/…/record.json"
if [ "$MATH64_ROOT_MODE" = registry ]; then
  python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" \
    --artifact-path "$ARTIFACT_PATH"
else
  python3 "$MATH64_RESOLVER" --package-root "$MATH64_PACKAGE_ROOT" \
    --artifact-path "$ARTIFACT_PATH"
fi
```

A successful resolution verifies locator and bytes only. It does not establish bibliographic identity, completeness of a reading, theorem validity, or mathematical applicability. Keep the shared PDF corpus and ledgers at the package root rather than duplicating them into each installed Skill.
