# Accessing H–I audit records and evidence

The shared audit ledger and PDF corpus live outside an installed Skill. The default commands locate the evidence package through the user-local registry at CODEX_HOME/math-skill-evidence-roots.json, or HOME/.codex/math-skill-evidence-roots.json when CODEX_HOME is unset. They work from any current working directory and do not depend on EVIDENCE_ROOT.

The registry route is explicit: every resolver call supplies --root-id and --registry, and runs through env -u EVIDENCE_ROOT. An inherited or stale EVIDENCE_ROOT therefore cannot redirect a lookup. Every shell block uses fail-fast handling; if an artifact or PDF lookup fails, later commands in that block do not run.

## Resolve the shared H–I evidence map, claims, and cases

Run this block as written. It discovers the unique math64-20261007 package root, checks the required bundle directories and this Skill's source index, then resolves the three shared H–I audit artifacts in order.

~~~sh
set -eu
MATH64_CODEX_HOME="${CODEX_HOME:-$HOME/.codex}"
MATH64_ROOT_ID="math64-20261007"
MATH64_REGISTRY="${MATH64_REGISTRY:-$MATH64_CODEX_HOME/math-skill-evidence-roots.json}"
MATH64_SKILL_ID="math-i04-random-carleson-occupancy-llogl"

MATH64_PACKAGE_ROOT="$(python3 - "$MATH64_REGISTRY" "$MATH64_ROOT_ID" <<'PYROOT'
import json, os, sys
from pathlib import Path
registry = Path(sys.argv[1]).expanduser()
try:
    doc = json.loads(registry.read_text(encoding="utf-8"))
    roots = doc.get("roots", [])
    matches = [row for row in roots if isinstance(row, dict) and row.get("id") == sys.argv[2]]
    if doc.get("schema") != "math-skill-evidence-root-registry-v1" or len(matches) != 1:
        raise ValueError("missing or ambiguous registry root ID")
    raw = os.path.expandvars(os.path.expanduser(str(matches[0].get("package_root", ""))))
    if not raw or "$" in raw:
        raise ValueError("package_root is empty or contains an unresolved variable")
    root = Path(raw).resolve()
    missing = [name for name in ("workspace_revised", "audit_current", "evidence") if not (root / name).is_dir()]
    if missing:
        raise ValueError("bundle root is missing " + ", ".join(missing))
    print(root)
except Exception as exc:
    raise SystemExit("configuration gap: " + str(exc))
PYROOT
)" || { echo "Could not resolve the selected registry root." >&2; exit 2; }

MATH64_RESOLVER="$MATH64_PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
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
  *) echo "Set MATH64_SKILL_DIR to this installed Skill's absolute directory." >&2; exit 2 ;;
esac
if [ "$(basename "$MATH64_SKILL_DIR")" != "$MATH64_SKILL_ID" ] || [ ! -f "$MATH64_SKILL_DIR/references/handoff-evidence-index.csv" ]; then
  echo "This Skill's handoff-evidence-index.csv was not found at $MATH64_SKILL_DIR" >&2
  exit 2
fi
if [ ! -f "$MATH64_RESOLVER" ]; then
  echo "Evidence resolver was not found: $MATH64_RESOLVER" >&2
  exit 2
fi
unset EVIDENCE_ROOT

env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" --artifact-path audit_current/hi/skill-step-evidence-map.md
env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" --artifact-path audit_current/hi/claims.json
env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" --artifact-path audit_current/hi/cases.json
~~~

## Resolve a paper or evidence PDF

Use the exact P-… paper ID or E… evidence ID shown in this Skill's references/handoff-evidence-index.csv. This block repeats registry discovery so it runs independently from any working directory. Set the identifier you need and leave the other empty.

~~~sh
set -eu
MATH64_CODEX_HOME="${CODEX_HOME:-$HOME/.codex}"
MATH64_ROOT_ID="math64-20261007"
MATH64_REGISTRY="${MATH64_REGISTRY:-$MATH64_CODEX_HOME/math-skill-evidence-roots.json}"
MATH64_SKILL_ID="math-i04-random-carleson-occupancy-llogl"

MATH64_PACKAGE_ROOT="$(python3 - "$MATH64_REGISTRY" "$MATH64_ROOT_ID" <<'PYROOT'
import json, os, sys
from pathlib import Path
registry = Path(sys.argv[1]).expanduser()
try:
    doc = json.loads(registry.read_text(encoding="utf-8"))
    roots = doc.get("roots", [])
    matches = [row for row in roots if isinstance(row, dict) and row.get("id") == sys.argv[2]]
    if doc.get("schema") != "math-skill-evidence-root-registry-v1" or len(matches) != 1:
        raise ValueError("missing or ambiguous registry root ID")
    raw = os.path.expandvars(os.path.expanduser(str(matches[0].get("package_root", ""))))
    if not raw or "$" in raw:
        raise ValueError("package_root is empty or contains an unresolved variable")
    root = Path(raw).resolve()
    missing = [name for name in ("workspace_revised", "audit_current", "evidence") if not (root / name).is_dir()]
    if missing:
        raise ValueError("bundle root is missing " + ", ".join(missing))
    print(root)
except Exception as exc:
    raise SystemExit("configuration gap: " + str(exc))
PYROOT
)" || { echo "Could not resolve the selected registry root." >&2; exit 2; }

MATH64_RESOLVER="$MATH64_PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
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
  *) echo "Set MATH64_SKILL_DIR to this installed Skill's absolute directory." >&2; exit 2 ;;
esac
if [ "$(basename "$MATH64_SKILL_DIR")" != "$MATH64_SKILL_ID" ] || [ ! -f "$MATH64_SKILL_DIR/references/handoff-evidence-index.csv" ]; then
  echo "This Skill's handoff-evidence-index.csv was not found at $MATH64_SKILL_DIR" >&2
  exit 2
fi
if [ ! -f "$MATH64_RESOLVER" ]; then
  echo "Evidence resolver was not found: $MATH64_RESOLVER" >&2
  exit 2
fi
unset EVIDENCE_ROOT

PAPER_ID=""
EVIDENCE_ID=""
if [ -n "$PAPER_ID" ]; then
  env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" --skill-dir "$MATH64_SKILL_DIR" --paper-id "$PAPER_ID" --verify-sha256
fi
if [ -n "$EVIDENCE_ID" ]; then
  env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" --skill-dir "$MATH64_SKILL_DIR" --evidence-id "$EVIDENCE_ID" --verify-sha256
fi
if [ -z "$PAPER_ID" ] && [ -z "$EVIDENCE_ID" ]; then
  echo "Set PAPER_ID or EVIDENCE_ID to an exact index value before running this block." >&2
  exit 2
fi
~~~

## Deliberate explicit package-root route

If you intentionally use a package path instead of the user-local registry, set MATH64_SOURCE_PACKAGE_ROOT to the absolute directory containing workspace_revised/, audit_current/, and evidence/. Do not use EVIDENCE_ROOT for this route. The resolver call passes --package-root and clears EVIDENCE_ROOT so an inherited value cannot silently win.

~~~sh
set -eu
: "${MATH64_SOURCE_PACKAGE_ROOT:?Set MATH64_SOURCE_PACKAGE_ROOT to an absolute evidence-package root}"
case "$MATH64_SOURCE_PACKAGE_ROOT" in
  /*) ;;
  *) echo "MATH64_SOURCE_PACKAGE_ROOT must be absolute." >&2; exit 2 ;;
esac
MATH64_SKILL_ID="math-i04-random-carleson-occupancy-llogl"
MATH64_SKILL_DIR="${MATH64_SKILL_DIR:-$MATH64_SOURCE_PACKAGE_ROOT/workspace_revised/$MATH64_SKILL_ID}"
MATH64_RESOLVER="$MATH64_SOURCE_PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
case "$MATH64_SKILL_DIR" in
  /*) ;;
  *) echo "MATH64_SKILL_DIR must be absolute." >&2; exit 2 ;;
esac
if [ "$(basename "$MATH64_SKILL_DIR")" != "$MATH64_SKILL_ID" ]; then
  echo "MATH64_SKILL_DIR does not name $MATH64_SKILL_ID." >&2
  exit 2
fi
if [ ! -f "$MATH64_RESOLVER" ] || [ ! -f "$MATH64_SKILL_DIR/references/handoff-evidence-index.csv" ]; then
  echo "Explicit package root is missing the resolver or Skill index." >&2
  exit 2
fi
unset EVIDENCE_ROOT
env -u EVIDENCE_ROOT python3 "$MATH64_RESOLVER" --package-root "$MATH64_SOURCE_PACKAGE_ROOT" --artifact-path audit_current/hi/skill-step-evidence-map.md
~~~

If the registry, root ID, complete bundle, resolver, or index is missing, stop and report the configuration gap; do not infer paths from the current directory or invent a registry entry. The shipped audit_current/delivery/source-root-registry.example.json is a template, not a user-local registry. Resolver success verifies path and current bytes only; it does not certify full-text reading, a theorem proof, or mathematical acceptance.
