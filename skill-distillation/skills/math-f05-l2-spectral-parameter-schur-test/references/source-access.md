# Indexed source resolution

The handoff index beside this Skill records each evidence ID, package-relative source path, and expected SHA-256. Resolve through the resolver shipped with the complete evidence package; do not construct paths relative to an installed Skill. The package root is the directory containing `workspace_revised/`, `audit_current/`, and `evidence/`.

For a Skill inside the extracted package, set `PACKAGE_ROOT` to that directory and use its Skill directory:

```sh
PACKAGE_ROOT="/path/to/extracted/math64-package"
python3 "$PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py" \
  --package-root "$PACKAGE_ROOT" \
  --skill-dir "$PACKAGE_ROOT/workspace_revised/math-fXX-skill-name" \
  --evidence-id E0000
```

Replace the example evidence ID and Skill directory with exact values from this Skill’s `references/handoff-evidence-index.csv`. These F indexes use `evidence_id`; use `--paper-id P-<16-hex-digits>` only with an index that has a `paper_id` column. An individual lookup hashes the source and compares it with the indexed SHA automatically; `--verify-sha256` is needed only for a bundle-wide `--audit-all` run.

An installed Skill still needs access to an extracted full package. Export `MATH64_SOURCE_PACKAGE_ROOT` with that package root and configure the user-local registry using the package template `audit_current/delivery/source-root-registry.example.json` (schema `math-skill-evidence-root-registry-v1`, root ID `math64-20261007`, package root `${MATH64_SOURCE_PACKAGE_ROOT}`). To exercise the registry fallback, unset `EVIDENCE_ROOT`, omit `--package-root`, and pass both `--root-id math64-20261007` and `--registry <registry-file>`:

```sh
unset EVIDENCE_ROOT
export MATH64_SOURCE_PACKAGE_ROOT="/path/to/extracted/math64-package"
python3 "$MATH64_SOURCE_PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py" \
  --root-id math64-20261007 \
  --registry "/path/to/math-skill-evidence-roots.json" \
  --skill-dir "/path/to/installed/math-fXX-skill-name" \
  --evidence-id E0000
```

Resolution precedence is `--package-root`, `EVIDENCE_ROOT`, the requested registry root ID, then package-internal resolver auto-detection. Auto-detection is not an installed-Skill setup. The resolver rejects unset root variables, unsafe or missing paths, missing sources, and SHA mismatches. A successful lookup verifies source identity by bytes and path only; it does not show that the source was read, that a theorem was proved, or that the Skill is installed.
