# J01 portable evidence lookup

Resolve all artifact entries below from any current working directory using [the portable access guide](portable-audit-access.md). After it defines `MATH64_RESOLVER`, `MATH64_ROOT_MODE`, `MATH64_PACKAGE_ROOT`, and (for registry mode) `MATH64_ROOT_ID` / `MATH64_REGISTRY`, use:

```sh
ARTIFACT_PATH="audit_current/…/record.json"
if [ "$MATH64_ROOT_MODE" = registry ]; then
  python3 "$MATH64_RESOLVER" --root-id "$MATH64_ROOT_ID" --registry "$MATH64_REGISTRY" --artifact-path "$ARTIFACT_PATH"
else
  python3 "$MATH64_RESOLVER" --package-root "$MATH64_PACKAGE_ROOT" --artifact-path "$ARTIFACT_PATH"
fi
```

The resolver returns the current artifact SHA-256. Replace `ARTIFACT_PATH` with the exact package-relative locator printed below; a path/SHA check does not establish mathematical validity.

## J01 E0785 full-read card

- Artifact: `audit_current/shared-reading/audit_jk/500a71307a38c036d97bca3c6b958c710c22abd9ca1479c201a1ed17713b424f.json`
- Expected SHA-256: `cd12ef9aff0a82aac45d1dc503cce2177f92b52ec95fd229fe4f6d0195ce68d2`
- Scope: exact-SHA full-read card; the bounded independent theorem review and current-use caveat remain separate evidence.
