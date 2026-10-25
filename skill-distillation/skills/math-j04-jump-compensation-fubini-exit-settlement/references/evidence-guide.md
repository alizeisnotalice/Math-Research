# J04 portable evidence lookup

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

## J04 P-497c3b054aafd391 full-read card

- Artifact: `audit_current/shared-reading/revise_hi/497c3b054aafd391e9dc7488c5e18fce50ce3dd5498d00807a4626d189fb50e9.json`
- Expected SHA-256: `34f4ab838ff9afaf906ac1caadf3355841a7c909db6d64fa8672d6ecce5a90e0`
- Scope: exact-SHA source reading card; transfer remains limited to the claim conditions recorded in the source and local J04 derivations.
