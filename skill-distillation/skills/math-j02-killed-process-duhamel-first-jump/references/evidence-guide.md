# J02 portable evidence lookup

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

## J02 step 6 local derivation

- Artifact: `audit_current/jk/derivations-j02-step06.md`
- Expected SHA-256: `b01c7485c86d484ed2c0eea43cc71d80e3040123f5f0cf28d8aa844e793e8b00`
- Scope: local measure-definition derivation for a specified reset kernel; see its assumptions before use.

## J02 steps 2–4 local derivation

- Artifact: `audit_current/jk/derivations-j02-steps-02-to-04.md`
- Expected SHA-256: `dbc5f9be9352a433f5c2e802cbd429a531299e926b77954557f8152112b20780`
- Scope: first-jump and semigroup Duhamel steps under their listed Markov, domain, stopping-time and integrability assumptions.

## J02 P-933ed79723ee3851 full-read card

- Artifact: `audit_current/shared-reading/revise_jk/933ed79723ee3851c089f1de04f32502c01bc34665dcef1593e83397c50cd2a1.json`
- Expected SHA-256: `6f97390df875669eba4e4f2ee72fb47552d43cd8b6e9aeb3450392fb36e2db1e`
- Relevance: adjacent diffusion/killing source only.

## J02 P-971c19d9f2f2fc25 full-read card

- Artifact: `audit_current/shared-reading/revise_jk/971c19d9f2f2fc25884618f556faed5a64fff3bb54d391e60ba5ee44f5db1d31.json`
- Expected SHA-256: `e6a0e067bbf8d31bfafe6bff7083ae707bfd60fc23952effa44f61b80d868585`
- Relevance: adjacent killed-process source only.
