# J03 portable evidence lookup

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

## J03 steps 1–2 local derivation

- Artifact: `audit_current/jk/derivations-j03-steps-01-to-02.md`
- Expected SHA-256: `cd65b6f24f9f17a3bac191533fdc32273b7e603a891056c18b5132c52737f2dc`
- Scope: right-continuous exit-time/overshoot examples and exit-measure mass checks; not a general theorem.

## J03 P-0178f07ee92816cf full-read card

- Artifact: `audit_current/shared-reading/revise_jk/0178f07ee92816cf4cba7861597f548038409453a8d2af0fbfb8b4490a9e19a0.json`
- Expected SHA-256: `90e98cbc2a3b153d1541ddddbb4d381218eb99bf265fda434545b8a5f0c36306`
- Scope: finite-state projection scheme for minimal time-homogeneous countable CTMCs under the recorded stability/integrability assumptions.

## J03 P-93cb5642ab2d6b8f full-read card

- Artifact: `audit_current/shared-reading/revise_jk/93cb5642ab2d6b8f1ee2d28dd699f0db690415fbf92dade9cbc38df6ff0b9d69.json`
- Expected SHA-256: `27c85dd8361865a3e7b17b00806ef1c62dfc4cf7a5805c64e6a3104ebcf1e387`
- Scope: telegraph-model exit distributions only; adjacent evidence.
