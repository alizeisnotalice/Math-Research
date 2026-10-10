# J02 portable evidence lookup

Evidence artifacts remain package-relative in the accepted verification bundle. From the bundle root (containing `workspace_revised/`, `audit_current/`, and `evidence/`), run `python3 audit_current/delivery/resolve_evidence.py --artifact-path <artifact-path>` to confirm a file exists and verify its SHA-256. If that bundle root is unavailable, report the artifact as unresolved; the installed Skill does not embed these audit files.

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
