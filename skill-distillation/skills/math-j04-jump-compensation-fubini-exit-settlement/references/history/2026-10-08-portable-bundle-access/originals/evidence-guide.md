# J04 portable evidence lookup

Evidence artifacts remain package-relative in the accepted verification bundle. From the bundle root (containing `workspace_revised/`, `audit_current/`, and `evidence/`), run `python3 audit_current/delivery/resolve_evidence.py --artifact-path <artifact-path>` to confirm a file exists and verify its SHA-256. If that bundle root is unavailable, report the artifact as unresolved; the installed Skill does not embed these audit files.

## J04 P-497c3b054aafd391 full-read card

- Artifact: `audit_current/shared-reading/revise_hi/497c3b054aafd391e9dc7488c5e18fce50ce3dd5498d00807a4626d189fb50e9.json`
- Expected SHA-256: `34f4ab838ff9afaf906ac1caadf3355841a7c909db6d64fa8672d6ecce5a90e0`
- Scope: exact-SHA source reading card; transfer remains limited to the claim conditions recorded in the source and local J04 derivations.
