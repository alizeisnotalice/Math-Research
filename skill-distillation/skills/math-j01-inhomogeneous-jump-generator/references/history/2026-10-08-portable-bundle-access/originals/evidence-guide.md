# J01 portable evidence lookup

Evidence artifacts are kept in the accepted verification bundle rather than copied into an installed Skill. Their paths are package-relative. From the bundle root (the directory containing `workspace_revised/`, `audit_current/`, and `evidence/`), resolve and hash a specific artifact with:

```sh
python3 audit_current/delivery/resolve_evidence.py --artifact-path audit_current/shared-reading/audit_jk/500a71307a38c036d97bca3c6b958c710c22abd9ca1479c201a1ed17713b424f.json
```

If the resolver or bundle root is unavailable, mark the external artifact unresolved; do not follow a relative path from the installed Skill directory or treat an omitted artifact as included.

## J01 E0785 full-read card

- Artifact: `audit_current/shared-reading/audit_jk/500a71307a38c036d97bca3c6b958c710c22abd9ca1479c201a1ed17713b424f.json`
- Expected SHA-256: `cd12ef9aff0a82aac45d1dc503cce2177f92b52ec95fd229fe4f6d0195ce68d2`
- Scope: exact-SHA full-read card; the bounded independent theorem review and current-use caveat remain separate evidence.
