# Historical Skill filename repair — 2026-10-11

A recursive Skill scan found 12 archived copies named `SKILL.md` below the J01–J05 and K01–K05 reference-history folders. Their basenames were changed to `SKILL.original.md` in the local bundle and public Skill tree. All 12 archived file contents and historical SHA-256 values are unchanged. This rename does not alter any current top-level Skill entrypoint, claim ledger, mathematical statement, or evidence grade.

The alias table `historical-skill-path-aliases.json` maps each prior path to its new path. Immutable historical audit records and the two archived manifests that mention old paths remain byte-for-byte unchanged. No active Markdown links targeted these archived copies.

Application evidence remains separated by provenance: the accepted r18 matrix contains 128 current slots and 100 exact target-claim bindings; 24 recent examples are explicitly manual and unblinded, with API/native execution recorded as zero for those 24 only; the older 104 records retain unknown/null external-call fields. A zero in the live-progress/index native-task evidence counter means no native-task evidence is recorded in that counter. It does not mean 128 executions were attempted and failed or that the Skills are ineffective. Installer hashes and resolver path/hash behavior are operational checks, separate from task-level mathematical behavior.

The prior ZIP and its acceptance records remain historical and unchanged. GitHub PR #1 was open/draft at the pre-repair commit; this repair is carried in a follow-up commit and updated PR, without merging the PR.
