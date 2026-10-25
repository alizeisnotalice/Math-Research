# Lower-bound archive whitelist: rounds 004--009

This is a retention whitelist, not a deletion instruction.  Everything outside
it is historical working material; in particular it must not be used as fresh
execution evidence.  A replay needs the listed round source plus
`lower/certify_rectangular.py` and, where stated, the locked complete C source
`lower_round003/a2_equals_a3_results.json`.

| round | conclusion | retain for replay | certificate/result hashes |
| --- | --- | --- | --- |
| 004 | Eight common-last-four rational offsets failed; closest `-1/256` remained below C. | `registration.json`, `common_tail_batch.py`, `common_tail_batch_results.json`, `common_tail_followup.md` | source `ac352dce7edde23f7fc3d369520c1e1778b0aed124e84c8e2dd2c9c1d2e82683`; result `371ff6ac7b91d613878f99c79160081415562b1eaad4db5b84eb559e22f147c9` |
| 005 | Full rational `l<=5` finite family did not improve C; five length-5 boxes were dominated. | `registration.json`, `rational_l5.py`, `rational_l5_results.json`, `report.md` | source `e2af3878efd7ffd5e8af232b4057490660bac5e59b71fb8c89570c30bad8fb98`; result `47c85d2c3950dc48734e83e9bb84f56790af9bbc90c24fbec27a58bca7907b7e` |
| 006 | Three alternating-first-axis **weight** finite periodic tests all failed; no transfer claim. | `registration.json`, `alternating_first_axis.py`, the three `alternating_theta_*_results.json`, `report.md`, `execution_log.md` | source `4154101294bd381dc98fc401f28f157bf27c1e6103de7d1919b8b827bac9bb15`; result hashes are recorded in `report.md` |
| 007 | The four nearest endpoint-order intervals on each side of C in the common-tail line have exact `P-C<=0` Bernstein exclusions. | `registration.json`, `common_tail_event_intervals.py`, `exactC_event_interval_0_results.json` through `exactC_event_interval_7_results.json`, `report.md`, `diagnostic_invalid_prepatch_interval.md` | source `7b9b688982af07cfc784c5ab8d16795697e84be4382913a17faac19ce17ac383`; eight final hashes are listed in `report.md` |
| 008 | Symmetric two-group tail split at three positive offsets failed against complete C. | `registration.json`, `two_group_tail.py`, three `two_group_delta_*_results.json`, `report.md` | source `32aad097557126865bdae072483e645c401371c5f7ca5caf75559125e1835dd8`; result hashes in `report.md` |
| 009 | Three equal-mass alternating-first-axis-gap periodic tests failed against complete C. | `registration.json`, `alternating_spacing.py`, three `alternating_spacing_delta_*_results.json`, `report.md` | source `1efeef8cec90cab8372f53f4ad90dd5b2582d86325eb12a67f309ac434558cc6`; result hashes in `report.md` |

## Round 007 compaction rule

Only the eight `exactC_event_interval_*_results.json` finalizations are valid
round-007 certificates.  Their final files embed the six exact dual-union
sample certificates required to replay interpolation, so the 48 standalone
corrected sample JSON files are not retained by this whitelist.  The old
prepatch sample set crossed omitted events and is replaced by
`diagnostic_invalid_prepatch_interval.md`; it is expressly not evidence.
Likewise, the provisional comparison to C's isolating upper endpoint is
withdrawn and is not whitelisted.
