# C055: sole sparse4 scalar probe completed; fixed condition excluded

The only authorized native call completed successfully at the registered point
`x0=(50,50,50,45,27,0)/50`, `beta=1/80`, `T=32/825`, using Arb512.
The saved interval satisfies `K_T.upper < T/4096 = 1/105600` exactly.
Its classification is `EXCLUDES_ONLY_THIS_FIXED_FRAME_TIME_NECESSARY_CONDITION`.
This excludes the tested necessary condition for this fixed frame and time.
It supplies no whole-cube mT certificate, global g/I certificate, or new C6 bound.

## Authorization and unchanged precheck

Main separately released one native after accepting B055. The release SHA256 is
`70ddaf10a385b49e9cb9c96db88ad02fd5303c2a8d041c77a158b54ffdd37491`;
B055 manifest SHA256 is
`c1171119bbefd07f9120633791d13d5c76dbc43a43034a5481462329f392e1ab`.
The precheck remains frozen at
`e5a61613fa0af7d235c5bd8017e8a9b1171b3577d06509c31ccb75f74088892c`.
Postrun saved-file verification found all18 frozen local files and22 frozen
external files unchanged, as well as the original binary. There was no new
preparation, S2 calculation, compile, source change, extra point/time, or
old053 scalar replay. The prior precheck status files are retained as historical
prelaunch records; this report and full_manifest.json record the new full phase.

## Actual execution and complete accounting

`python3 launch_once.py` was invoked once, after confirming that execution,
token and checkpoint files did not yet exist and no prior native was live.
Tool session70248, initial chunk `f43ca9`, terminal chunk `ea43f8`, native
PID59613 all refer to that same attempt. It exited0 after
91.596721540991s; native kernel time was 91.5061s.
The timeout remained1800s. Darwin maximum child RSS was
225787904bytes. Stderr is empty.

All1920 paired image representatives and13229 cones were processed:

- paired image/cone tuples:25399680;
- logical contributions for each field G and H:609592320;
- G DD5 evaluations:25399680; H DD5 evaluations:25399680;
- total actual DD5 field evaluations:50799360;
- zero-node paired terms:259937, logical terms per field:6238488;
- confluent paired terms:84987, logical terms per field:2039688;
- unique exact cached nodes:303744.

The61 checkpoint records contain60 successive32-image completion records,
ending at1920images, followed by one complete computed record. The exact
counts match the raw native output. There was no skipped term or restart.
The native PID was absent in the postrun process check.
The checkpoint is a persistent progress ledger; it is not a serialized partial
sum. Recovery of this completed job consists of reading its saved output and
receipt. A second native invocation remains forbidden by the attempt marker.

## Saved interval result

The following decimal bounds are rounded outward to50 places from the exact
binary rational endpoints in native_full.stdout.json:

| Quantity | Lower bound | Upper bound |
|---|---:|---:|
| G | 0.00001681892544263902088994249411166539840647790147 | 0.00001681892544263902088994249411166539840647790148 |
| H_T | 0.00000789486793799146978863105303299382098826704766 | 0.00000789486793799146978863105303299382098826704767 |
| K_T=G-H_T | 0.00000892405750464755110131144107867157741821085380 | 0.00000892405750464755110131144107867157741821085381 |
| 1/105600-K_T | 0.00000054563946504941859565825589102539227875884316 | 0.00000054563946504941859565825589102539227875884317 |

The postrun Python verifier used exact Fraction arithmetic for the dyadic
endpoints, verified that K encloses subtraction of the saved G/H intervals,
and proved the saved K upper endpoint is strictly below the exact threshold.
This is a checked interval computation for the registered fixed point/time,
with analytic eligibility supplied by the previously accepted off-true-S2
certificate and F050/F052 dependencies. It is not a new global proof.

## Reproduction and final artifacts

Read-only result verification:

```bash
python3 verify_full_saved.py
```

The verifier reads the saved records and their frozen inputs; it never launches
native, calls the DD engine, recomputes S2, or rebuilds the binary.
full_validation.json is the actual successful output of that verifier.
Raw native stdout, empty stderr, full_execution.json, full_checkpoint.jsonl,
point_comparison.json, release, tool receipt and this report are hashed in
full_manifest.json. Native stdout SHA256 is
`b0a9a0adc8702c58db72b079e6024109eb31315cbd8408298370612b06b9515f`.
No stage, commit or push was performed by this agent. No automatic retry or
additional precision/parameter exploration is authorized.
