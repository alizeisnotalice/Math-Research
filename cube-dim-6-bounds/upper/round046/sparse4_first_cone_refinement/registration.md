# A046 sparse4 first-cone exact refinement

Task started UTC 2026-10-25T21:28:35Z; hard stop 21:48:35Z. Registered UTC 2026-10-25T21:29:56.556564+00:00.
Only this directory writable. Read-only B044 report/manifest/frame hashes matched parent values. No delegation, MATLAB, kernel/native/G/H/K/p, optimizer, full fan, old2160 baseline, commit or push.

Fixed scope: old parent cone1 only, beta=1/80 true152-column frame. Independently test 120 new forms on its six rays (720 dots), identify exactly crossings (2,10,26,50,82), then handle all32 closed sign patterns on the slice sum u=1,u>=0. Eleven halfspaces per cell. Enumerate at most14784 active-set solves; cache the sign-independent boundary systems so at most462 distinct 6x6 rational systems are needed, with exact feasibility rechecked for all32 patterns.

Use exact FLINT rational vertices/ranks; recursive pulling from the actual11 halfspaces and one global vertex order. No float hull as proof. Check each output simplex lies in parent and its cell, all152 forms have constant weak sign, exact determinant-volume sum is1 relative to parent, primitive physical rays/rates use true152 columns. Drop lower-rank sign cells only after retaining their vertices and proving they are covered by closures of full cells. Global coverage is the halfspace arrangement partition plus pulling induction, not merely a volume check.

Limits:32 patterns,14784 active sets,100000 output full-dimensional simplices,1200-second wall limit from task start. If limits or a correctness check fails, stop and preserve actual gap; no expansion. Outputs planned: check_refinement.py, local_fan.json, results.json, report.md/theorem.md, run receipt and manifest. This certifies only one parent cone, no global heat kernel or Cn.
