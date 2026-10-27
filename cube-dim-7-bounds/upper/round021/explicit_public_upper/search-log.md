# Bounded original-source check for an explicit public upper bound

Date2026-10-27, Asia/Shanghai. Existing academic-literature-search skill applied. Scope: Stein–Strömberg1983 original plus at most two papers directly giving explicit constants; this run needed only the first paper and did not expand to an additional literature sweep. Fixed target model: Lebesgue, centered axis-parallel7cubes, one common radius, allL1, allr>0. No numeric scalar optimization or MATLAB used by this literature/proof subtask. Original sources read-only; own clone only; no commit/push.

Actual executed searches (ordinary web.run, no cited-by database):

1. `Stein Stromberg Behaviour of maximal functions in Rn for large n 1983 pdf`
2. `Stein Stromberg maximal weak 1 1 explicit constant n log n 1983`

Returned counts: not provided. Original paper located in Tsinghua YMSC archive. Secondary search results were navigation only. Direct open of the original PDF succeeded. A guessed ProjectEuclid PDF endpoint was inaccessible via web tool; no claim of having read that endpoint. Python urllib successfully downloaded427259bytes to the outside-repository cache. Full relevant proof, Section2 printedpp260–263, read via web original PDF and then local pypdf extraction; not merely its statement or abstract. Section2 pp264ff strongLp/Euclidean part is not used as an endpoint theorem for cubes.

Original URL: https://archive.ymsc.tsinghua.edu.cn/pacm_download/116/7398-11512_2006_Article_BF02384314.pdf

Exact source SHA-256: `ecb7c3bbbf804545fd68795f5628d242ba8d03f2a3ebeab67806af9c8da6c615`.

Read scope and status: original fulltext of complete relevant theorem/lemma proof; original Theorem1, p260, asserts existence of an unspecified absolute c with c*n*logn. No value of c is printed there. The original proof fixed star dilationn and used geometric radius bins; keeping all its finite factors gives the explicit extractable but unoptimized expression

`[1+((n+2)/n)^n] [1+ceil(logn/log(1+1/n))*((n+1)/n)^(2n)]`.

No numeric conclusion follows merely from the statement's bigO. The draft quantitative refinement now under independent audit is

`C_n(a)=[1+((a+2)/a)^n][1+((a+1)/a)^n(1+nlog a)]`, a>1.

It retains the original finite selection rule with parametera and replaces the geometric-bin bound by an exact finite-sum layer-cake identity. This refinement is **a new derivation in this audit of the old proof method**, not a claimed previously printed explicit numerical theorem. See frozen PROOF-DRAFT.md and source-receipt.json for exact mathematical steps and draft hash. Independent review must complete before adopting a new bound. No additional direct-constant paper was needed/claimed read; unused allowance of two papers was not filled with background references.
