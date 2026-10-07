#!/usr/bin/env python3
"""Build the single source of truth: a unified assertion register.

Resolves the 20-vs-28 counting inconsistency by enumerating every verification
object from the three batch artifacts plus the two sharp-constant addenda,
assigning each a unique ID, and classifying along TWO independent dimensions:

  verification_status (mutually exclusive, strongest achieved):
    independent_proof   -- this audit produced a complete argument incl. branches,
                           degeneracy, attainability
    conditional_proof   -- complete argument under premises taken as given
                           (standard theorem / stated hypothesis)
    instance_check_only -- numeric or finite-instance verification
    pending             -- insufficient evidence or partially checked
    refuted             -- found false (none so far)

  evidence_basis (possibly several):
    own_derivation | source_comparison | standard_theorem_citation |
    numeric_experiment | literature_check

Counts are reported at three levels: original assertions, sub-assertions
(when one original object splits into independently-graded parts), and tier
membership. The register also records where each earlier ledger entry was
dropped (f07 upper bound, e01 CVaR) so the correction is traceable.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
from collections import Counter

R = pathlib.Path(__file__).resolve().parents[1]

# ---------------------------------------------------------------- registry
# id: (skill, statement, status, bases, parent_original, notes)
A: list[dict] = []

def add(aid, skill, stmt, status, bases, parent, note=""):
    A.append({"id": aid, "skill": skill, "statement": stmt,
              "verification_status": status, "evidence_basis": bases,
              "parent_original": parent, "note": note})

# ---- B1 case numerics (audit_r8/numeric_assertion_verification.json) ----
add("J01-CASE-01", "math-j01-inhomogeneous-jump-generator",
    "P(X_1=0)=exp(-int_0^1(1+s)ds)=e^{-3/2}=0.22313 for rate lambda_01(t)=1+t",
    "instance_check_only", ["numeric_experiment", "standard_theorem_citation"],
    "J01", "closed form exact; time-varying-rate Poisson formula is a standard result not independently proved")
add("J02-CASE-01", "math-j02-killed-process-duhamel-first-jump",
    "survival e^-2=0.135335, cemetery 1-e^-2=0.864665, density 2e^{-2t} at rate 2, T=1",
    "instance_check_only", ["numeric_experiment", "standard_theorem_citation"], "J02",
    "three quantities mutually consistent")
add("I02-CASE-01", "math-i02-filtration-projections-compensators",
    "EN_1=EA_1=1-e^{-1}=0.63212 for A_t=min(t,gamma)",
    "instance_check_only", ["numeric_experiment", "standard_theorem_citation"], "I02",
    "the structural claim (compensator changes once F_0 reveals gamma) NOT verified")
add("I03-CASE-01", "math-i03-fourth-moment-suffix-maximal-comparison",
    "E|M_2|^4=8, E[max_{k<=2}|M_k|]^4=8.5, 8.5<=(4/3)^4*8",
    "independent_proof", ["own_derivation", "numeric_experiment"], "I03",
    "four-path complete enumeration -- exact")
add("H04-CASE-01", "math-h04-exact-rate-function-variance-strata",
    "I(3/4)=0.130812 for Bernoulli(1/2); endpoints I(0)=-log(1-p), I(1)=-log p",
    "conditional_proof", ["own_derivation", "standard_theorem_citation"], "H04",
    "substitution and endpoint extension algebra complete; Cramer's theorem itself standard, not proved here")
add("H05-CASE-01", "math-h05-chernoff-hoeffding-bernstein",
    "E e^X = cosh(1)=1.5431 <= e^{1/2}=1.6487 (strict) for fair Rademacher",
    "instance_check_only", ["numeric_experiment", "standard_theorem_citation"], "H05",
    "strictness at this point is a computation; Hoeffding lemma general form not proved here")
add("K04-CASE-01", "math-k04-anti-concentration-empirical-process",
    "Var(|G|)=1-2/pi; Q=2Phi(0.1)-1=0.07966; 0.0478<=Q<=0.574",
    "instance_check_only", ["numeric_experiment", "standard_theorem_citation"], "K04",
    "Theorem 2 bound itself not verified")
add("A01-CASE-01", "math-a01-fractional-normalization",
    "on X=[0,1], max x/(1+x)=1/2 attained at x=1",
    "independent_proof", ["own_derivation"], "A01", "monotonicity + endpoint")

# ---- sharp constants (B1 PART_2 + B2 additional_results) ----
add("A02-SHARP-01", "math-a02-outer-relaxation-budget",
    "(x-kappa)_+ <= x^2/(4kappa) for all real x, constant 1/4 not reducible; equality iff x=0 or x=2kappa",
    "independent_proof", ["own_derivation", "numeric_experiment"], "A02",
    "algebraic identity + global critical-point analysis + 4000-point exact rational search")
add("F07-SHARP-01a", "math-f07-marcinkiewicz-zygmund-vector-extension",
    "MZ constant LOWER bound >= 1 (equality attained)",
    "conditional_proof", ["own_derivation"], "F07",
    "equality example: m=1, T=identity, p=q=2, f_k co-directed -> Minkowski equality; PREVIOUS LEDGER DROPPED the upper-bound half")
add("F07-SHARP-01b", "math-f07-marcinkiewicz-zygmund-vector-extension",
    "MZ constant UPPER bound <= 1",
    "instance_check_only", ["numeric_experiment", "standard_theorem_citation"], "F07",
    "300 random positive multi-linear families, max ratio 1.000000; general duality proof NOT rewritten. PREVIOUS LEDGER DROPPED this half entirely")
add("B02-SHARP-01", "math-b02-semi-infinite-exchange",
    "min x s.t. x>=t for all t in [0,1] has optimum 1",
    "independent_proof", ["own_derivation"], "B02", "")
add("E01-SHARP-01", "math-e01-occupation-measure-flow-lp",
    "Cantelli radius r=sqrt(1/eps-1) cannot be reduced",
    "conditional_proof", ["own_derivation", "numeric_experiment"], "E01",
    "two-point construction: mean 0, variance exactly 1, P(X>=r)=eps exact; Cantelli inequality itself standard, not proved here")

# ---- B2 batch (audit_r9/batch_b2_verification_r9.json) ----
add("E01-CTX-01", "math-e01-occupation-measure-flow-lp",
    "VP validity range 'unimodal and 0<eps<=1/6' is the exact threshold",
    "pending", ["literature_check"], "E01",
    "RECLASSIFIED from R: literature confirmation is the BASIS, but as a verification status this is a cited fact, not my proof; threshold sqrt(8/3)<=>1/6 exact; my initial 4/81 suspicion refuted by literature")
add("E01-BND-01", "math-e01-occupation-measure-flow-lp",
    "VP extremizer (mass 1-4/(3lam^2) at mean + centered uniform) attains 4/(9lam^2)",
    "conditional_proof", ["own_derivation", "numeric_experiment"], "E01",
    "variance and tail computed exactly at lambda=1.7/2.0/3.0; VP theorem premise standard")
add("D03-CTX-01", "math-d03-parameterized-free-energy",
    "source Theorem 1.1 (paper.md:93,198) really writes +infty outside the OPEN interval",
    "instance_check_only", ["source_comparison"], "D03",
    "direct textual check of the frozen conversion; RECLASSIFIED: this is a source-comparison fact; PREVIOUS LEDGER folded it into prose rather than listing it")
add("D03-IDEN-01", "math-d03-parameterized-free-energy",
    "psi*(+/-||t||_1)=m ln2; |alpha|>||t||_1 gives +inf; infinite support endpoint = inf; t=0 gives psi*(0)=0",
    "conditional_proof", ["own_derivation", "numeric_experiment"], "D03",
    "asymptotic expansion complete to 4.5e-14 with stable logcosh; attainability INSIDE the open interval not verified")
add("E01-CEX-01", "math-e01-occupation-measure-flow-lp",
    "CVaR atom-inconsistency: Eq(3)=2 vs Eq(4)=Eq(5)=1 for V=1 w.p. 1/2, eps=1/4",
    "conditional_proof", ["own_derivation", "numeric_experiment"], "E01",
    "recomputed exactly with rationals; PREVIOUS LEDGER DROPPED this entry entirely")
add("E01-CEX-02", "math-e01-occupation-measure-flow-lp",
    "Theorem 5.5 literal fixed-scalar no-gap is false: fixed-t optimum 2/3 vs measure LP 1",
    "conditional_proof", ["own_derivation", "numeric_experiment"], "E01",
    "every pivotal step independently re-derived (ES 2q/3 attainable, domination <=1, chain-rule residuals<1e-14); scope limited to the literal scalar-t* semantics; A7 uses C1 polynomial density, not source [20]")
add("F06-BND-01", "math-f06-menshov-rademacher-maximal-partial-sums",
    "||max_{k<=N}|sum f_j|||_2^2 <= (m+1)^2 sum||f_j||^2, m=ceil(log2 N)",
    "instance_check_only", ["numeric_experiment"], "F06",
    "1400 orthogonal trials, worst ratio 1.76 vs bound 49; NOT a proof; skill itself labels it 'sufficient, not claimed optimal'")
add("E04-CASE-01", "math-e04-snell-envelope-optimal-stopping",
    "two-period positive case S_0=1, tau*=1; negative case S_0=0 beats forced stop -1",
    "independent_proof", ["own_derivation"], "E04", "finite-state enumeration")
add("N01-CASE-01", "math-n01-john-loewner-ellipsoid",
    "||x||_2 <= sqrt(n)||x||_inf <= sqrt(n) on [-1,1]^n, equality at vertices",
    "independent_proof", ["own_derivation", "numeric_experiment"], "N01", "")

# ---- B3 batch (audit_r10/batch_b3_verification_r10.json) ----
add("D02-IDEN-01", "math-d02-gibbs-entropy-duality",
    "KL(P||Q_V)=KL(P||Q)-E_P V+log Z; D(mu_C||mu)=ln(1/mu(C))",
    "instance_check_only", ["numeric_experiment", "standard_theorem_citation"], "D02",
    "5 instances at 1e-12 + exact special case; general algebraic proof not rewritten; P<<Q, integrability, extended-real conventions NOT independently audited")
add("D04-IDEN-01", "math-d04-entropy-chain-rule",
    "D(P_XY||Q_XY)=D(P_X||Q_X)+E_{P_X}D(P_{Y|x}||Q_{Y|x})",
    "instance_check_only", ["numeric_experiment", "standard_theorem_citation"], "D04",
    "4 instances at 1e-12; RN-decomposition proof not rewritten; zero-mass handling not audited")
add("M01-BND-01a", "math-m01-convex-hull-extreme-points",
    "support at most d+1 (LP basic-solution theory)",
    "pending", ["standard_theorem_citation"], "M01",
    "cited, NOT independently proved; register as pending until the LP argument is written out or the citation is pinned")
add("M01-BND-01b", "math-m01-convex-hull-extreme-points",
    "sharpness: instances requiring d+1 exist (150-trial max support = d+1 for d=1..8)",
    "instance_check_only", ["numeric_experiment"], "M01",
    "numeric observation; 'ALL instances need d+1' is false and not claimed")
add("N04-IDEN-01", "math-n04-summable-error-entropy-descent",
    "sum d_k <= H_0 - H_inf + sum eps_k^+",
    "independent_proof", ["own_derivation"], "N04",
    "one-line telescoping argument, upgraded from numeric in R11")
add("N04-IDEN-02", "math-n04-summable-error-entropy-descent",
    "necessity of SKILL item-3 warning on signed errors",
    "pending", ["own_derivation"], "N04",
    "my check has NO discriminating power (sum eps <= sum eps^+ identically); the warning is conservative-correct but its necessity is unproven")
add("C01-CASE-01", "math-c01-semidefinite-lifting",
    "two-ball QCQP: x*=(-1,0) feasible, objective exactly -0.54",
    "instance_check_only", ["own_derivation", "numeric_experiment"], "C01",
    "optimality of x*, Beta exactness theorem, l1^T W l2=0 sufficiency NOT verified")
add("C02-IDEN-01a", "math-c02-completely-positive-copositive",
    "CP subset COP and CP subset COP* (one-directional inclusions)",
    "independent_proof", ["own_derivation", "numeric_experiment"], "C02",
    "one-line algebra each: (v^T x)^2>=0 for x,v>=0; <A,sum w w^T>=sum w^T A w >=0 for A in COP, w>=0")
add("C02-IDEN-01b", "math-c02-completely-positive-copositive",
    "CP = COP* (equality, reverse inclusion COP* subset CP)",
    "pending", ["standard_theorem_citation"], "C02",
    "Diananda-type deep result, cited not verified")
add("C02-IDEN-01c", "math-c02-completely-positive-copositive",
    "PSD basis size binom(n+d,d)",
    "independent_proof", ["own_derivation"], "C02", "combinatorial dimension of the symmetric-polynomial basis")
add("H02-IDEN-01", "math-h02-exponential-tilting-saddlepoint",
    "Gaussian saddlepoint prefactor equals the Mills leading term exp(-x^2/2)/(x sqrt(2pi))",
    "instance_check_only", ["own_derivation", "numeric_experiment"], "H02",
    "1e-15 agreement, relative error decreasing; general-K prefactor, error bounds, Lugannani-Rice branch NOT verified (skill has its own guard)")
add("B03-IDEN-01", "math-b03-generalized-moments",
    "sup E[X|X>=t] = mu + sigma^2/(mu-t) for t<mu, sigma>0; sup=inf for t>=mu; NOT attained",
    "independent_proof", ["own_derivation", "numeric_experiment"], "B03",
    "two-component variance decomposition + simultaneous boundary equality; branches t<mu/=mu/>mu, degeneracy sigma=0, p=1, attainability (X=t atom belongs to Xi) all covered; upgraded in R11")

# ------------------------------------------------------------- reporting
by_status = Counter(a["verification_status"] for a in A)
parents = sorted({a["parent_original"] for a in A})
# count parents properly: distinct (skill,parent) pairs where parent==skill code
orig_ids = sorted({a["parent_original"] for a in A})
print("=" * 72)
print("统计口径（三层，互不混用）")
print("=" * 72)
print(f"  原始断言（distinct parent 对象）: {len(orig_ids)}")
print(f"     {orig_ids}")
print(f"  子断言（本登记簿条目）        : {len(A)}")
print(f"  验证状态分布                  : {dict(by_status)}")
print()
print("此前两个错误计数的解释：")
print("  『20 条』= B1案例8 + sharp4 + B3推导8，遗漏了 B2 的 9 条中未重复的部分")
print("  『28 条』= ledger 条目数，但丢了 F07 上界、E01-CVaR、D03-指控 三条，")
print("             又把 VP极值、n04警告 拆出单列 —— 拆分无映射、丢失无记录")
print()
print("=" * 72)
print("拆分映射（原始 → 子断言）")
print("=" * 72)
from collections import defaultdict
m = defaultdict(list)
for a in A:
    m[a["parent_original"]].append(a["id"])
for p in sorted(m):
    ids = m[p]
    if len(ids) > 1:
        print(f"  {p} → {len(ids)} 条: {', '.join(ids)}")
    else:
        print(f"  {p} → 1 条: {ids[0]}")

out = R / "audit_r12" / "assertion_register.json"
pass  # dir exists
payload = {
    "schema": "assertion-register-r12",
    "generated": "2026-10-06",
    "counting_units": {
        "original_assertion": "one verification object as executed in B1/B2/B3 (distinct parent_original)",
        "sub_assertion": "register entry; an original splits when its parts have DIFFERENT verification strengths",
        "tier_count": "membership of sub-assertions in a status class; tiers are the status values themselves, mutually exclusive by construction",
        "original_count": len(orig_ids), "sub_count": len(A),
        "status_distribution": dict(by_status),
        "correction_history": {
            "claim_of_20": "B1(8 cases)+sharp(4)+B3(8) -- omitted 5 distinct B2 objects",
            "ledger_28": "dropped F07 upper bound, E01-CVaR, D03-CTX; split VP-extremizer and N04-warning without a recorded mapping"}},
    "dimension_definitions": {
        "verification_status": {
            "independent_proof": "complete argument incl. branches, degeneracy, attainability produced by this audit",
            "conditional_proof": "complete argument under premises taken as given (standard theorem / stated hypothesis)",
            "instance_check_only": "numeric or finite-instance verification only",
            "pending": "insufficient evidence, partial check, or unresolved necessity",
            "refuted": "found false (none so far)",
            "mutually_exclusive": True,
            "priority_rule": "a claim whose parts differ in strength is SPLIT into separate sub-assertions rather than assigned a compound label"},
        "evidence_basis": {
            "values": ["own_derivation", "source_comparison", "standard_theorem_citation",
                        "numeric_experiment", "literature_check"],
            "mutually_exclusive": False,
            "note": "independent dimension: HOW the verification was performed, not how strong it is"}},
    "assertions": A}
out.write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
print(f"\nwritten {out}")
print("sha256", hashlib.sha256(out.read_bytes()).hexdigest())