# R001 theory review notes

Date: 2026-10-27. These notes accompany `proof.md` and `source_manifest.json` in this directory. The proof is a reconstruction after reading the two identified R027 sources, not a blind audit. Their prior PASS/complete labels are not premises. No commit or push was made.

## Result and mathematical status

The proposed bound is justified for every positive integer n, t>1, fixed axis-parallel cubes with one common radius, all positive scales, general L¹ input, and constant Lebesgue background. Its specialization to n=8 is explicit in proof.md §1. The complete proof retains finite a>2t until after establishing the same full-operator bound for every finite a. I found no unresolved mathematical gap in the reconstructed mechanism.

This author assessment is not a replacement for the assigned fresh project auditor. It makes no numerical claim about the chosen scalar parameter, no optimized value, and no public novelty claim.

## Checks that matter

| Issue | Resolution in proof.md |
|---|---|
| The smallest-radius contributor may have later equal-radius terms | Choose the last selected index among minimum-radius positive contributors, §4. Then every other positive contributor is earlier. |
| The eligibility test does not pay the selected cube itself | Keep the separate self term a^(-n), §§2,4,7. It disappears only in the scalar a→∞ limit. |
| Eligibility selection might skip a larger rejected radius | Prove rejection by the first smaller selection, or by termination if none exists, §3. This covers zero increments and all radius ties. |
| Rejecting centers is insufficient to cover the rejected cubes | Expand b to b+1 and integrate its charge, §3. |
| Cap removal could wrongly discard parts of an increment | Reverse triangle uses b−1; every point of the actual increment is in the larger cap domain, §5. |
| Unstarred increments versus stars | Only the unstarred D_j are disjoint; the proof never requires disjoint star cubes, §§3,5. |
| Potentially wrong shell derivative or factor | Use coordinate-face volume parametrization with determinant 1, §6; compute one captured face per coordinate at the corner to get n/2, §7. |
| Corner maximum might be guessed from a volume derivative | Prove each shell summand is coordinatewise nonincreasing using its face count and side lengths, §6. |
| Boundaries and ties in shell coordinates | Largest-coordinate ties have ambient Lebesgue measure zero; a cap-face equality affects at most one radius since c>1, §6. |
| Finite witnesses might not imply a full maximal bound | Cover each compact subset of the actual open level set with actual bad small cubes; then use inner regularity, §8. |
| Full L¹ and all scales | No density, source, compact-support, product, or radius restriction is made; |f| is used throughout, §8. |
| Sending a to infinity might change the operator or selector | Take the limit only of the already proved scalar constant for the same operator, §9. |
| Constant background | All v_j, cube volumes, dilation Jacobians, shell integrals, and operator denominators use Lebesgue measure. No variable-density extension is claimed, §§1,5,6,8. |
| Positivity of the square root argument | B_n is an explicitly positive inner term plus an outer integral bounded below by (n/2)log t, §9. |
| Dimension one | Empty products are 1; the shell and positivity formulas remain valid, §§6,7,9. |

## Actual local mathematical computation

The default Python 3 and the bundled Python reported `ModuleNotFoundError: No module named 'sympy'` when tested. No dependency was installed. Instead Python 3 standard-library `fractions.Fraction` was actually used for exact rational-function arithmetic over formal indeterminates a,t,ρ. Five cleared-denominator numerators were identically zero:

1. 1+c−α(2+1/t).
2. c−1−(a/t−2)/(a+1).
3. (1+c−α/ρ)/2−α[1+1/(2t)−1/(2ρ)].
4. (1+c−α)/2−α(1+1/t)/2.
5. α/(c−1)−t−2t²/(a−2t).

The actual command exited 0 and printed `exact formal numerator = 0` for all five identities. This exact symbolic check helps detect parameter transcription errors. The geometric and full-operator proof is written explicitly and does not depend on sampled numerical checks. No MATLAB run, floating-point scan, numerical optimizer, or subagent was used by this role.

The following minimal code reproduces the formal check using only Python's standard library. Sparse polynomial keys are powers of (a,t,ρ), and coefficients are exact fractions; no numeric values are substituted for the indeterminates.

```python
from fractions import Fraction

def add(p, q):
    r = p.copy()
    for k, v in q.items():
        r[k] = r.get(k, Fraction(0)) + v
    return {k: v for k, v in r.items() if v}

def mul(p, q):
    r = {}
    for k, v in p.items():
        for l, w in q.items():
            m = tuple(x + y for x, y in zip(k, l))
            r[m] = r.get(m, Fraction(0)) + v*w
    return {k: v for k, v in r.items() if v}

class F:
    def __init__(self, n=0, d=None):
        self.n = n if isinstance(n, dict) else (
            {(0, 0, 0): Fraction(n)} if n else {})
        self.d = d if d is not None else {(0, 0, 0): Fraction(1)}
    def __add__(self, o):
        o = o if isinstance(o, F) else F(o)
        return F(add(mul(self.n, o.d), mul(o.n, self.d)),
                 mul(self.d, o.d))
    __radd__ = __add__
    def __neg__(self):
        return F({k: -v for k, v in self.n.items()}, self.d)
    def __sub__(self, o):
        return self + (-o if isinstance(o, F) else F(-o))
    def __rsub__(self, o):
        return F(o) + -self
    def __mul__(self, o):
        o = o if isinstance(o, F) else F(o)
        return F(mul(self.n, o.n), mul(self.d, o.d))
    __rmul__ = __mul__
    def __truediv__(self, o):
        o = o if isinstance(o, F) else F(o)
        return F(mul(self.n, o.d), mul(self.d, o.n))
    def __rtruediv__(self, o):
        return F(o) / self

A = F({(1, 0, 0): Fraction(1)})
T = F({(0, 1, 0): Fraction(1)})
R = F({(0, 0, 1): Fraction(1)})
alpha = A / (A + 1)
c = (A*(1 + 1/T) - 1) / (A + 1)
checks = {
    'c_plus_1': 1 + c - alpha*(2 + 1/T),
    'c_minus_1': c - 1 - (A/T - 2)/(A + 1),
    'corner_kernel': (1 + c - alpha/R)/2
                     - alpha*(1 + 1/(2*T) - 1/(2*R)),
    'corner_inner': (1 + c - alpha)/2 - alpha*(1 + 1/T)/2,
    'corner_transition_gap': alpha/(c - 1) - T - 2*T*T/(A - 2*T),
}
for name, value in checks.items():
    assert value.n == {}, (name, value.n)
    print(name + ': exact formal numerator = 0')
```

## Clone and filesystem provenance

The independent local clone is `/Users/zhengzhihao/Desktop/T/Math-Research-cube8-sol-theory`, created with

```text
git clone --no-hardlinks --no-checkout /Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7 /Users/zhengzhihao/Desktop/T/Math-Research-cube8-sol-theory
git sparse-checkout init --no-cone
git sparse-checkout set /AGENTS.md /cube-dim-8-bounds/upper/round001/theory/
git checkout
```

All four tool-returned exit codes were 0. However, `git checkout` printed twice:

```text
error: invalid object 100644 48d436008eaeecca94834fc826350c6e94e26490 for 'centered-cube-upper-bound/adversarial_shell_pressure_20261007/adversarial_shell_search.py'
```

It also printed that the branch was up to date with origin/main. This is an observed checkout diagnostic despite exit code 0, not a verified clean full checkout. Root `AGENTS.md` was materialized and read, and only the assigned new n8 theory directory was authored. The base HEAD resolves to `79f8acb20c418e60c5196a3c830103253e34c900`. The unrelated missing object was not repaired or modified. The original project was only read.

The three deliverables are exactly `proof.md`, `source_manifest.json`, and `review_notes.md`. Source hashes and absolute paths are in the manifest. No cached experiments, temporary scripts, unrelated files, stage changes, commits, or pushes belong to this role's output.
