# B023 stage1: true denominator degree417, with an exact Green coefficient

**The exact partial-fraction compression completed and stopped.** At the new
point x0=(50,50,50,45,27,0)/50 for the original body32 Cauchy frame, the true
reduced density denominator in z=T^2 has degree417, with359 distinct positive
pole parameters. This is below the earlier512-degree gate. The raw889-degree
candidate was not the reduced degree and would have incorrectly rejected the
route if treated as one.

No global numerator, common-denominator expansion, polynomial gcd, Sturm
sequence or root isolation was formed. No positive time was sampled. The
work is explicitly **one full symbolic kernel computation**, zero numerical
kernel calls, zero new positive times, and eight formal scalar identity
cases. It is not a global cube certificate, all-time inequality, or C6 bound.
Any next rational-function or all-time stage requires separate authorization.

## Exact reduction and retained tables

|Quantity|Exact count|
|---|---:|
|Frozen input tuples processed|40320|
|Represented full-group logical terms|967680|
|Represented point-orbit/cone terms|80640|
|Weighted local partial contributions|241920|
|Signed `(b,order)` slots touched|1148|
|Nonzero signed coefficients A|736|
|Nonzero positive-fold coefficients C|417|
|Nonzero ordinary z-principal coefficients|417|
|Raw squared-pole factors / degree|792 /889|
|True reduced positive factors / degree|359 /417|

The true exponent histogram is311 factors of order1,38 of order2,10 of
order3. In particular the raw order4 factor does not survive with that order.
The complete three nonzero tables and true factor exponents are retained in
`partial_tables.json`; no dense large numerator is hidden in that file.
Each ordinary real principal term has the form `d_(q,j)/(z+q)^j` with q>0.
All coefficients and q are exact rationals.

The zero-node channels were accumulated before taking real parts. They are

```
A_(0,1) = 238211881100772294265625 / 338966153899070810112,
A_(0,2) = A_(0,3) = A_(0,4) = A_(0,5) = A_(0,6) = 0.
```

Thus all even zero-pole channels are exactly zero, as required. The nonzero
order1 channel is purely imaginary after multiplication by i^(1-6); its real
part is identically0. It was not discarded from the input or overwritten.
The chosen global-negation representative is not itself a conjugate-symmetric
complex sum; only its real part is the density. This explains why a nonzero
imaginary zero-pole coefficient can remain without a physical singularity.

## Complete identity and genuine reduced degree

For each tuple, group equal signed nodes b with multiplicities m_b. In exact
fmpq arithmetic, expand

`prod_(c!=b)(b-c+u)^(-m_c)`

through degree m_b-1. Its coefficient a_k contributes
`i^(m_b-k-6)*a_k/(T-i*b)^(m_b-k)` to the complex partial fraction. This rule
retains zero b and all repeated nodes. The frozen tuple weight multiplies
each contribution exactly once; it already includes stabilizer12 and both
global signs. The tuple weights sum to the independently known full alpha
checksum520066141/48769479475200, but that checksum is only accounting and
is not used as a substitute for the partial identity.

After complete accumulation, the real fold is explicitly

`C_(b,j)=A_(b,j)+(-1)^j*A_(-b,j)`, b>0.

The real numerator over `(z+b^2)^j` is

`N_(b,j)(z)=sum_r (-1)^(j-r-3) binom(j,2r)b^(j-2r)z^r`,

where0<=r<=floor(j/2). At its pole,

`N_(b,j)(-b^2)=(-1)^(j+1)2^(j-1)b^j !=0`.

Therefore the highest nonzero folded coefficient at each b determines the
actual pole order. Lower orders cannot cancel it, and other b values have
coprime distinct factors. `checked_results.json` contains every highest-pole
witness. The exact degree417 follows locally; it does not rely on counting
nonzero signed A before folding or on estimating cancellation numerically.

`theorem.md` proves the complex identity by matching principal parts and
vanishing at infinity, the real folding parity, the local numerator formula,
the reduced-denominator statement and the Green integration formula. The
proof includes all zero/confluent cases, not just the eight test examples.

## New exact Green coefficient from this same symbolic pass

The complete rational

`Gamma_G=sum_(b>0) C_(b,1)`

and `G(x0)*pi^5=Gamma_G/128` are stored exactly in `results.json`. Their
numerators and denominators are retained in full, not decimal fits. This is
a new full symbolic kernel result at x0, not a value inherited from B015 or
a claim of zero full-kernel work.

Only nonzero-pole order1 contributes to the real time integral:

`Re[-i*A/(T-i*b)]=A*b/(T^2+b^2)`,

whose integral is `(pi/2)sign(b)A`. For j>=2, integrating the complex power
and multiplying by i^(j-6) gives a purely imaginary endpoint, so its real
integral is0. All even zero poles have already cancelled; odd zero channels
are identically0 in the real function. These steps integrate in dT, not dz.
The outside normalization `(2*pi)^-6` then gives

`G(x0)=Gamma_G/(128*pi^5)`.

The saved-data checker independently obtains the SAME exact Gamma_G from
three expressions: the signed order1 table, the positive folded order1
table, and the ordinary z-principal table integrated by the rational formula

`integral_0^infinity (T^2+b^2)^(-j)dT
 =(pi/2)binom(2j-2,j-1)/(4^(j-1)b^(2j-1))`.

No numerical time evaluation or additional full tuple traversal is used for
this check. The finite real partial terms are absolutely integrable after
the zero-pole cancellation, justifying termwise integration.

## What was actually checked

Before the full symbolic pass, the eight registered FORMAL scalar cases
were run once. They include distinct same-sign and mixed-sign nodes, one and
two zeros, double/triple/quadruple repetitions, and six zeros. For each case,
an independent degree-at-most6 complex polynomial multiplication verifies
that all proposed partial fractions sum to the reciprocal product exactly.
An independent real polynomial check verifies the numerator
`z^3-e2*z^2+e4*z-e6` against the folded expression. No T was sampled.

The full aggregation then processed all40320 frozen tuples once. Logical
terms with zero nodes total88320 and with repeated nodes total201408; all
were retained. Its accounting totals, final even-zero cancellations and
local denominator exponents passed exact checks.

Afterward, `check.py saved` used Python Fraction independently of the fmpq
aggregator to rebuild the positive fold, the shifted real principal table,
all highest-pole coefficients, the denominator degree, and the three Green
coefficient expressions. It verified the pinned input and output hashes.
This saved-table checker does **not** claim to have independently re-summed
all40320 tuples; a second full symbolic pass was not authorized or performed.
The identity proof, exact primary aggregation and its scalar tests remain
the basis of that full sum, pending the main agent's stage audit.

## Frozen source, actual execution and stop receipt

Source directory is the read-only C022
`/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round022/reduced_point_preflight/`.
Its archive is `merged_tuples.json.gz`, SHA
`562e3703286fc9a9712efdad30bbebb71433c4f3f15f867c985ee41b80c6f719`;
decompressed content SHA
`7deb058449e11fa74d7e1a4658311cb738a38fbbb57ab03118b494138cd39390`.
The source report SHA is
`fad56fada09c08e96311fa7bf25f8d5326ca67947415688febb964d47bac72c5`.
The original geometry, point-orbit enumeration and tuple merger were NOT
rerun. The archive was read losslessly and all stored multiplicities used.

Registration was PID37533, tool chunkf9ed15, exit0. The input-manifest SHA is
`b285f4ddccaa92cb86a8d6ac860de4a1d8a6303061bdc294d7f379254d62d623`.
Eight scalar cases ran at PID37525, about0.00115s, exit0. The sole full
symbolic pass had PID37554 and tool session31877, start chunkec9122 and final
chunk30e3c0, exit0. Its recorded elapsed time was1.541427917s and largest
sampled RSS87916544 bytes, within the registered1200s/768MiB budget. The
session is closed; no native engine or background computation was started.

Saved-table checking was PID37578, tool chunk08645f, about0.04268s, exit0,
zero additional full-symbolic or numerical kernel passes. The environment
was Python3.15.0 with python-flint0.9.0. Full exact inputs, code/proof hashes,
timing/count receipts and output hashes are retained in this directory.
The primary nonzero table SHA is
`8446b3fe84136c005d66501a64a7978e1c1adc452bc1993acd25545ddc8ff299`.

Only this B023 directory was written. C022/B015/B018 stayed frozen; there was
no MATLAB, subagent, native C014 call, git stage, commit or push. Stage1 is
**complete and stopped**, despite meeting the earlier reduced-degree gate.
