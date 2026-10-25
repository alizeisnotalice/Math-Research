# B023 stage1: exact partial fractions only

Released by the main agent for one exact symbolic compression at the NEW point
`x0=(50,50,50,45,27,0)/50`, original body32 frame. This is not B015's point.
The sole frozen input is C022's40320 global-negation-merged unordered tuples.
No C022 group enumeration will be repeated, and no working/unfrozen input is
read. Input archive:

`/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round022/reduced_point_preflight/merged_tuples.json.gz`

Compressed SHA256:
`562e3703286fc9a9712efdad30bbebb71433c4f3f15f867c985ee41b80c6f719`.
Decompressed JSON SHA256:
`7deb058449e11fa74d7e1a4658311cb738a38fbbb57ab03118b494138cd39390`.
C022 report SHA256:
`fad56fada09c08e96311fa7bf25f8d5326ca67947415688febb964d47bac72c5`.

Each row's weight is the complete summed alpha, ALREADY containing the
stabilizer12 and both global signs. Do not multiply either again. The sole
outside normalization is `(2*pi)^-6`. Logical terms total967680; represented
point-orbit/cone terms total80640. Zero and repeated nodes are retained.

Finite budget: one traversal of at most40320 tuples for full symbolic partial
aggregation;967680 represented logical terms;1200 seconds from registration;
at most8 scalar rational partial-identity cases. A1 symbolic G coefficient
will also be summed from this same completed table, so the work is explicitly
**one symbolic full kernel computation**, zero numerical kernel calls and
zero evaluated positive times. No C014, G/H numerical evaluator, new point,
cutoff, frame, MC, grid, Sturm/root isolation or global common denominator is
authorized. Raw degree889 is not a reason to stop before cancellation.

For a distinct node b of multiplicity m, expand
`prod_(c!=b)(b-c+u)^(-m_c)` through degree m-1. Its coefficient a_k contributes
`i^(m-k-6)*a_k/(T-i*b)^(m-k)`. Multiply only the stored row weight and sum in
fmpq by `(b,order)`. Exact zero-b coefficients at all even orders must cancel;
odd orders have identically zero real part and must be reported, not skipped.

For b>0, explicitly fold `C_(b,j)=A_(b,j)+(-1)^j*A_(-b,j)`. Generate the
degree-at-most3 local real numerator polynomial over `(z+b^2)^j`, z=T^2.
Its value at z=-b^2 is `(-1)^(j+1)*2^(j-1)*b^j`, nonzero for b>0. Thus each
pole's highest nonzero folded order gives its TRUE denominator exponent;
different b^2 factors are coprime. Local polynomial shifts may also produce
ordinary real principal-fraction coefficients, but no global numerator.

From the same folded table compute `Gamma_G=sum_(b>0)C_(b,1)` and the exact
coefficient `G*pi^5=Gamma_G/128`. A proof must justify all nonzero-pole orders
>=2 having real integral zero, zero-pole cancellation before integration,
and integration in dT rather than dz. This is new symbolic kernel evidence.

Eight predeclared scalar node multisets, used solely as FORMAL rational
identities, never as times at x0:

1. `(1,2,3,4,5,6)`;
2. `(-3,-2,-1,1,2,3)`;
3. `(0,1,2,3,4,5)`;
4. `(0,0,-2,-1,1,3)`;
5. `(1,1,-1,-1,2,3)`;
6. `(-2,-2,-2,1,3,4)`;
7. `(1,1,1,1,-1,2)`;
8. `(0,0,0,0,0,0)`.

Stop with a checkpoint on source/hash/format/count mismatch, unsupported
arithmetic, nonzero final even zero-pole coefficient, timeout, or resident
memory above768MiB. Store partial progress on failure, never a partial table
as a complete result. No automatic retry or increased budget. Source is read
only; all writes stay in this directory. Other agents are active; no edits to
022/015/018 or their outputs, no MATLAB/subagents/stage/commit/push. Completion
of this first stage stops even if the reduced degree is small.
