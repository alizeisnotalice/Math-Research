# A020: a valid pair-orbit frame and exact updated fan, no point certificate

The proposed body32 plus30 pair-direction Cauchy frame is a genuine positive
Markov model. Its added symbol is linear on the entire sorted chamber, so
the frozen21 canonical cones remain valid. A new exact beta=1 fan is
exported; no cone-linearity failure was found. This is an alternative-model
theory/fan delivery, not a new g, m_T or C6 certificate.

Only beta=1 was fixed numerically. The actual finite run passed exactly
12 registered structural checks, including three finite symbol vectors.
No kernel, MC, optimization, grid, MATLAB, agents, commit or push was run.
Every earlier source directory and the old A015 checker remain unchanged.

## Real frame, symbol and common-scale interface

The62 physical columns are the32 body representatives and
beta(e_i+e_j),beta(e_i-e_j), i<j. Each parameter is an independent Cauchy
variable of the same common scale t. Thus

    psi_beta=psi_body+beta psi_pair,
    psi_pair=2 sum_(i=1)^6(6-i)z_i, z=sort_desc(|xi|),
    A_beta A_beta^T=(32+10beta^2)I.

The pair identity follows from |a+b|+|a-b|=2max(|a|,|b|). Its sorted
linear coefficient is **(10,8,6,4,2,0)**. The spectral increment is linear
in beta and the Gram increment quadratic; their coincidence42 at beta=1
must not obscure the distinction.

Independent convolution supplies positivity, probability mass, symmetry,
self-adjointness, semigroup structure and L1/L-infinity contractions. The
Gram gives coercivity, hence a continuous strictly positive density and
strong L1 continuity. The generator has H1 domain and H^(1/2) form domain.
For the same new generator and every T>0,

    mass(K_T)=T,   L K_T=p_T-delta0,
    (K_T)_r=r^-1 K_(rT),
    L(K_T)_r=(p_(rT)-delta0)/r.

The endpoint in normalized spatial dilation is rT. The separate inward
point identity is K_T(r y)=r^-5 K_(T/r)(y)>=K_T(y),0<r<=1.
`frame_theory.md` verifies these facts and gives the positive Markov-average
weak-L1 argument for general L1 inputs and every common positive cube
radius. If a new-model global m>T/4096 is proved, the genuine conditional
cost is T/(64m)<64. That premise has not been established here.

## Exact fan export and proof binding

On xi1>=...>=xi6>=0, every new sum and difference has fixed sign. No new
hyperplane cuts an old cone interior. Consequently each old integer ray r
and old rate q acquires the exact rate

    q_new=q+10r1+8r2+6r3+4r4+2r5.

All126 ray-rate entries were transformed; rays and determinants were
preserved. This mandatory artifact transformation was not a126-point
experiment. The frozen geometry coverage audit and exact integer-scaling
bridge were hash-verified instead of rerunning old fan enumeration.

The compatible exported files are:

* `pair_orbit_beta1_integer_fan.txt`:21 cones in the original integer
  reader format, with new rates42 through230.
* `pair_orbit_beta1_fan.json`: old/new rates, exact pair increments,
  preserved determinant/rays and new per-cone det/product(q_new) weights.

Their SHA-256 values are respectively

    4202b6c685567e51fbd53275add62b6ed978176c834b9fb3f3aaac1f1e5f8a0e
    913f40e375e95dcd57ae94ea6b30ea8f68c0e53450d49484f3ff5f6d61ff9ed7

For example the first old rate72 at ray(2,1,1,1,1,0) becomes112; the
all-ones ray rate60 becomes90. The whole signed-permutation group size
remains46080. No full spectral integral or kernel value was computed.

The fan's geometry being unchanged does not make the model unchanged:
previously artificial Weyl facets may now be genuine pair-symbol kinks.
Source records, genericity/range checks, divided-difference collision forms,
data hashes and native evaluators must be rebound to the new rates. Old
point certificates and old global/singular lower constants are not reused.

## Actual S2 increment, including all dependencies

The geometric classification is proved in `S2_and_rank.md`:

|Class|Direction pairs|Distinct planes added to old S2|
|---|---:|---:|
|Body-body|496|0|
|Sign-matching body-pair|480|0; aliases240 old planes|
|Orthogonal body-pair|480|480|
|Same-support pair-pair|15|15 coordinate2-planes|
|One-overlap pair-pair|240|80 triangle planes|
|Disjoint-support pair-pair|180|180 planes|

There are1891 direction pairs and **1251 distinct S2 planes**, namely the
old496 plus **755 new planes**. Counts concern planes as objects; their
tubes and their point sets intersect and are not disjoint spatial events.

The one-overlap triangle has three pair lines, e.g.
(e1+e2)-(e1+e3)-(e2-e3)=0. A matching body plane also contains its pair
restriction, giving a dependent body/body/pair triple. The classification
derives exactly320 dependent triples without an exhaustive triple scan.
Large-variable arguments must select an actual independent subset by rank,
not assume every triple is independent.

The proof also establishes the new true singular set qualitatively, without
old tube constants. On a selected pair plane, extend that pair to an
invertible six-column basis. In the positive convolution restrict the
remaining sum to a ball of radiusLt, a fixed positive-probability event
after Cauchy scaling. The solved pair parameters stay bounded, giving two
factors at least c t, while the other four basis parameters are O(t),
giving four factors at least c/t. The exact basis determinant is retained,
so p_t(x)>=c' t^-2 for small t and K_T(x)=G(x)=infinity on S2. This is an
ordinary positive-event argument, not zero-fiber conditioning.

Away from S2, a large-set rank-greedy selection provides three independent
large variables. Three tail factors t^3 cancel a three-dimensional quotient
peak t^-3 after coarea; the finite sum over independent triples uniformly
bounds p_t near zero on compact sets. The large-time t^-6 peak is integrable.
Hence G and K_T are finite and continuous off S2. These qualitative facts
do not supply a competitive global minimum or an explicit new tube radius.

## Time normalization and Taylor/checker limits

The coordinate marginal is exactly Cauchy of scale(32+10beta)t, so use
tau=(32+10beta)T as the comparison time. At beta=1, T=tau/42. Matching the
old base T=.048 marginal time tau=192/125 would give the single exact
normalized cutoff **T0=32/875**. The exact new/old symbol ratios are21/16
on an axis and3/2 on the all-ones ray, so this is not a uniform rescaling
of the old model. T0 is the single first-probe normalization recommended
here, not an optimum or a successful cutoff. No value was evaluated there.

The base T48 and beta1 axis T=2/45 models have already failed their own
point tests (the latter point cost is about69.109 as reported by main).
Neither that failure nor the improved Gram predicts the pair-frame result.
The new source/fan requires its own subsequent point certification by main/
Astra, followed by genuine global coverage if still viable.

The A015 negative-Hessian argument extends mathematically to this frame:

    D_h^2 B>=-2||h||2^2 J2_B/(32+10beta^2),
    Taylor loss=D U2/(32+10beta^2).

Phi is still available by translating at most six original body parameters
and leaving the30 pair parameters fixed. All band data and source identities
must be new-model matched. The frozen A015 code/schema only permits body
and axes; it was not modified and no operational pair-checker support is
claimed. No old tube, g or source-data numerical certificate is inherited.

## Execution and stop status

`build_frame.py` exited0, script time **0.017090125 seconds** (tool chunk
e25074), Python PID **36464**. It ran exactly the
12 frozen structural checks: two symbol/Gram identities, the mandatory
fan transformation, one rational time normalization, five canonical pair
structures, and three finite symbol vectors. No new fan arrangement,
kernel, MC, optimization, grid or MATLAB was run.

The nine retained files are the registration, two proof documents, this
report, exact builder, two new fan files, result record and input manifest.
All source and artifact hashes are in `input_manifest.json`/`results.json`.
The finite theory batch is complete and stopped. A successful point,
competitive lower bound or C6 improvement remains unproved.

Reproduce only these structural operations from this directory:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python build_frame.py
