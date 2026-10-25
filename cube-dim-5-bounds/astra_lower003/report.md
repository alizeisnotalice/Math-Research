# Astra lower003: genuine halo-hole insertion and free-threshold optimization

This round gives a new, delimited obstruction. At one satellite displacement
strictly inside a genuine baseline halo hole, the first-order ratio remains
negative for every fixed real threshold slope h in λ=1+hε. It does not
improve C, and does not exclude all displacements or finite-mass inputs.

Round001 covered only h=0 and the common codilation corresponding to h=1.
Those tests did not settle free h. This round derives the actual profile

    D(z,h)=B0+1−C+h(C−K)−Σ_j A_j[α_j h−β_j]_+,

with exact source period mass1+ε and ratio λU/(1+ε). The+1 is the isolated
satellite halo gain, valid because the selected position is outside the
baseline finite halo. The threshold and denominator term is(h−1)C.
The common-radius speed is r(n_sat/N−h)/5. Full derivation and smooth-L1
transfer are in `threshold_proof.md`.

The exact primary endpoint masks select one uncovered open cell with a
deterministic symmetry-first rule, then refine it by satellite-count events.
The tested centre is approximately

    (0.5370890252561094,
     0.005487417743144651,0.005487417743144651,
     0.3126971573211914,0.3126971573211914).

Exact coefficients and cell boundaries are in `selected_cell.json`. The
open cell is disjoint from all42 retained primary receiver boxes. Only its
centre's derivative was evaluated, not the whole cell's sign. The first
selection launch stopped on an overly strong mask-symmetry assertion before
selecting a candidate; it was corrected to compare coordinate cuts while
retaining the distinct per-coordinate box masks. No result from that failed
launch was used.

The finite baseline has no zero-width raw boxes and no unequal-radius
same-oriented-plane coincidences. Consequently its free-h derivative is
concave piecewise affine: exposed boundaries contribute an affine term,
while internal contacts open nonnegative hinge gaps. These assertions were
checked on all51 raw boxes, not assumed from the picture.

The18 reflection classes, whose multiplicities sum to32, produce333 nonzero
contact-count pieces and20 exact h events. Their normalized result is:

| Quantity | Decimal display only |
|---|---:|
| Global first-order maximizer |h=0|
| Maximum derivative |−0.1317150987176058|
| Left h slope at0 |+0.0905528353751288|
| Right h slope at0 |−0.263635885566582|
| Adjacent event h |0.089445370485576|
| Derivative at adjacent event |−0.15529610817540165|

Only these two h candidates were evaluated. All h-event ordering and sign
decisions use Q(u,v) cross products and rational isolating intervals. A
particularly short global proof follows because the zero-event hinge weight
already exceeds the initial positive h slope. Therefore the response
increases up to h=0 and decreases thereafter; later hinges only decrease it
further. This certifies sup_(h∈R)D(z,h)<0 for the selected z.

Each reflection's exposed-boundary calculation took about0.93--1.71 seconds;
contact integration took about0.09--0.16 seconds. A separate coordinate-mask
union engine reproduced all333 nonzero contact areas. Literal one-dimensional
interval unions reproduced333 gap losses at the positive h event, and the
normalized derivative numerator matched exactly. These checks share scalar
and baseline/contact geometry and are not described as an independent
geometry audit.

The attempted direct5D mixed-receiver jet audit is **incomplete**. Its first
reflection hit the registered120-second checkpoint inside recursive union,
session5003, exit1. No direct5D result was produced; the remaining17
reflections were not started. `direct_jet_checkpoint.json` records this,
and the source is retained for a later registered chunked continuation.
No resumable union cache was saved. The shorter mask/gap path is separately
registered and does not turn the failed5D attempt into a claimed pass.

The successful verification entry point is:

```sh
python3 -B cube-dim-5-bounds/astra_lower003/verify_results.py
```

To regenerate geometry and the threshold profile, run `select_geometry.py`,
`threshold_kernel.py geometry`, then `threshold_kernel.py reflection i` for
i=0,…,17, followed by `threshold_kernel.py finalize` and `contact_check.py`.
The pending full5D attempt is deliberately not part of the successful replay.
All code uses local Python and exact arithmetic. MATLAB and the main clone
were untouched; no commit or push was made.

The remaining lower-bound question is whether another hole or multi-coordinate
cell has a sufficiently larger mixed-boundary gain, or whether a different
finite family avoids the contact losses. This centre's failure does not
answer that question. There is no new record lock or operator upper bound.
