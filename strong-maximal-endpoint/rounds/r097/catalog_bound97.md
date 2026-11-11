# r097 finite catalog density bound

This note gives a strict rational lower bound for every source indicator on each of the twelve finite periodic grids in the r097 search. It covers only these grid-aligned rectangles and receiver-cell catalogs. It does not imply a bound for the continuous strong maximal operator, a general dimension, or an asymptotic regime.

The values below were recomputed by [`catalog_bound97.py`](../../code/catalog_bound97.py) using Python `Fraction`. The listed candidate upper bounds come from the exact `source_density_fraction` fields in the existing r097 per-case certificates. Those are realized periodic indicators; the HiGHS dual bounds are not used here.

| dimension | grid | p | theta | minimizing h | s_h | minimizing C | catalog lower bound | certified candidate density upper bound |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 2 | dyadic6 | 19/20 | 9/10 | 23 | 5/16 | 1 | 31/160 | 101712187/134217728 |
| 2 | dyadic6 | 39/40 | 9/10 | 21 | 1/4 | 9/10 | 221/400 | 888985765/1073741824 |
| 2 | dyadic6 | 99/100 | 99/100 | 23 | 5/16 | 99/100 | 16621/20000 | 2058064999/2147483648 |
| 2 | uniform 8 | 19/20 | 9/10 | 32 | 1/2 | 1 | 7/40 | 823763289/1073741824 |
| 2 | uniform 8 | 39/40 | 99/100 | 32 | 1/2 | 99/100 | 1171/2000 | 497446033/536870912 |
| 2 | uniform 8 | 99/100 | 3/4 | 24 | 3/8 | 3/4 | 21/32 | 782659393/1073741824 |
| 3 | uniform 3 | 19/20 | 99/100 | 14 | 14/27 | 99/100 | 17119/27000 | 1363153/1572864 |
| 3 | uniform 3 | 39/40 | 3/4 | 11 | 11/27 | 3/4 | 2789/4320 | 166094843/226492416 |
| 3 | uniform 3 | 99/100 | 9/10 | 13 | 13/27 | 9/10 | 7567/9000 | 51170635/56623104 |
| 3 | uniform 4 | 19/20 | 3/4 | 32 | 1/2 | 1 | 7/40 | 21181315/33554432 |
| 3 | uniform 4 | 39/40 | 9/10 | 29 | 29/64 | 9/10 | 1807/3200 | 884162851/1073741824 |
| 3 | uniform 4 | 99/100 | 99/100 | 32 | 1/2 | 99/100 | 4141/5000 | 32595081/33554432 |

In every row the certified candidate upper density is at least the catalog lower bound, as required. The interval between them is not claimed sharp.

## Proof for one finite grid catalog

Let the torus be partitioned into finitely many grid cells with positive rational volumes (v_1,ldots,v_M). The catalog consists of boxes that are unions of whole grid cells. Let (f) be any measurable indicator (the proof also permits any measurable (0\leq f\leq1)) and write (d=\int f).

Suppose a collection of catalog boxes, each with average of (f) strictly greater than (p), covers receiver cells of total volume at least (\theta). Discard boxes until the remaining collection is inclusion-minimal for covering those receiver cells. Let its size be (h), and let (U) be the union of these (h) boxes, with (C=|U|). Then (C\geq\theta).

Each box has a private receiver cell: otherwise that box could be removed while preserving the receiver-cell cover. These (h) private cells are distinct. If their volumes sum to (S), then (S\geq s_h), where (s_h) is the sum of the (h) smallest grid-cell volumes. The multiplicity with which the (h) boxes cover (U) is one on each private cell and at most (h) elsewhere. Therefore

\[
 \sum_{j=1}^h |R_j|
 =\int_U \sum_{j=1}^h 1_{R_j}
 \leq S+h(C-S)
 \leq hC-(h-1)s_h.
\]

Since each box has average greater than (p), its deficit satisfies \(\int_{R_j}(1-f)<(1-p)|R_j|\). The union deficit is at most the sum of the box deficits. Thus

\[
 d\geq\int_U f
 > C-(1-p)\sum_j |R_j|
 \geq [1-(1-p)h]C+(1-p)(h-1)s_h.
\]

Also (C\leq\sum_j|R_j|\leq h\max_j|R_j|), so at least one selected box has volume at least (C/h). Its average condition gives (d>pC/h). Using weak inequalities for the final reported lower bound, for each possible (h) we have

\[
 d\geq \max\left\{[1-(1-p)h]C+(1-p)(h-1)s_h,\;pC/h\right\},
 \qquad C\in[\max(\theta,s_h),1].
\]

For fixed (h), the maximum of these two affine functions is minimized at an interval endpoint or at their crossing. The script checks those rational candidates exactly, then takes the minimum over (h=1,\ldots,M). The actual minimal cover has one such (h), so the minimum over all (h) is a valid lower bound. The strict box-average hypotheses can make the underlying inequality strict; the table reports the safe non-strict rational lower bound.

## Scope and limitations

- Cell volumes, all minimizations, and the reported lower fractions are exact rationals; no floating-point solver output enters the lower-bound computation.
- The result applies to arbitrary indicator sources for the stated finite grid-aligned catalog, not only to the cellwise-constant sources used to produce the candidate upper bounds.
- The candidate upper densities are finite witnesses. The gap between a witness and the catalog lower bound is not a solver gap or a sharpness claim.
- This is a finite periodic catalog statement only. It gives no continuous-box transfer and no global (H_N\) estimate.
