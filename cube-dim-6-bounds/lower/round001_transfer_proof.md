# Round 001: six-dimensional rectangular-lattice certificate and transfer

## Finite periodic certificate

Let `u,v>0` satisfy `u^6=2` and `v^6=3`, and put

```
a = (1+u^4/2, (2+u^5 v)/4, (1+u)/2, (1+u)/2, (1+u)/2, (1+u)/2).
```

The periodic atomic measure is one unit atom at every point of
`Lambda=a_1 Z x ... x a_6 Z`.  For a count vector
`ell in {1,2,3,4}^6`, write `m=prod ell`, `r=m^(1/6)/2`, and choose the
center block whose coordinate starts are `-floor(ell_j/2)` for odd
`ell_j` and `1-ell_j/2` for even `ell_j`.  The verifier constructs the
intersection of all its radius-`r` center slabs with the positive half
fundamental cell.  It retains a box only if it is nonempty, checks both
center endpoints against both source endpoints in each coordinate (which
implies all `2^6` source-block corners against that same `r`), removes
contained boxes, and takes the exact union.

There are 127 nonempty candidate boxes and 96 retained boxes.  The exact
number-field calculation proves that their reflected union has measure
`R`, with the rational enclosure

```
206513594276568 / 10^14 < R < 206513594276569 / 10^14.
```

This is an interval emitted from rational isolating intervals of width
`10^-45` for both roots, not a floating-point display.  The center cube
of each retained box contains the selected `m` unit atoms and has volume
`(2r)^6=m`, so its atomic average is at least one; additional lattice
atoms can only raise it.  On every relatively open
subbox it is therefore in `{M mu > lambda}` for every `lambda<1`.
Reflection copies lie in disjoint half-cells up to null boundaries, hence
the `2^6` factor in the certificate.

The source verifier is `round001_rectangular_certificate.py`; its output
also contains all 96 count vectors and the 36 exact coefficients in the
monomial basis `u^i v^j (0<=i,j<6)`.  It makes no claim that the 96-box
family exhausts the periodic bad set.

## Passage to a finite Lebesgue L1 input

Fix `delta>0`.  In every retained witness box replace it by its compact
subbox at distance at least `delta` from every non-null face; the finite
set of subboxes has a positive uniform geometric margin.  Choose a
nonnegative `C_c^infty` bump of mass one and support in the infinity-ball
of radius smaller than that margin.  Replacing each lattice atom by its
translate preserves unit mass and makes every selected common-radius cube
contain the full mass `m`; its average remains one.  Thus every selected
subbox lies in the strict superlevel set at any fixed `lambda<1`.

For integer `K`, retain receiver fundamental cells labelled
`[-K,K]^6` and source atoms in `[-K-2,K+2]^6`.  A selected block has each
relative lattice coordinate between `-2` and `2`, including reflected
half-cell copies, so every selected receiver cube uses only retained
atoms.  The source mass is `(2K+5)^6`; the selected receiver volume is
`(2K+1)^6 (R-o_delta(1))`.  Consequently

```
lambda |{M f_{K,delta}>lambda}| / ||f_{K,delta}||_1
 >= lambda ((2K+1)/(2K+5))^6 (R-o_delta(1)).
```

First let `K` tend to infinity, then let `delta` tend to zero and
`lambda` increase to one.  This proves `C_6>=R` for the stated Lebesgue,
centered, axis-parallel, common-radius operator and ordinary finite
smooth `L1` inputs.  No assertion about independent rectangle scales or
non-Lebesgue backgrounds is used.

## Comparison and scope

The exact four-dimensional inherited certificate is the `dimension:4`
entry of `baseline/original_certificates.json`.  Comparison is with its
complete rational upper endpoint, not a displayed decimal.  The lower
rational endpoint above is strictly larger, so this six-dimensional
instance strictly improves the inherited R4 seed.  A separate exact
`Fraction` comparison also returned true for this lower endpoint against
the complete current five-dimensional internal certificate upper endpoint
in `Math-Research-cube-dim-5/cube-dim-5-bounds/results/bounds_ledger.json`
(read-only source SHA-256
`6314b41ac929b0fde404ea3687ea408fc0c744b06f6f9aa6ecf29c03d78b3b60`).
This is a new evaluated six-dimensional
instance of a known construction family, not a new construction mechanism
and not a claim about the public optimum until the literature ledger is
resolved.

## Provenance and actual run

- Reference `baseline/algebraic_bounds.py` SHA-256:
  `8b3cb2dd042332e52ef9ed3360c9cc8fce210ea37f2112584eb9f5a8d0c31030`.
- Reference `baseline/original_certificates.json` SHA-256:
  `030ce0d945b1864e70a6dee24259f0145ed6a71650925524987cf299ef4d6e69`.
- Reference transfer note SHA-256:
  `7ec33f05308b887cfb9923590d6ea727c6cf0ac66e4897356da551d2f04dbfda`.
- This verifier SHA-256 after the recorded run:
  `9072275d47bd66c3ba9ec7de0c435efdba0785c39a02b122ea1ebb960403ec82`.
- Actual command: `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python
  cube-dim-6-bounds/lower/round001_rectangular_certificate.py`.
  It completed successfully in about 10.5 seconds using Python 3.15.0.

The prior packet's code and logs were read-only inputs.  The timing and
certificate output above were produced in this six-dimensional task.
