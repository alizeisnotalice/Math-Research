# Period-two half-cell compression: transfer-interface audit

This note audits the stated interface only.  It does not certify the
floating period-two screen or claim a new lower bound.

## Periodic normalization and reflection

Let the atoms be at `(k_1 a_1,...,k_6 a_6)` with mass
`w_{k_1}=1+delta(-1)^{k_1}`, where `|delta|<1`.  A coarse fundamental cell is

```
C=[0,2a1] x product_{j>1}[0,aj],
```

and contains the two first-axis atom indices 0 and 1, of total mass
`(1+delta)+(1-delta)=2`.  The proposed half cell has volume `|C|/64`:

`H=[0,a1] x product_{j>1}[0,aj/2]`.

Reflection in the first coordinate, `x1 -> 2a1-x1`, maps a source index
`k` to `2-k`, preserving `(-1)^k`; the five remaining reflections map
`k_j` to `1-k_j` and do not alter the mass.  Thus the reflected copies of a
valid witness subfamily have equal witness radii and occupy the 64 half-cell
copies of `C` up to null faces.  If `U_H` is the witness union in `H`, its
periodic receiver volume is `64|U_H|`, so at level one its atomic mass ratio
is

```
(64 |U_H|) / 2 = 32 |U_H|.
```

For a first-axis consecutive block starting at `s` with length `ell1`, its
mass is exactly

```
prod_{j>1} ell_j * (ell1 + delta sum_{q=s}^{s+ell1-1} (-1)^q).
```

The alternating sum vanishes for even `ell1` and is `(-1)^s` for odd
`ell1`, agreeing with the stated radius formula.  Under the first reflection
the start becomes `s'=3-s-ell1`; for odd `ell1`, `s'` has the same parity as
`s`, and for even `ell1` mass is start-independent.  Hence every reflected
block has exactly the same mass and common radius.

## Collar audit

The independent `start_range_check.py` re-enumerated all `ell in {1,...,4}^6`
and direct first-axis intersection conditions for `delta=0,1/20,1/10,1/5`.
It found starts exactly `{-1,0,1}`, original first-axis source indices
`[-1,2]`, and, after first-axis reflection, `[-1,3]`.  The latter indices
belong only to coarse cells at offsets `-1,0,1` from the receiver coarse
cell.  In the other five coordinates the old centered starts give source
indices `[-1,2]`; reflection `k_j ->1-k_j` preserves that range, again only
coarse offsets `-1,0,1`.

There is also a uniform elementary bound behind the finite check for the
registered `0<=delta<=1/5`: every first-axis witness radius is below `2`,
while `a1=1+2^(-1/3)>1`.  The literal intersection inequalities
`(s+ell1-1)a1-r<a1` and `s a1+r>0`, with `ell1<=4`, force the displayed
start/atom ranges (the length-one case has the sharper `r<a1`).

Therefore source coarse cells `[-K-1,K+1]^6` already suffice for receivers
`[-K,K]^6`; the proposed conservative collar `[-K-2,K+2]^6` is sufficient.
It contains `2(2K+5)^6` total source mass.  With the 64 reflected receiver
copies, the finite selected witness ratio is consequently bounded below by

```
lambda * 32 |U_H| * ((2K+1)/(2K+5))^6.
```

As usual, choose compact inner witness boxes and smooth each atom by a
unit-mass bump whose infinity support is smaller than the finite positive
geometric margin.  Then the same common-radius cubes contain all selected
bump masses, giving a strict superlevel statement for each fixed
`lambda<1`.  This produces ordinary finite Lebesgue `L1` inputs if a future
period-two union is certified; it does not turn the current floating screen
into such a certificate.
