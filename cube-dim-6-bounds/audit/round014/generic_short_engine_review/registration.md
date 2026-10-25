# F014 independent short-potential engine audit

Only this directory is writable. C013/C014 and F013 remain frozen.
Budget: 20 minutes, one saved-data run, the same six registered scalar
configurations (five nonzero, one zero rejection), no full kernel,
reference sweep, binary launch, compilation, MC, MATLAB or agents.

Reuse the F013 exact-node, integer-range, permutation/chain and ABI audit.
Independently check f_T=atan(T/b), its required f' barycentric term,
Hermite jets through order6, copied-node derivatives, density rational
identity, signs/normalization and joint Euler identity. Reconstruct all
saved endpoints and source comparisons with the correct inclusion
directions. Build scalar references using rational polynomial Hermite
interpolation, symbolic rational jets and Arb1536 values, not C's DD
recurrence. Read all C evidence without executing its scripts.
