# Independent audit batch 002 — 2026-10-25

Inputs frozen for this batch: exact spacing formula and coefficients in baseline/original_certificates.json; candidate claim concerning the fixed one-orbit rectangular lattice at threshold 1; dimension 3 and 4 only. Read-only comparison source: main clone src/extension_certificate.py. No calls/imports of its verifier or inherited arithmetic.

Mechanisms: (1) analytic derivation of nearest-block completeness in the positive half-cell and strict-level boundary treatment; (2) independently implemented rational root bisection to classify all 1..8 count vectors and tail inequality; (3) independent exact volume by signed intersections/inclusion–exclusion using SymPy quotient-polynomial reduction; (4) MATLAB floating-point cross-check of tail and finite classification only.

Finite budget: one initial Python run and at most two corrections for implementation faults, at most 1,000,000 finite tuple/subset operations per dimension; one MATLAB evaluation. Terminate this batch after exact classification and coefficient comparison, or record a mismatch/unresolved interval instead of claiming success. No optimization, no new parameter search and no high-dimensional expansion. User has removed previous machine-compute hard caps; this explicit batch still bounds the experiment.

Independence disclosure: shared mathematical input/formula, Python standard-library Fraction and inherited certificate data; no inherited checker import. SymPy is a separate polynomial engine. Numerical MATLAB evidence does not constitute proof.
