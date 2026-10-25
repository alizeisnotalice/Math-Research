# lower005 — block-span pruning audit

This is a bounded 20-minute analytic/computational audit of length vectors
only.  It evaluates no box union and searches no period-two mechanism.
For `ell_j in {1,...,m}`, `m=5,6`, it uses rational root intervals to test
the necessary common-radius span inequalities and the necessary product
(geometric-mean) inequality obtained by multiplying them.  The acceptance
object is a certified classification into excluded by a coordinate span,
additionally excluded by the product condition, or not excluded by these
necessary tests.  The last class is not claimed nonempty, maximal, contained,
or useful for a lower bound.
