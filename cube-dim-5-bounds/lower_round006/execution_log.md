# Execution log

The first attempted `theta=1/8` launch stopped before constructing any box:
Python evaluates `(-1)**j` as a floating value when `j` is negative, so the
mass lost its `Fraction` type.  This was an implementation failure, not a
numerical or mathematical candidate failure.  The script now obtains the sign
from integer parity as an exact `Fraction`, before the first successful run.
