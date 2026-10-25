# A008: finite chamber-convexity falsification batch

Test the new hypothesis that body32 G is convex on the ordered positive-face
chamber Delta={x1=1> x2> x3> x4> x5> x6>0}. Do not assume it. This does not
rerun time-radius tests or edit any round007 material.

Budget: at most two strict interior centers, five direction definitions,
15 new generic rational G points, and 20 minutes. Start with just one center
and the last free coordinate direction, h=1/100. The three-point segment
stays strictly in the chamber. Existing C++ direct-kernel code is a numerical
diagnostic only. If its Jensen gap is negative, immediately certify those
same three points by independent exact fmpq arithmetic, then stop. Do not
continue to other directions after a certified negative gap. If no negative
gap is found, any nonnegative finite results remain diagnostics only.

The center, with denominator 1000000007, is
(1000000007,853271829,691287433,517389127,307281763,117361927).
Multiply numerator/denominator by 100 so that +/-1/100 is represented exactly.
All kernel inputs and original source hashes are locked. Only this owned
directory in the independent theory clone is writable. No MATLAB, agents,
commits, pushes, or changes to other workers' files.

Adaptive continuation if the first gap is nonnegative: use directions
e4-e5, e4, e2-e3, e3-e4, with x1 fixed (indices refer to the six-vector).
These five direction definitions in total require at most 11 points at the
first center. Only if all five gaps are nonnegative, use a second strict
center near (1,.9514,.7913,.5914,.3913,.1914) and reuse e4-e5, requiring
three additional points. Thus the actual batch cap is 14 points. Stop on
the first negative diagnostic and certify only that triple.
