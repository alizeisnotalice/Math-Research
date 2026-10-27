# R024 optimized threshold parameter family

The registered exact candidate gives the strict operator upper bound

    C7 < 46644289559722 / 10^12 = 46.644289559722.

Its 512-bit expression interval is approximately 46.6442895597215742010 and is strictly below the full R023 winner interval, approximately 46.6504462333024540772. The simpler exact choice a=36 with τ=√(A(36)B(36)) independently gives C7<46.644289619643.

The frozen theorem proves that, after eliminating τ, the particular parameter family

    F(a) = [√((1+2/a)^7) + √((1+1/a)^7(1+7 log a))]^2

has one global minimizing a>1. Strict 512-bit signs followed by 80 exact dyadic bisections place that family minimizer in

    36.0044448859592342254287786028
       < a* <
    36.0044448859592342254287802572,

with exact bracket width 2⁻⁷⁹. The precise rational endpoints and candidate midpoint are in result.json. Evaluating F at that midpoint supplies the displayed legal upper bound because the proven maximal-operator inequality holds for every admissible parameter. The candidate's value interval is **not** a lower bound on C7 or an interval enclosure of the exact family-minimum value; no such enclosure is claimed. Uniqueness is only for this proved family of bounds.

The finite run performed exactly 84 registered evaluations: F(36), two endpoint signs, 80 midpoint signs and one candidate F. Every sign interval was strictly separated from zero. The two valid upper expressions were compared against R023's complete hash-verified interval and strictly rounded upward on a 10⁻¹² rational grid. No old a=36,τ=7 row, extra parameter, optimizer, phase or Monte Carlo sample was evaluated.

The run used 83 log calls and 170 square-root calls including the candidate's algebraic cross-check; cached log36 was shared. Wall time was 0.01249 s, CPU time 0.00853 s and peak RSS 31,195,136 bytes, in one process. The checker repeats only the same 84 fixed inputs, using 84 log and 166 square-root calls, and verifies every bisection update, dyadic interval and ceiling, and compares the full prior interval. No MATLAB, commit or push was used. The all-L¹/all-positive-radius operator conclusion depends on the separately frozen theorem.
