% Independent floating-point check of the cutoff interval in norm_geometry.tex.
fprintf('NORM-GEOMETRY MATLAB CHECK\n');
for n = 1:12
    root = (n + 4 - sqrt(n^2 + 16)) / 2;
    T = min(1, root);
    tau = 0.8 * T;
    lower = tau / (2*n);
    upper = (n - 2*tau) / (n*(n - tau));
    assert(tau > 0 && tau < T && lower < upper);
    fprintf('n=%d T=%.16g interval=(%.16g, %.16g)\n', n, T, lower, upper);
end
fprintf('PASS\n');
