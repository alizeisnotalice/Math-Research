% Deterministic MATLAB cross-check of the finite layer-cake shell formula.
fprintf('MATLAB %s layer-cake corner expansion\n', version);
cases = [2, 0.5; 2, 1.25; 3, 0.75; 4, 0.5];
for row = 1:size(cases, 1)
    n = cases(row, 1);
    tau = cases(row, 2);
    kappa = 2^n*n/tau;
    expected = -tau*(n-tau)/2;
    for power = [6, 8, 10]
        delta = 2^(-power);
        ratio = cornerMass(n, tau, delta)/kappa;
        quotient = (ratio - 1)/delta^2;
        fprintf('n=%d tau=%.8g delta=2^-%d quotient=%.14g expected=%.14g\n', ...
            n, tau, power, quotient, expected);
        assert(abs(quotient-expected) < 0.1);
    end
end

fprintf('Cutoff compatibility for tau=1/2:\n');
for n = 2:8
    bcrit = (n + 4 - sqrt(n^2 + 16))/(2*n);
    b = 0.5/n;
    lo = b/2;
    hi = (1 - 2*b)/(n - 0.5);
    assert(b < bcrit && lo < hi);
    fprintf('n=%d bcrit=%.14g a in (%.14g, %.14g)\n', n, bcrit, lo, hi);
end

function value = cornerMass(n, tau, delta)
lo = 1-delta;
hi = 1+delta;
first = 2^n*lo^tau/tau;
middle = integral(@(s) (s+1-delta).^n .* s.^(tau-n-1), lo, hi, ...
    'AbsTol', 1e-13, 'RelTol', 1e-13);
tail = 2^n*hi^(tau-n)/(n-tau);
value = (n-tau)*(first + middle + tail);
end
