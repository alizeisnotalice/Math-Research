% Exact content of the MATLAB Symbolic Math Toolbox verification executed
% on 2026-10-27, reformatted into lines. No path reset or clear operation.
syms asy_t asy_r positive;
asy_ok = true;
for asy_n = 1:6
    asy_ell = 1 + 1/(2*asy_t) - 1/(2*asy_r);
    asy_h = (1 + 1/asy_t)/2;
    asy_I = int(expand(asy_ell^(asy_n-1)/asy_r), asy_r, 1, asy_t);
    asy_B = 1 - asy_h^asy_n + asy_n*log(asy_t) - asy_n*asy_I/2;
    asy_D = asy_n/(2*asy_t) + asy_n*asy_h^(asy_n-1)/(2*asy_t^2);
    if asy_n >= 2
        asy_D = asy_D + asy_n*(asy_n-1)/(4*asy_t^2) * ...
            int(expand(asy_ell^(asy_n-2)/asy_r), asy_r, 1, asy_t);
    end
    asy_pass = isAlways(simplify(diff(asy_B, asy_t) - asy_D) == 0);
    fprintf('n=%d exact derivative check=%d\n', asy_n, asy_pass);
    asy_ok = asy_ok && asy_pass;
end
fprintf('ALL_EXACT_CHECKS_PASSED=%d\n', asy_ok);
