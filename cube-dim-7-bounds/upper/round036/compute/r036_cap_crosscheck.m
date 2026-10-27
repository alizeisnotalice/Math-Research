% R036 independent cross-check. Uses no path, cwd, or session resets.
r036_zsym = sym('r036_zsym','positive');
r036_poly = r036_zsym^4/sym(24)-sym(4096)*r036_zsym;
r036_dpoly = diff(r036_poly,r036_zsym);
fprintf('poly128_exact=%s\n',char(subs(r036_poly,r036_zsym,128)));
fprintf('deriv128_exact=%s\n',char(subs(r036_dpoly,r036_zsym,128)));
fprintf('second_derivative=%s\n',char(diff(r036_dpoly,r036_zsym)));
for r036_zvalue = [128 256 512]
    r036_zexact = sym(r036_zvalue);
    % log((L*d/32)/(2*z)) and log((51200 L^(-5/4))/(50z L^(-7/4))).
    r036_log_gap = r036_zexact/2-log(sym(64))-log(r036_zexact)/2;
    r036_log_capgap = log(sym(1024))-log(r036_zexact)+r036_zexact/2;
    fprintf('z=%d log_threshold_gap=%s log_cap_envelope_gap=%s\n',r036_zvalue,char(vpa(r036_log_gap,40)),char(vpa(r036_log_capgap,40)));
    assert(isAlways(r036_log_gap>0));
    assert(isAlways(r036_log_capgap>0));
end
assert(isAlways(subs(r036_poly,r036_zsym,128)>0));
assert(isAlways(subs(r036_dpoly,r036_zsym,128)>0));
fprintf('PASS: R036 alternate exponent checks at registered points; derivatives support monotonicity.\n');
