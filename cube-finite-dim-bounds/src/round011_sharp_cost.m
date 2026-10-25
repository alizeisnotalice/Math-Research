% Finite symbolic/decimal calibration only; coefficient proof is separate.
u11q_ratios=[sym(1)/24,sym(17)/640,sym(482191)/13088544];
assert(isAlways(u11q_ratios(1)>u11q_ratios(2)));
assert(isAlways(u11q_ratios(3)>u11q_ratios(2)));
assert(isAlways(sym(17)/640>sym(21)/25000*(sym(22)/7)^3));
u11q_I=sym(92072351)/500000000000;
u11q_g=sym(17)/(640*sym(pi)^3);
u11q_cost=(sym(256)*u11q_I^3/u11q_g^4)^(sym(1)/3);
disp(vpa(u11q_g,40));disp(vpa(u11q_cost,40));
fprintf('11Q ratio ordering PASS; decimal cost is calibration, not independent full kernel proof.\n');
