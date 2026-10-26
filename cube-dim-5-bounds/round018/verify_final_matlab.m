% New018 exact arithmetic crosscheck in a private MATLAB batch session.
cube18_dir = fileparts(mfilename('fullpath'));
cube18_read = @(n) jsondecode(fileread(fullfile(cube18_dir,n)));
cube18_q = @(s) str2sym(char(s));
cube18_a = cube18_read('primary_integral156_acceptance.json');
cube18_b = cube18_read('primary_curvature155_acceptance.json');
cube18_c = cube18_read('primary_integral150_acceptance.json');
cube18_g = cube18_q(cube18_a.qualified_g);
cube18_m = cube18_q(cube18_b.minimum_margin110);
cube18_poly = sym(1); cube18_term = sym(1);
for cube18_k = 1:8
    cube18_term = cube18_term*cube18_m/cube18_k;
    cube18_poly = cube18_poly+cube18_term;
end
cube18_exactfloor = sym(11)/100000*cube18_poly;
cube18_I = cube18_q(cube18_a.complete_I_upper);
cube18_H = cube18_q(cube18_a.new_UH);
cube18_J = cube18_q(cube18_a.same_complete_J139_upper);
cube18_fourth = cube18_q(cube18_a.best_upper.fourth_power);
cube18_lo = cube18_q(cube18_a.best_upper.interval{1});
cube18_hi = cube18_q(cube18_a.best_upper.interval{2});
cube18_prev4 = cube18_q(cube18_b.best_upper.fourth_power);
cube18_prevlo = cube18_q(cube18_b.best_upper.interval{1});
cube18_simple = cube18_q(cube18_a.simple_rational_upper);
cube18_ok = [isAlways(cube18_m>0 & cube18_g>cube18_q(cube18_b.linear_floor)), ...
    isAlways(cube18_g<=cube18_exactfloor & cube18_exactfloor-cube18_g<sym(2)^(-192)), ...
    isAlways(cube18_g*sym(2)^192==floor(cube18_g*sym(2)^192)), ...
    isAlways(cube18_I==cube18_H+cube18_J & cube18_I<cube18_q(cube18_c.complete_I_upper)), ...
    isAlways(cube18_fourth*cube18_g^5==3125*cube18_I^4), ...
    isAlways(cube18_lo^4<=cube18_fourth & cube18_fourth<=cube18_hi^4), ...
    isAlways(cube18_hi<=cube18_simple & cube18_simple<cube18_prevlo & cube18_prev4>cube18_fourth), ...
    isAlways(cube18_simple<32 & cube18_q(cube18_a.exact_I_gain_vs155)>0)];
assert(all(cube18_ok));
fprintf('CUBE18_CERT=%s\n',jsonencode(double(cube18_ok)));
