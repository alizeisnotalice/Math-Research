% One new exact rational cross-check for the completed global upper.
% No figures, pool, path changes, clear, kernels or old scalar rerun.
cube17_dir = fileparts(mfilename('fullpath'));
cube17_now = jsondecode(fileread(fullfile(cube17_dir,'primary_global_acceptance.json')));
cube17_old = jsondecode(fileread(fullfile(cube17_dir,'..','round016','primary_upper139_acceptance.json')));
cube17_g = str2sym(cube17_now.qualified_g);
cube17_g0 = str2sym(cube17_now.nominal_g);
cube17_gp = str2sym(cube17_old.g);
cube17_m = str2sym(cube17_now.minimum_saved_log_margin);
cube17_i = str2sym(cube17_now.complete_I_upper);
cube17_ip = str2sym(cube17_old.complete_I_upper);
cube17_c4 = str2sym(cube17_now.best_upper.fourth_power);
cube17_old4 = str2sym(cube17_old.best_upper.fourth_power);
cube17_nom4 = str2sym(cube17_now.nominal_g110_cost.fourth_power);
cube17_l = str2sym(cube17_now.best_upper.interval{1});
cube17_u = str2sym(cube17_now.best_upper.interval{2});
cube17_s = str2sym(cube17_now.simple_rational_upper);
cube17_checks = false(1,5);
cube17_checks(1) = isAlways(cube17_m>0) && isAlways(cube17_g==cube17_g0*(1+cube17_m)) && isAlways(cube17_g>cube17_g0) && isAlways(cube17_g0>cube17_gp);
cube17_checks(2) = isAlways(cube17_i==cube17_ip) && isAlways(cube17_c4*cube17_g^5==3125*cube17_i^4);
cube17_checks(3) = isAlways(0<cube17_l) && isAlways(cube17_l^4<=cube17_c4) && isAlways(cube17_c4<=cube17_u^4);
cube17_checks(4) = isAlways(cube17_c4<cube17_nom4) && isAlways(cube17_nom4<cube17_old4) && isAlways(cube17_old4<sym(32)^4);
cube17_checks(5) = isAlways(cube17_u<=cube17_s) && isAlways(cube17_s^4*cube17_g^5>=3125*cube17_i^4) && isAlways(cube17_s<32);
fprintf('CUBE17_CERT=[%d,%d,%d,%d,%d]\n',cube17_checks);
assert(all(cube17_checks),'The new exact global upper cross-check failed');
