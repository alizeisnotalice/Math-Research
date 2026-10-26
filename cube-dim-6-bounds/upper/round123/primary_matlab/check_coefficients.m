cube6r123_check = struct('scope','symbolic normalization only; no expectation or C6 certification','matlab_version',version,'symbolic_available',~isempty(ver('symbolic')));
if cube6r123_check.symbolic_available
    cube6r123_pi = sym('pi');
    cube6r123_c = sym(3)/(16*cube6r123_pi^(sym(5)/2));
    cube6r123_k = simplify(cube6r123_c^6*(cube6r123_pi^3/384)^5);
    assert(isAlways(cube6r123_k == sym(3)/sym(2)^59));
    cube6r123_t = (sym(8)/9)^6*sym(24576)^5/sym(65000)^6;
    assert(isAlways(cube6r123_t == sym('36893488147419103232')/sym('629200415790081024169921875')));
    assert(isAlways(cube6r123_t == sym(2)^83/(sym(3)^7*sym(65000)^6)));
    cube6r123_check.coefficient_exact = char(cube6r123_k);
    cube6r123_check.threshold_sixth_power_exact = char(cube6r123_t);
    cube6r123_check.threshold_display_only = char(vpa(cube6r123_t^(sym(1)/6),30));
    cube6r123_check.status = 'PASS_3_SYMBOLIC_NORMALIZATION_IDENTITIES';
else
    cube6r123_check.status = 'SYMBOLIC_TOOLBOX_UNAVAILABLE';
end
cube6r123_check.utc = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss Z'));
assert(~isfile('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round123/primary_matlab/results.json'));
cube6r123_fid = fopen('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round123/primary_matlab/results.json','w');
assert(cube6r123_fid>=0);
fprintf(cube6r123_fid,'%s\n',jsonencode(cube6r123_check));
fclose(cube6r123_fid);
disp(jsonencode(cube6r123_check));
