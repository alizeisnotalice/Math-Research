% Parent serial4D normalization/geometry/fixed-cost cross, not optimizer/full2440.
% No local functions, clear/path changes, author engine or cache calls.
aud035_R={sym([1 1 2 3;0 1 1 1;0 0 1 1;0 0 0 1]), ...
    sym([2 1 1 1;1 1 1 1;1 0 1 1;0 0 0 1]), ...
    sym([2 1 1 3;1 1 1 1;1 0 1 1;0 0 1 1])};
aud035_rates={sym([8 8 16 24]),sym([16 8 12 12]),sym([16 8 12 24])};
aud035_volume=sym(0);aud035_physical=sym(zeros(1,3));aud035_expected_det=[1,1,2];
for aud035_i=1:3
    aud035_det=abs(det(aud035_R{aud035_i}));
    assert(isAlways(aud035_det==aud035_expected_det(aud035_i)));
    aud035_volume=aud035_volume+384*aud035_det/(24*prod(aud035_rates{aud035_i}));
    aud035_tet=aud035_R{aud035_i}./aud035_R{aud035_i}(1,:);
    aud035_physical(aud035_i)=abs(det(aud035_tet(2:4,2:4)-aud035_tet(2:4,1)))/6;
end
assert(isAlways(aud035_volume==sym(11)/4608));
assert(isAlways(24*aud035_volume/16==sym(11)/3072));
assert(isequal(aud035_physical,[sym(1)/36,sym(1)/12,sym(1)/18]));
aud035_t=sym('aud035_t','positive');
assert(isAlways(int(1/((aud035_t^2+1)*(aud035_t^2+4)),aud035_t,0,Inf)==sym(pi)/12));
assert(isAlways(int(1/(aud035_t^2+1)^2,aud035_t,0,Inf)==sym(pi)/4));
aud035_g=sym(21)/25000;aud035_I=sym(92072351)/500000000000;
aud035_cube=256*aud035_I^3/aud035_g^4;aud035_target=sym(14753)/1000;
aud035_gap=aud035_target^3-aud035_cube;
assert(isAlways(aud035_cube==sym('2275587701862626480057')/sym('708750000000000000')));
assert(isAlways(aud035_gap==sym('212329869322269943')/sym('708750000000000000')));
assert(isAlways(aud035_gap>0));
fprintf('035 PASS4D spectral volume11/4608,peak11/3072,physical tets1/36,1/12,1/18,scalar pi calibrations\n');
fprintf('035 fixed cost cube=%s target=14.753 strict_gap=%s cost=%s\n', ...
    char(aud035_cube),char(aud035_gap),char(vpa(aud035_cube^(sym(1)/3),40)));
fprintf('035 source2440 residuals and full proof are independent Python/audit;MATLAB finite cross only.\n');
