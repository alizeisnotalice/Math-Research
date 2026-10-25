% Parent serial cap-radial-integral cross only; exact determinants in030_check.py.
% No local functions; eval(fileread(...)) is valid. No clear/path/search.
aud030_Cnum=[5892,36800/3,38400,187200,5632];
aud030_cmin=[.1,1/3,1,1,1];
aud030_eps=[1e-6,1e-8];
aud030_expected=[.0024161972247023651177,4.1975153421455868684e-7; ...
    .0001985785386497470137,3.4305629472958446223e-8; ...
    4.1335389526102769874e-5,7.0973070190146339733e-9; ...
    .0004449220695322015830,7.6393342647093405320e-8; ...
    9.9278480495151934895e-5,4.6080988648434128865e-8];
aud030_maxerr=0;
for aud030_i=1:5
    if aud030_i<5,aud030_n=3;aud030_m=10;aud030_det=4;
    else,aud030_n=4;aud030_m=8;aud030_det=16;end
    aud030_C=aud030_Cnum(aud030_i)/pi^aud030_n;
    aud030_cm=aud030_cmin(aud030_i);
    aud030_D=1/(aud030_det*(aud030_n-1)*pi^aud030_n);
    for aud030_j=1:2
        aud030_e=aud030_eps(aud030_j);
        if aud030_n==3
            aud030_A=aud030_D+aud030_C/(2*aud030_cm^2)* ...
                log1p(4*aud030_m^2*aud030_n^2*aud030_cm^2);
            aud030_B=aud030_C/aud030_cm^2;aud030_L=-log(aud030_e);
            % r=e*exp(-u), so integrate the radial envelope directly on u>=0.
            aud030_value=2*pi*sqrt(2)*aud030_e^2*integral(@(aud030_u) ...
                exp(-2*aud030_u).*(aud030_A^1.5+aud030_B^1.5*(aud030_L+aud030_u).^2), ...
                0,Inf,'RelTol',1e-12,'AbsTol',1e-12);
        else
            aud030_B=aud030_C*pi*aud030_m*aud030_n/aud030_cm;
            % r=e*u, integrate the two radial powers rather than endpoint formula.
            aud030_value=2^(1/3)*4*pi*integral(@(aud030_u) ...
                aud030_D^(4/3)*aud030_e^3*aud030_u.^2+ ...
                aud030_B^(4/3)*aud030_e^(5/3)*aud030_u.^(2/3), ...
                0,1,'RelTol',1e-12,'AbsTol',1e-20);
        end
        aud030_relerr=abs(aud030_value/aud030_expected(aud030_i,aud030_j)-1);
        assert(aud030_relerr<1e-10);aud030_maxerr=max(aud030_maxerr,aud030_relerr);
        fprintf('030 cap case=%d eps=%.0e value=%.17g relerr=%.3g\n', ...
            aud030_i,aud030_e,aud030_value,aud030_relerr);
    end
end
aud030_tangent=[1;1;1];aud030_P=3*eye(3)-aud030_tangent*aud030_tangent.';
assert(min(eig(aud030_P))>-1e-12);
aud030_tangent=[1;1];aud030_P=2*eye(2)-aud030_tangent*aud030_tangent.';
assert(min(eig(aud030_P))>-1e-12);
fprintf('030 PASS ten radial cap integrals max_relerr=%.3g; projection PSD n3/n4\n',aud030_maxerr);
