% Parent serial symbolic conversion/error-constant cross only, not fresh nodes.
% No local functions, clear/path changes, cache or Green calls.
aud033_u=sym('aud033_u','real');aud033_n=20;
for aud033_k=[0,1,7,20]
    aud033_p=expand(chebyshevT(aud033_k,2*aud033_u-1));
    aud033_power=fliplr(coeffs(aud033_p,aud033_u,'All'));
    aud033_reconstructed=sym(0);
    for aud033_j=0:aud033_n
        aud033_b=sym(0);
        for aud033_l=0:min(aud033_j,aud033_k)
            aud033_b=aud033_b+aud033_power(aud033_l+1)* ...
                sym(nchoosek(aud033_j,aud033_l))/nchoosek(aud033_n,aud033_l);
        end
        aud033_reconstructed=aud033_reconstructed+aud033_b*nchoosek(aud033_n,aud033_j)* ...
            aud033_u^aud033_j*(1-aud033_u)^(aud033_n-aud033_j);
    end
    assert(isAlways(expand(aud033_reconstructed-aud033_p)==0));
end
aud033_eps=sym(1)/50;aud033_h=sym(1)/32;aud033_rho=sym(3);
aud033_s=aud033_h*(aud033_rho-1/aud033_rho)/2;
aud033_delta=aud033_eps-aud033_s/3;
assert(isAlways(aud033_s==sym(1)/24));
assert(isAlways(aud033_delta==sym(11)/1800));
aud033_M=15/(1088*sym(pi)^3*aud033_delta^2);
aud033_D=15/(1088*sym(pi)^3*(aud033_eps-aud033_eps/3)^2*aud033_eps);
aud033_nodeerror=2*aud033_D*aud033_h/sym(2)^128;
aud033_E=84*aud033_M/sym(3)^20;
assert(double(aud033_E)<2.868313237155e-7);
assert(double(aud033_nodeerror)<2.296921e-38);
fprintf('033 PASS Bernstein degree0/1/7/20 exact identities;M=%.17g D=%.17g E=%.17g nodeerror=%.17g\n', ...
    double(aud033_M),double(aud033_D),double(aud033_E),double(aud033_nodeerror));
fprintf('033 full59976-node certificate is independent Python;this MATLAB script only finite algebra cross.\n');
