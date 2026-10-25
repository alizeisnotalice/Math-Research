% Analytic budget calibration only; does not certify I or rerun boxes.
cb11_integral_root=fileparts(mfilename('fullpath'));
cb11_integral_data=jsondecode(fileread(fullfile(cb11_integral_root,'cb11_integral_results.json')));
cb11_integral_J=integral(@cb11_integral_radial,0,1,'AbsTol',1e-11,'RelTol',1e-11);
cb11_integral_A=55200/pi^3;
cb11_integral_K=1600*pi*cb11_integral_A^(3/2)*cb11_integral_J;
cb11_integral_cond=((9/16000)^(2/3)+(cb11_integral_K*2^(-52))^(2/3))^(3/2);
assert(abs(cb11_integral_cond-.0005659300191215578)<1e-12);
assert(cb11_integral_cond<57/100000);
fprintf('cb11_integral budget calibration PASS; J=%.17g, conditionalI=%.17g; actual I NOT certified\n',cb11_integral_J,cb11_integral_cond);
function y=cb11_integral_radial(w)
y=zeros(size(w));m=w>0&w<1;u=w(m)./(1-w(m));y(m)=u.*log1p(1./u.^2).^(3/2)./(1-w(m)).^2;y(w==1)=1;
end
