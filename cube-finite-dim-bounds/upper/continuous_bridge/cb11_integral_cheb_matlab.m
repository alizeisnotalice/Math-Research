% Floating calibration of locked continuous pilot/chord formulas, not rigorous audit.
cb11_integral_cheb_root=fileparts(mfilename('fullpath'));
cb11_integral_cheb_data=jsondecode(fileread(fullfile(cb11_integral_cheb_root,'cb11_integral_cheb_four_result.json')));
cb11_integral_cheb_eps=9/5000;cb11_integral_cheb_h=1/256;cb11_integral_cheb_rho=5/2;
cb11_integral_cheb_delta=cb11_integral_cheb_eps-cb11_integral_cheb_h*(cb11_integral_cheb_rho-1/cb11_integral_cheb_rho)/6;
cb11_integral_cheb_M=15/(1088*pi^3*cb11_integral_cheb_delta^2);
cb11_integral_cheb_error=66*4*cb11_integral_cheb_M*cb11_integral_cheb_rho^(-32)/(cb11_integral_cheb_rho-1);
assert(abs(cb11_integral_cheb_error-7.7062935059281956e-8)<1e-20);
for cb11_integral_cheb_k=1:numel(cb11_integral_cheb_data.patches)
 cb11_integral_cheb_p=cb11_integral_cheb_data.patches(cb11_integral_cheb_k);
 cb11_integral_cheb_L=cb11_integral_cheb_dyadic(cb11_integral_cheb_p.L);
 cb11_integral_cheb_U=cb11_integral_cheb_dyadic(cb11_integral_cheb_p.U);
 cb11_integral_cheb_m=cb11_integral_cheb_dyadic(cb11_integral_cheb_p.mean_H_upper);
 cb11_integral_cheb_s=(cb11_integral_cheb_U+sqrt(cb11_integral_cheb_U*cb11_integral_cheb_L)+cb11_integral_cheb_L)/(sqrt(cb11_integral_cheb_U)+sqrt(cb11_integral_cheb_L));
 cb11_integral_cheb_v=(cb11_integral_cheb_L^(3/2)+cb11_integral_cheb_s*(cb11_integral_cheb_m-cb11_integral_cheb_L))*cb11_integral_cheb_p.weight/128^2;
 cb11_integral_cheb_saved=cb11_integral_cheb_dyadic(cb11_integral_cheb_p.IH_contribution_upper);
 assert(abs(cb11_integral_cheb_v-cb11_integral_cheb_saved)<1e-18);
end
cb11_integral_cheb_weights=0;
for cb11_integral_cheb_i=0:127
 for cb11_integral_cheb_j=0:cb11_integral_cheb_i
  cb11_integral_cheb_weights=cb11_integral_cheb_weights+1+(cb11_integral_cheb_i~=cb11_integral_cheb_j);
 end
end
assert(cb11_integral_cheb_weights==128^2);
cb11_integral_cheb_Etime=2^(-4/3)*(3/4)*((256/3)/pi^3)^(1/3)*((16/3)/pi)^(2/3)*cb11_integral_cheb_eps^(4/3);
cb11_integral_cheb_Ballowed=((29/50000)^(2/3)-cb11_integral_cheb_Etime)^(3/2);
fprintf('cb11_integral_cheb fourpatch chord/area PASS; M=%.17g,error=%.17g,Etime=%.17g,Ballowed=%.17g; full I requires all8256\n',cb11_integral_cheb_M,cb11_integral_cheb_error,cb11_integral_cheb_Etime,cb11_integral_cheb_Ballowed);
function out=cb11_integral_cheb_dyadic(v)
if iscell(v)
 a=v{1};b=v{2};
else
 a=v(1);b=v(2);
end
if isnumeric(a),a=num2str(a);end
if ischar(b)||isstring(b),b=str2double(b);end
out=str2double(a)*2^double(b);
end
