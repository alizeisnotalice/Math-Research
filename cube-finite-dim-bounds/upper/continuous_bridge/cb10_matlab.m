% cb10 one fixed mu_c(e1), numerical closed-form cross-check only.
% No nested Green quadrature, search, clear, or path changes.
cb10_root=fileparts(mfilename('fullpath'));
cb10_c=(4*log(2)-1)/(24*pi^3);
cb10_plus=integral(@(cb10_r)max(cb10_c-cb10_Q(cb10_r,1),0),0,3/2,'Waypoints',[1/4096,1/4,1/2,3/4,1,5/4],'AbsTol',1e-12,'RelTol',1e-11);
cb10_minus=integral(@(cb10_r)max(cb10_c-cb10_Q(cb10_r,-1),0),0,3/4,'Waypoints',[1/4096,1/8,1/4,1/2,5/8],'AbsTol',1e-12,'RelTol',1e-11);
cb10_mu=4*(cb10_plus+cb10_minus)/pi;
cb10_lo=329537969/125000000000;cb10_hi=2636761891/1000000000000;
cb10_central_margin=15/(176*pi^3)-cb10_c;
assert(cb10_central_margin>0 && cb10_mu>cb10_lo && cb10_mu<cb10_hi);
fprintf('cb10 PASS fixed e1: mu %.17g inside rigorous rational interval [%.12f,%.12f], central margin %.17g; no inner G quadrature.\n',cb10_mu,cb10_lo,cb10_hi,cb10_central_margin);
function cb10_val=cb10_Q(cb10_r,cb10_sg)
 cb10_val=zeros(size(cb10_r));cb10_i=cb10_r>0;cb10_x=cb10_r(cb10_i);
 cb10_val(cb10_i)=(-cb10_x.^2.*(3+cb10_sg*cb10_x)./(2+cb10_sg*cb10_x).^2.*log(cb10_x/2)+(cb10_sg*cb10_x-1).*log(1+cb10_sg*cb10_x/2)+cb10_sg*cb10_x./(2+cb10_sg*cb10_x))/(8*pi^3);
 % Q_sigma(0)=0 is the true continuous limit of r^2 G near the body line.
end
