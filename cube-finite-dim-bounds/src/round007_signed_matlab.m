u7signed_t=sym('u7signed_t','real');
u7signed_b=(4-6*u7signed_t^2+3*u7signed_t^3)/6;
u7signed_phi=u7signed_t^2-sym(11)/18*u7signed_t^3;
u7signed_ib=2*int(u7signed_b,u7signed_t,0,1);
u7signed_ip=2*int(u7signed_phi,u7signed_t,0,1);
u7signed_jb=2*subs(diff(u7signed_b,u7signed_t),u7signed_t,1);
u7signed_jp=2*subs(diff(u7signed_phi,u7signed_t),u7signed_t,1);
for u7signed_n=3:4
 u7signed_integral=u7signed_n*u7signed_jp*u7signed_ib^(u7signed_n-1)+u7signed_n*(u7signed_n-1)*u7signed_ip*u7signed_jb*u7signed_ib^(u7signed_n-2);
 if u7signed_n==3, u7signed_scale=sym(3)/8;u7signed_expected=sym(1079)/1024;else,u7signed_scale=sym(27)/64;u7signed_expected=sym(13135)/12288;end
 u7signed_D=1-u7signed_scale*u7signed_integral/2^u7signed_n;
 assert(isAlways(u7signed_D==u7signed_expected));
 fprintf('signed n=%d integral=%s D=%s PASS\n',u7signed_n,char(u7signed_integral),char(u7signed_D));
end
