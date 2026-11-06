B=fzero(@(b)log(b)-1-1/b,[3 4]); mu=log(B); sigma=sqrt(mu);
for p=[4 8 16 32 64 128 256]
 L=p*log(10); fun=@(z)log(erfcx(z/sqrt(2))/2)-z.^2/2-log(sigma)-log(z)+L/2;
 z=fzero(fun,[.1 sqrt(L)+1],optimset('TolX',1e-13)); H=log(L)+log(sigma*sqrt(2*pi));
 z3=sqrt(L)-H/sqrt(L); z5=z3+(2*H-1-H^2/2)/L^(3/2);
 fprintf('10^%3d z=%.15g residual=%+.3g 3term=%+.6g 5term=%+.6g\n',p,z,fun(z),z-z3,z-z5);
end
