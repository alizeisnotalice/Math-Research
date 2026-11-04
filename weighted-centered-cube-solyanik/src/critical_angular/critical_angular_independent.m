function critical_angular_independent()
% Independent double/log-domain bisection check.  No direct exp(-H) is used.
ns=[2,3]; gs=[.25,.5,.75,.9]; ks=[8,16,32,64,128,256,512,1024,4096]; maxRes=0; count=0; geomCount=0;
for n=ns
 for gamma=gs
  beta=1-gamma; tau=n*beta; c=(n-1)/n; b=(n-1)*beta/(n-gamma);
  for k=ks
   H=k*log(2); L=localL(H,c); kappa=b*c^(-beta); C=H^(1/(2*n)); delta=1/H;
   Hup=H-n*log(C); Hlo=H-n*log(1-2*delta)+n*log(delta);
   lognorm=-beta*H+gamma*log(H); local=exp(localLogF(H,c,b,beta)-lognorm);
   upper=exp(-tau*log(C)+localLogF(Hup,c,b,beta)-lognorm);
   lower=exp(tau*log(delta)+localLogF(Hlo,c,b,beta)-lognorm);
   assert(Hup>0 && Hlo>0 && all(isfinite([local,upper,lower,kappa])));
   geometric=(C>2 && delta<1/3); geomCount=geomCount+geometric;
   maxRes=max(maxRes,abs(L-log(1+c*L)-H)); count=count+1;
  end
 end
end
assert(maxRes<1e-11); fprintf('critical-angular MATLAB: %d rows, %d geometry-valid, max residual %.3g\n',count,geomCount,maxRes);
end
function L=localL(H,c)
lo=0; hi=H+log(1+c*(H+2))+2;
for j=1:800
 mid=(lo+hi)/2;
 if mid-log(1+c*mid)-H<0, lo=mid; else, hi=mid; end
end
L=(lo+hi)/2;
end
function y=localLogF(H,c,b,beta)
L=localL(H,c); y=-beta*L+log(1+b*L);
end
