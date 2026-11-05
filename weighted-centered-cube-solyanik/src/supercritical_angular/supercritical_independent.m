function supercritical_independent()
% Independent log-domain double bisection for the round31 stated formula.
ns=[2,3];gs=[.25,.5,.75,.9];cs=[.01,.05,.1,.25,.5,.75,.95];baseK=[16,64,256,1024];extraK=[4096,16384];count=0;maxResid=0;maxRel=0;interior=0;conservative=0;
for n=ns
 for gamma=gs
  beta=1-gamma;
  for c0=cs
   d=n*gamma*c0;tau=n*beta+d;e=n-tau;P=d/gamma;tailPower=d/e;p=n*gamma/e;coeff=(n-1)/P;mass1=tau*gamma/e;mass2=tau*(n-1)/(P*(n-gamma));A=tau/(n-gamma)*coeff^gamma;
   assert(n*beta<tau&&tau<n&&all([e,P,tailPower,p,coeff,mass1,mass2,A]>0));
   if d/tau<1/3,delta=d/tau;interior=interior+numel(baseK)+2*(c0==.01);else,delta=.25;conservative=conservative+numel(baseK)+2*(c0==.01);end
   assert(delta>0&&delta<1/3);
   for kk=1:(numel(baseK)+2*(c0==.01))
    if kk<=numel(baseK),k=baseK(kk);else,k=extraK(kk-numel(baseK));end
    H=k*log(2);L=solveL(H,p,coeff,tailPower);logV=logVolume(L,p,coeff,tailPower);logM=logMass(L,beta,mass1,mass2,tailPower);ratio=exp(logM+beta*H);assert(isfinite(ratio)&&abs(logV+H)<1e-10);maxResid=max(maxResid,abs(logV+H));maxRel=max(maxRel,abs((ratio-A)/A));count=count+1;
   end
  end
 end
end
B10=(2*sqrt(5)/(3*11^(1/10)))^10;assert(abs(B10-3200000/649539)<1e-13&&B10>(7/6)^10);fprintf('round31 MATLAB: %d rows, interior %d, conservative %d, max log residual %.3g, max relative error %.3g\n',count,interior,conservative,maxResid,maxRel);
end
function L=solveL(H,p,coeff,tailPower)
lo=0;hi=H+log(coeff)+10;assert(logVolume(lo,p,coeff,tailPower)>-H&&logVolume(hi,p,coeff,tailPower)<-H);
for j=1:800
 mid=(lo+hi)/2;if logVolume(mid,p,coeff,tailPower)>-H,lo=mid;else,hi=mid;end
end
L=(lo+hi)/2;
end
function y=logVolume(L,p,coeff,tailPower),y=logadd(-p*L,log(coeff)-L+log1mexp(tailPower*L));end
function y=logMass(L,beta,mass1,mass2,tailPower),y=logadd(-mass1*L,log(mass2)-beta*L+log1mexp(tailPower*L));end
function y=log1mexp(x),if x==0,y=-Inf;else,y=log(-expm1(-x));end,end
function z=logadd(x,y),hi=max(x,y);lo=min(x,y);if isinf(lo)&&lo<0,z=hi;else,z=hi+log1p(exp(lo-hi));end,end
