function flat_concentration_independent()
% Double precision check: unresolved endpoint cases are explicitly skipped.
betas=[.1,.25,.5,.75,.9];ks=[4,8,16,32,64,128,256];solved=0;skipRound=0;skipNear=0;skipSmall=0;maxHalf=0;
for beta=betas
 for kind=1:2
  for k=ks
   e=2^(-k);if kind==1,u=e;else,u=1-e;end
   if kind==2&&u==1,skipRound=skipRound+1;continue,end
   [lo,hi,~,~,ok]=solveB(u,beta);if ~ok,if kind==2,skipNear=skipNear+1;else,skipSmall=skipSmall+1;end;continue,end
   b=(lo+hi)/2;a=u*b;omega=M(a,u,beta);
   if beta==.5,r=sqrt(2*u);ac=u/(2*(1+r+u));oc=(u+r)/(1+r);maxHalf=max(maxHalf,max(abs(a-ac),abs(omega-oc)));end
   solved=solved+1;
  end
 end
end
fprintf('round32 MATLAB: solved %d, skipped-near-one-rounding %d, indeterminate-near-one %d, indeterminate-small %d, max half mismatch %.3g\n',solved,skipRound,skipNear,skipSmall,maxHalf);
% Independent finite-k grid lower diagnostic (not an exact halo computation).
gridCases=0;maxLimitDefect=0;
for nn=[2,3]
 for gamma=[.05,.2,.5]
  beta=1-gamma;ell=gamma/beta;eta=.01;h=(1-eta)*(2/ell)^(nn-1);
  Jlimit=(1-eta)^beta*(2^beta*ell^(1-beta)/(1+ell))^(nn-1);
  Jlo=((10*ell)/(10*(1+ell)+1))^(nn-1)*h^beta;
  Jhi=((10000*ell)/(10000*(1+ell)+1))^(nn-1)*h^beta;
  assert(abs(Jhi/Jlimit-1)<abs(Jlo/Jlimit-1));
  maxLimitDefect=max(maxLimitDefect,abs(Jhi/Jlimit-1));gridCases=gridCases+1;
 end
end
fprintf('round32 MATLAB grid diagnostic: %d cases, max k=10000 relative defect %.3g\n',gridCases,maxLimitDefect);
end
function [lo,hi,Klo,Khi,ok]=solveB(u,beta)
lo=eps;hi=.49;Klo=K(u*lo,u,beta);Khi=K(u*hi,u,beta);ok=Klo>0&&Khi<0;if ~ok,return,end
for j=1:200,mid=(lo+hi)/2;if mid==lo||mid==hi,break,end;km=K(u*mid,u,beta);if km>0,lo=mid;elseif km<0,hi=mid;else,break,end,end
Klo=K(u*lo,u,beta);Khi=K(u*hi,u,beta);ok=Klo>0&&Khi<0;
end
function y=M(a,u,beta),y=(a^beta+(u-a)^beta)/(a^beta+(1-a)^beta);end
function y=K(a,u,beta),m=a^beta+(u-a)^beta;d=a^beta+(1-a)^beta;mp=beta*(a^(beta-1)-(u-a)^(beta-1));dp=beta*(a^(beta-1)-(1-a)^(beta-1));y=mp*d-m*dp;end
