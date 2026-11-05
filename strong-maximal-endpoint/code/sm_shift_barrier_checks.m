function result=sm_shift_barrier_checks
B=fzero(@(b)log(b)-1-1/b,[3 4]);mu=log(B);sig=sqrt(mu);K=100;C0=1e10;A=10;
pde=cell(0,1);
for ell=[10 100 1000 1e6]
 smallest=Inf;
 for ratio=[1 2 100 1e6]
  tt=C0*ell*ratio;zz=linspace(-A*sqrt(ell),A*sqrt(ell),401);xx=sig*sqrt(tt)*zz;
  hh=K*(zz.^2+2*log(ratio));
  coef=-K/tt-K*xx.^2/(mu*tt^2)+2*K*xx.*hh/(mu*tt^2)+2*K^2*xx.^2.*(xx-hh)/(mu^2*tt^3);
  smallest=min(smallest,min(-coef*tt./(K*(1+zz.^2))));
 end
 pde{end+1}=struct('log_horizon',ell,'minimum_normalized_heat_margin',smallest); %#ok<AGROW>
 fprintf('L=%g nonlinear PDE margin=%.12g\n',ell,smallest);
end
[ys64,ws64]=laguerre(64);[ys128,ws128]=laguerre(128);
rows=cell(0,1);
for ell=[10 20 40]
 t0=C0*ell;tt=t0+1;
 for z=[-1 0 1]*sqrt(ell)
  xx=sig*sqrt(tt)*z;
  fun=@(c)step(c,ys64,ws64,xx,tt,t0,ell,K,A,mu);
  [cc,val,flag]=fminbnd(fun,.05,5,optimset('TolX',1e-8));
  [val,idx]=min([val fun(0) fun(50) fun(Inf)]);cs=[cc 0 50 Inf];cc=cs(idx);
  fine=step(cc,ys128,ws128,xx,tt,t0,ell,K,A,mu);
  fixed=step(mu,ys128,ws128,xx,tt,t0,ell,K,A,mu);
  rows{end+1}=struct('log_horizon',ell,'z',z,'candidate_c',cc, ...
   'normalized_step_margin_64',val,'normalized_step_margin_128',fine, ...
   'fixed_mu_margin_128',fixed,'quadrature_difference',fine-val,'exitflag',flag); %#ok<AGROW>
  fprintf('L=%g z=%.5g step margin=%.9g quadratureDelta=%.3g\n',ell,z,fine,fine-val);
 end
end
result=struct('matlab_version',version,'C0',C0,'K',K,'A',A, ...
 'method','exact PDE; Gauss-Laguerre 64/128; 4-point Gaussian integral for stable CDF difference', ...
 'status','floating-point diagnostics, not uniform certificates','pde',{pde},'step_rows',{rows});
root=fileparts(fileparts(mfilename('fullpath')));fid=fopen(fullfile(root,'rounds','r021','checks.json'),'w');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
function out=step(c,ys,weights,x,t,t0,ell,K,A,mu)
sig=sqrt(mu);hc=height(x,t);yc=x-hc;u=yc/(sig*sqrt(t));
if isinf(c),p=0;else,p=1/(1+exp(c));end
out=(1-p)*sum(weights.*difference(mu-ys));
if p>0,out=out+p*sum(weights.*difference(mu-ys-c));end
 function h=height(xx,tt)
  zz=abs(xx)/(sig*sqrt(tt));aa=A*sqrt(ell);ss=zz/aa-2;
  cap=zz.^2;transition=ss>0 & ss<1;v=ss(transition);
  cap(transition)=aa^2*(4+v.^2+4*v-10*v.^4+8*v.^5+v.^6-12*v.^7/7);
  cap(ss>=1)=44*aa^2/7;
  h=K*(cap+2*log(tt/t0));
 end
 function d=difference(shift)
  hp=height(x+shift,t-1);
  du=(shift-(hp-hc))/(sig*sqrt(t-1))+yc/sig/(sqrt(t)*sqrt(t-1)*(sqrt(t)+sqrt(t-1)));
  nodes=[.0694318442029737 .330009478207572 .669990521792428 .930568155797026];
  ww=[.173927422568727 .326072577431273 .326072577431273 .173927422568727];
  d=zeros(size(du));
  for j=1:4,vv=du*nodes(j);d=d+ww(j)*exp(-u*vv-vv.^2/2);end
  d=sig*sqrt(t)*du.*d;
 end
end
function [nodes,weights]=laguerre(n)
j=(1:n)';matrix=diag(2*j-1)+diag((1:n-1)',1)+diag((1:n-1)',-1);
[v,d]=eig(matrix,'vector');[nodes,idx]=sort(d);weights=v(1,idx)'.^2;
end
