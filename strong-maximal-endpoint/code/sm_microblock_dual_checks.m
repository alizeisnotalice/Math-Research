function result=sm_microblock_dual_checks(gridN,out)
% Positive quadrature diagnostics for the capped microblock convex program.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,gridN=1024;end
if nargin<2,out=fullfile(root,'rounds','r026','checks.json');end
mu=fzero(@(x)x-1-exp(-x),[1 2]);rows={};tic;
for n=[2 4 8 16 32 64 128 256]
 c=mu-sqrt(mu)/2*sqrt(log(n)/n);
 for m=1:2
  b=n/m; N=gridN; if m==1,N=8*gridN;end
  p=1/(1+exp(c));meanZ=1+c*p;varZ=1+c*c*p*(1-p);
  cut=b*meanZ+16*sqrt(b*varZ);splits=[];if b<=4,splits=c*(0:b);end;[t,quad]=legendre(N,cut,splits);
  [lf,ell]=logdensity(b,c,t);pw=quad.*exp(lf);
  if m==2
   ell=ell+ell';pw=pw.*pw';
  end
  ell=ell(:);pw=pw(:);valid=isfinite(ell)&pw>0;ell=ell(valid);pw=pw(valid);
  lo=min(ell)-n;hi=max(ell);
  for it=1:60
   tau=(lo+hi)/2;[h,~]=profile(ell-tau,n);expect=sum(pw.*h);
   if expect>1,lo=tau;else,hi=tau;end
  end
  tau=lo;[h,~]=profile(ell-tau,n);expect=sum(pw.*h);
  % Resolve an atom of the quadrature measure at the 0/1 threshold.
  partial=false;
  if expect-1>1e-10
   [~,ix]=min(abs(ell-tau));tau=ell(ix);[h,~]=profile(ell-tau,n);
   tie=ell==tau;h(tie)=0; fraction=(1-sum(pw.*h))/sum(pw(tie));
   assert(fraction>=-1e-7&&fraction<=1+1e-7);h(tie)=min(max(fraction,0),1);partial=true;
  end
  expect=sum(pw.*h);assert(abs(expect-1)<1e-8);y=log(max(h,realmin));
  logPhi=log(max(h,realmin))+(n-1)*log1p(max(y,0));
  d=ell-tau; interior=h>1+1e-12 & h<2-1e-12;
  kktResidual=0;
  if any(interior)
   yi=y(interior);kktResidual=max(abs((n-2)*log1p(yi)+log(n+yi)-d(interior)));
  end
  assert(kktResidual<1e-8);
  k=0:b;lp=gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1)+k*log(p)+(b-k)*log1p(-p);
  blockTail=sum(exp(lp).*gammainc(max(cut-c*k,0),b,'upper'));
  omittedProbabilityBound=m*blockTail;
  active=h>0;z=log(pw(active))-ell(active)+logPhi(active);peak=max(z);logJ=peak+log(sum(exp(z-peak)));
  rows{end+1}=struct('n',n,'m',m,'c',c,'nodes_per_axis',numel(t),'cut',cut,'probability_mass',sum(pw),'tau',tau,'expectation',expect,'partial_boundary_atom',partial,'log_score',-logJ,'max_height_on_grid',max(h),'height_cap',2,'inverse_derivative_residual',kktResidual,'omitted_probability_bound',omittedProbabilityBound,'omitted_modular_relative_bound',2*omittedProbabilityBound*exp(-tau-logJ),'solver','60 multiplier bisections; seven inverse-derivative Newton steps; quadrature only'); %#ok<AGROW>
  fprintf('n%d m%d N%d logScore%.10f mass%.15g rawExpect%.12g maxH%.6g sec%.2f\n',n,m,N,-logJ,sum(pw),expect,max(h),toc);
 end
end
result=struct('matlab',version,'gridN',gridN,'status','positive quadrature diagnostics, not directed continuum certificate; analytic tail formulas evaluated in double','seconds',toc,'rows',{rows});fid=fopen(out,'w');assert(fid>=0);cl=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
function [h,y]=profile(d,n)
h=zeros(size(d));y=zeros(size(d));active=d>=0;h(active)=1;soft=d>log(n);
y(soft)=max(expm1((d(soft)-log(n))/(n-1)),0);
for k=1:7
 ys=y(soft);res=(n-2)*log1p(ys)+log(n+ys)-d(soft);der=(n-2)./(1+ys)+1./(n+ys);
 y(soft)=max(0,ys-res./der);
end
y=min(y,log(2));h(soft)=exp(y(soft));
end
function [lf,logR]=logdensity(b,c,t)
% e^t f_b(t)=sum binom(b,k)(t-ck)_+^(b-1)/[Gamma(b)(1+e^-c)^b].
k=0:b;lw=gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1);
x=t-c*k;lg=lw+(b-1)*log(max(x,realmin));lg(x<=0)=-Inf;
pk=max(lg,[],2);logR=pk+log(sum(exp(lg-pk),2))-b*log1p(exp(-c))-gammaln(b);
lf=logR-t;
end
function [t,w]=legendre(N,cut,splits)
% Composite GL16 avoids a large dense eigenproblem.
q=16;ii=(1:q-1)';beta=ii./sqrt(4*ii.^2-1);[V,D]=eig(diag(beta,1)+diag(beta,-1));[x,ix]=sort(diag(D));v=2*(V(1,ix)').^2;panels=ceil(N/q);edges=unique([linspace(0,cut,panels+1),splits]);edges=edges(edges>=0&edges<=cut);mid=(edges(1:end-1)+edges(2:end))/2;half=diff(edges)/2;t=x*half+mid;w=v*half;t=t(:);w=w(:);
end
