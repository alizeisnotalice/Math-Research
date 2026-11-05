function result=sm_quantile_derivative_checks
B=fzero(@(b)log(b)-1-1/b,[3 4]); mu=log(B); sigma=sqrt(mu); ps=(mu-1)/mu;
rows=cell(0,1);
for n=[256 1024 4096 1e4 1e5 1e6 1e7]
 L=log(n); rho=sqrt(L/n); beta=sigma*rho; cs=mu-sigma*rho/2;
 [q0,omit0]=quant(n,mu,beta,mu); [q1,omit1]=quant(n,cs,beta,mu);
 fprintf('n=%g explicit gain/L=%.10g limit=%.10g omitted<~%.3g\n',n,(q1-q0)/L,(mu-1)/8,max(omit0,omit1));
 row=struct('n',n,'q_mu',q0,'q_shifted',q1,'gain_over_log_n',(q1-q0)/L, ...
 'shifted_c',cs,'omitted_mass_estimate',max(omit0,omit1));
 if n<=4096
  [~,f0]=dist(n,mu,q0); F=dist(n,mu,q0); Fprev=dist(n-1,mu,q0-mu);
  exact=n*ps*(Fprev-F)/f0; hh=1e-4;
  fp=dist(n,mu+hh,q0); fm=dist(n,mu-hh,q0);
  row.quantile_derivative=exact;
  row.cdf_derivative_identity_residual=(fp-fm)/(2*hh)-n*ps*(F-Fprev);
  row.derivative_over_sqrt_n_log_n=exact/sqrt(n*L);
 end
 rows{end+1}=row; %#ok<AGROW>
end
result=struct('matlab_version',version,'target_gain_coefficient',(mu-1)/8, ...
 'method','12 sd binomial truncation; exact gamma CDF; 55 quantile bisections', ...
 'status','double-precision diagnostics, not certificates','rows',{rows});
root=fileparts(fileparts(mfilename('fullpath')));
fid=fopen(fullfile(root,'rounds','r016','checks.json'),'w');
assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
function [q,omit]=quant(n,c,beta,mu)
lo=0;hi=n*mu+12*sqrt(n);
for j=1:55
 mid=(lo+hi)/2; val=dist(n,c,mid); if val<beta,lo=mid;else,hi=mid;end
end
q=(lo+hi)/2;[~,~,omit]=dist(n,c,q);
end
function [F,f,omit]=dist(n,c,q)
p=1/(1+exp(c)); sd=sqrt(n*p*(1-p)); center=floor((n+1)*p);
k=max(0,floor(n*p-12*sd)):min(n,ceil(n*p+12*sd)); pivot=center-k(1)+1;
w=zeros(size(k));w(pivot)=1;
for j=pivot+1:numel(k),w(j)=w(j-1)*(n-k(j)+1)/k(j)*p/(1-p);end
for j=pivot-1:-1:1,w(j)=w(j+1)*(k(j)+1)/(n-k(j))*(1-p)/p;end
w=w/sum(w); x=max(q-c*k,0);F=sum(w.*gammainc(x,n,'lower'));
if nargout>1
 pos=x>0; dens=zeros(size(x));dens(pos)=exp((n-1)*log(x(pos))-x(pos)-gammaln(n));f=sum(w.*dens);
end
if nargout>2,omit=binocdf(k(1)-1,n,p)+binocdf(k(end),n,p,'upper');end
end
