function result=sm_soft_shift_probe
B=fzero(@(b)log(b)-1-1/b,[3 4]);mu=log(B);sigma=sqrt(mu);rows=cell(0,1);
for n=[256 1024 4096 1e4 1e5 1e6 1e7]
 rho=sqrt(log(n)/n);c=mu-sigma*rho/2;
 q0=softroot(n,mu,mu);q1=softroot(n,c,mu);
 logI=.5*log(2*pi*n)-1-n*log1p(-1/n)+log(gammainc(n-1,n,'upper'))+1/(12*n)-1/(360*n^3);
 fprintf('n=%g soft gain/L=%.10g target=%.10g logI=%.9g\n',n,(q1-q0)/log(n),(mu-1)/8,logI);
 rows{end+1}=struct('n',n,'fixed_mu_q',q0,'shifted_q',q1,'shifted_c',c, ...
 'gain_over_log_n',(q1-q0)/log(n),'log_I_Stirling',logI,'shifted_log_score',q1-logI, ...
 'fixed_log_average_residual',logaverage(n,mu,q0),'shifted_log_average_residual',logaverage(n,c,q1)); %#ok<AGROW>
end
result=struct('matlab_version',version,'status','numerical boundary averages; strictification/truncation and proof still required','rows',{rows});
fid=fopen('/Users/zhengzhihao/Desktop/T/work/soft_shift_probe.json','w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
function q=softroot(n,c,mu)
lo=0;hi=n*mu+12*sqrt(n);
for j=1:55
 mid=(lo+hi)/2;v=logaverage(n,c,mid);if v>0,lo=mid;else,hi=mid;end
end
q=(lo+hi)/2;
end
function a=logaverage(n,c,q)
p=1/(1+exp(c));eta=1/n;rate=1-eta;
logmgf=-log1p(-eta)+log1p(p*expm1(eta*c));pt=1/(1+exp(rate*c));sd=sqrt(n*pt*(1-pt));
k=max(0,floor(n*pt-12*sd)):min(n,ceil(n*pt+12*sd));center=floor((n+1)*pt);pivot=center-k(1)+1;
w=zeros(size(k));w(pivot)=1;
for j=pivot+1:numel(k),w(j)=w(j-1)*(n-k(j)+1)/k(j)*pt/(1-pt);end
for j=pivot-1:-1:1,w(j)=w(j+1)*(k(j)+1)/(n-k(j))*(1-pt)/pt;end
w=w/sum(w);failure=sum(w.*gammainc(rate*max(q-c*k,0),n,'lower'));
a=n*logmgf-q/n+log1p(-failure);
end
