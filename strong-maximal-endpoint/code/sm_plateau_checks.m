function result=sm_plateau_checks
root=fileparts(fileparts(mfilename('fullpath')));
old=jsondecode(fileread(fullfile(root,'rounds','r017','checks.json')));
B=fzero(@(b)log(b)-1-1/b,[3 4]);mu=log(B);rows=cell(0,1);
for r=old.rows'
 n=r.n;L=log(n);w=.5*L+log(L);p=1/(1+exp(mu));rate=1-1/n;
 [k0,w0]=weights(n,p);[kt,wt]=weights(n,1/(1+exp(rate*mu)));
 mgf=-log1p(-1/n)+log1p(p*expm1(mu/n));
 lo=r.fixed_mu_q-10;hi=r.fixed_mu_q+1;
 assert(averageMinusOne(lo)>0 && averageMinusOne(hi)<0);
 for j=1:40
  q=(lo+hi)/2;v=averageMinusOne(q);if v>0,lo=q;else,hi=q;end
 end
 q=(lo+hi)/2;factor=1-exp(-w)+exp(-w+r.log_I_Stirling);
 row=struct('n',n,'width',w,'plateau_root',q,'root_minus_soft',q-r.fixed_mu_q, ...
 'log_modular_factor',log(factor),'log_score_minus_soft',q-log(factor)-r.fixed_mu_q, ...
 'strict_threshold',q-1/n,'strict_average_margin',averageMinusOne(q-1/n), ...
 'root_average_residual',averageMinusOne(q));
 rows{end+1}=row; %#ok<AGROW>
 fprintf('n=%g plateau q-qsoft=%.9g score-qsoft=%.9g target=-.5\n',n,row.root_minus_soft,row.log_score_minus_soft);
end
result=struct('matlab_version',version,'status','double-precision platform root diagnostics; not directed certificates', ...
 'parameters',struct('binomial_sd_cutoff',12,'bisection_steps',40,'randomness','none'),'rows',{rows});
fid=fopen(fullfile(root,'rounds','r019','checks.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
function val=averageMinusOne(q)
 qq=q+w;ft=sum(wt.*gammainc(rate*max(qq-mu*kt,0),n,'lower'));
 logA=n*mgf-qq/n+log1p(-ft);
 F0=sum(w0.*gammainc(max(q-mu*k0,0),n,'lower'));
 F1=sum(w0.*gammainc(max(qq-mu*k0,0),n,'lower'));
 val=expm1(logA)+F1-F0;
end
end
function [k,w]=weights(n,p)
sd=sqrt(n*p*(1-p));k=max(0,floor(n*p-12*sd)):min(n,ceil(n*p+12*sd));pivot=floor((n+1)*p)-k(1)+1;
w=zeros(size(k));w(pivot)=1;
for j=pivot+1:numel(k),w(j)=w(j-1)*(n-k(j)+1)/k(j)*p/(1-p);end
for j=pivot-1:-1:1,w(j)=w(j+1)*(k(j)+1)/(n-k(j))*(1-p)/p;end
w=w/sum(w);
end
