function out=sm_ramp63(Ns,offs,workers,outdir)
% Actual capped ramp checks with finite Gamma--binomial mixtures.
if nargin<1, Ns=2.^(8:20);end
if nargin<2, offs=[-1 0 1];end
if nargin<3, workers=12;end
if nargin<4, outdir=fullfile(tempdir,'strong-maximal-r063');end
if ~isfolder(outdir),mkdir(outdir);end
p=parpool('threads',workers);clean=onCleanup(@()delete(p));
mu=fzero(@(x)x-1-exp(-x),[1 2]);[ii,jj]=ndgrid(1:numel(Ns),1:numel(offs));tasks=[ii(:),jj(:)];rows=cell(size(tasks,1),1);started=tic;
parfor k=1:size(tasks,1)
 n=Ns(tasks(k,1));off=offs(tasks(k,2));L=log(n);c=mu+off*sqrt(L/n);pc=1/(1+exp(c));rad=sqrt(110*n);ks=(max(0,floor(n*pc-rad)):min(n,ceil(n*pc+rad)))';w=exp(gammaln(n+1)-gammaln(ks+1)-gammaln(n-ks+1)-c*ks-n*log1p(exp(-c)));
 z=fzero(@(z)erfc(z/sqrt(2))/2-sqrt(mu/n)*z,[0 sqrt(L)+2]);d=sqrt(mu*n)*z;
 q=fzero(@(q)balance63(q,n,c,ks,w),[n*mu-d-20*L,n*mu-d+20*L]);
 integrand=@(v)exp(log1p(v/sqrt(n))+(n-1)*log1p(log1p(v/sqrt(n)))-v*sqrt(n));
 [ival,err]=quadgk(integrand,0,min(sqrt(n),20),'RelTol',1e-11,'AbsTol',1e-12);I=sqrt(n)*ival+exp(log(2)+(n-1)*log1p(log(2))-n);s=q-log(I);
 gq=max(q-c*ks,realmin);gs=max(s-c*ks,realmin);fq=sum(w.*exp((n-1)*log(gq)-gq-gammaln(n)));fs=sum(w.*exp((n-1)*log(gs)-gs-gammaln(n)));Fq=sum(w.*gammainc(gq,n,'lower'));Qt=sum(w.*gammainc(max(q+n-c*ks,0),n,'upper'));
 raw=q-c*ks;raw2=q+n-c*ks;gg=max(raw,0);gg2=max(raw2,0);
 Q=gammainc(gg,n,'upper');Q2=gammainc(gg2,n,'upper');
 Pq=exp((n-1)*log(max(gg,realmin))-gg-gammaln(n));Pq(gg==0)=0;
 P2=exp((n-1)*log(max(gg2,realmin))-gg2-gammaln(n));P2(gg2==0)=0;
 dp=zeros(size(gg)); good=gg>0;dp(good)=(n-1)./gg(good)-1;
 pos=n*gammainc(gg,n+1,'upper')-raw.*Q;
 pos2=n*gammainc(gg2,n+1,'upper')-raw2.*Q2;
 hv=Q+(pos-pos2)/n;aa=n*pc-ks;
 ac=ks.*Pq+ks/n.*(Q-Q2);
 acc=-(ks.^2).*Pq.*dp+(ks.^2)/n.*(Pq-P2);
 curvature=sum(w.*((aa.^2-n*pc*(1-pc)).*hv+2*aa.*ac+acc));
 normalBracket=z*exp(-z*z/2)/sqrt(2*pi)/(erfc(z/sqrt(2))/2)+1;
 rr=struct('common_curvature',curvature,'normal_curvature_bracket',normalBracket,'curvature_normal_residual',curvature+(mu-1)/mu*normalBracket,'N',n,'offset',off,'q',q,'s',s,'canonical_s_slack_over_logN',(s-(n*mu-d))/L,'mean_residual',balance63(q,n,c,ks,w),'N_location_response',n*fq+1-Fq-Qt,'N_density_at_actual_budget',n*fs,'refined_response_residual',n*(fq-fs)-Fq-Qt,'budget_quadrature_error_estimate',sqrt(n)*err,'finite',all(isfinite([fq fs q s])));
 rows{k}=rr;
end
out=struct('status','complete','workers',p.NumWorkers,'seconds',toc(started),'rows',[rows{:}],'scope','Actual capped scalar ramp with mean1; finite Gamma-binomial and floating quadrature, not interval certification');fid=fopen(fullfile(outdir,'ramp63.json'),'w');fwrite(fid,jsonencode(out));fclose(fid);fprintf('Ramp checks completed %.3f seconds\n',out.seconds);clear clean
end
function r=balance63(q,n,c,ks,w)
 raw=q-c*ks;raw2=q+n-c*ks;g=max(raw,0);g2=max(raw2,0);
 F=sum(w.*gammainc(g,n,'lower'));
 pos=sum(w.*(n*gammainc(g,n+1,'upper')-raw.*gammainc(g,n,'upper')));
 pos2=sum(w.*(n*gammainc(g2,n+1,'upper')-raw2.*gammainc(g2,n,'upper')));
 r=F-(pos-pos2)/n;
end
