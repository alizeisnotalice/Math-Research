function sm_actions65(outdir,workers)
if nargin<1,outdir=fullfile(tempdir,'strong-maximal-r065');end
if ~isfolder(outdir),mkdir(outdir);end
if nargin<2,workers=12;end
pool=parpool('threads',workers);clean=onCleanup(@()delete(pool));mu=fzero(@(x)x-1-exp(-x),[1 2]);Ns=2.^(8:2:20);infos=cell(size(Ns));
for j=1:numel(Ns)
 n=Ns(j);L=log(n);rho=sqrt(L/n);ca=mu-2*rho;cd=mu+2*rho;rad=sqrt(110*n);ks=(max(0,floor(n/(1+exp(cd))-rad)):min(n,ceil(n/(1+exp(ca))+rad)))';z=fzero(@(z)erfc(z/sqrt(2))/2-sqrt(mu/n)*z,[0 sqrt(L)+2]);infos{j}=struct('n',n,'L',L,'rho',rho,'mu',mu,'a',ca,'d',cd,'q0',n*mu-sqrt(mu*n)*z,'ks',ks,'choose',gammaln(n+1)-gammaln(ks+1)-gammaln(n-ks+1),'lgamma',gammaln(n));
end
[ii,jj]=ndgrid(1:numel(Ns),1:6);tasks=[ii(:),jj(:)];rows=cell(size(tasks,1),1);start=tic;
parfor ix=1:size(tasks,1)
 I=infos{tasks(ix,1)};profile=tasks(ix,2);kk=profile-3;alpha=1+kk*I.rho;
 if profile==6,alpha=1;q=I.q0;kind='canonical-reference';else,q=fzero(@(q)pick65(I.mu,I,q,alpha,1),[I.q0-20*I.L,I.q0+20*I.L],optimset('TolX',1e-10));kind='mean-one-frozen-ramp';end
 cg=linspace(I.a,I.d,13);rr=zeros(numel(cg),4);for k=1:numel(cg),rr(k,:)=response65(cg(k),I,q,alpha);end
 candidates=[I.a,I.d];for k=1:numel(cg)-1,if rr(k,2)*rr(k+1,2)<0,candidates(end+1)=fzero(@(c)pick65(c,I,q,alpha,2),[cg(k),cg(k+1)],optimset('TolX',1e-11));end,end
 vals=zeros(size(candidates));for k=1:numel(candidates),vals(k)=pick65(candidates(k),I,q,alpha,1);end;[best,jbest]=max(vals);cc=candidates(jbest);at=response65(cc,I,q,alpha);atmu=response65(I.mu,I,q,alpha);
 limit=min(sqrt(I.n)/alpha,20+4*sqrt(I.L));integrand=@(v)exp(log1p(alpha*v/sqrt(I.n))+(I.n-1)*log1p(log1p(alpha*v/sqrt(I.n)))-v*sqrt(I.n));[integral,err]=quadgk(integrand,0,limit,'RelTol',1e-11,'AbsTol',1e-12);budgetI=sqrt(I.n)*integral+exp(log(2)+(I.n-1)*log1p(log(2))-I.n/alpha);s=q-log(budgetI);
 r=struct('N',I.n,'profile',profile,'kind',kind,'slope',alpha,'q',q,'mean_at_mu_minus_one',atmu(1),'c_candidate',cc,'c_shift_over_rho',(cc-I.mu)/I.rho,'mean_at_candidate_minus_one',best,'gradient_at_candidate',at(2),'curvature_at_candidate',at(3),'grid_curvature_max',max(rr(:,3)),'root_candidate_count',numel(candidates)-2,'grid_gradient_increasing_violations',sum(diff(rr(:,2))>1e-9),'max_prior_weight_mass_residual',max(abs(rr(:,4))),'budget_s',s,'budget_s_slack_over_logN',(s-I.q0)/I.L,'budget_quadrature_error_estimate',sqrt(I.n)*err,'finite',all(isfinite([cc q s at])));
 rows{ix}=r;save65(fullfile(outdir,sprintf('case_%03d.mat',ix)),r);
end
out=struct('status','complete','version',version,'workers',pool.NumWorkers,'seconds',toc(start),'rows',[rows{:}],'scope','Numerical constrained action candidates for frozen capped ramps; no interval/global optimizer certificate');fid=fopen(fullfile(outdir,'actions65.json'),'w');fwrite(fid,jsonencode(out));fclose(fid);fprintf('Completed %d frozen-profile cases in %.3f seconds\n',numel(rows),out.seconds);clear clean
end
function v=pick65(c,I,q,alpha,k),v=response65(c,I,q,alpha);v=v(k);end
function out=response65(c,I,q,alpha)
 n=I.n;k=I.ks;p=1/(1+exp(c));w=exp(I.choose-c*k-n*log1p(exp(-c)));raw=q-c*k;raw2=q+n/alpha-c*k;g=max(raw,0);g2=max(raw2,0);F=gammainc(g,n,'lower');Q=gammainc(g,n,'upper');Q2=gammainc(g2,n,'upper');pos=n*gammainc(g,n+1,'upper')-raw.*Q;pos2=n*gammainc(g2,n+1,'upper')-raw2.*Q2;hm1=-F+alpha/n*(pos-pos2);
 pdf=exp((n-1)*log(max(g,realmin))-g-I.lgamma);pdf(g==0)=0;pdf2=exp((n-1)*log(max(g2,realmin))-g2-I.lgamma);pdf2(g2==0)=0;dp=zeros(size(g));mask=g>0;dp(mask)=(n-1)./g(mask)-1;
 aa=n*p-k;ac=k.*pdf+alpha*k/n.*(Q-Q2);acc=-(k.^2).*pdf.*dp+alpha*(k.^2)/n.*(pdf-pdf2);
 out=[sum(w.*hm1),sum(w.*(aa.*hm1+ac)),sum(w.*((aa.^2-n*p*(1-p)).*hm1+2*aa.*ac+acc)),sum(w)-1];
end
function save65(p,r),save(p,'r');end
