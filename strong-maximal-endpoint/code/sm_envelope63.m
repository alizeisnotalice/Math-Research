function out=sm_envelope63(Ns,offs,workers,outdir,nodes)
% Numerical envelope moments, not interval-certified.
if nargin<1, Ns=2.^(8:20); end
if nargin<2, offs=[-1 0 1]; end
if nargin<3, workers=12; end
if nargin<4, outdir=fullfile(tempdir,'strong-maximal-r063'); end
if nargin<5, nodes=1025; end
if ~isfolder(outdir), mkdir(outdir); end
p=parpool('threads',workers); cleanup63=onCleanup(@()delete(p));
mu=fzero(@(u)u-1-exp(-u),[1 2]);
[ii,jj]=ndgrid(1:numel(Ns),1:numel(offs)); tasks=[ii(:),jj(:)]; rows=cell(size(tasks,1),1); start=tic;
parfor ix=1:size(tasks,1)
 n=Ns(tasks(ix,1)); off=offs(tasks(ix,2)); L=log(n); rho=sqrt(L/n); c=mu+off*rho; pc=1/(1+exp(c)); mc=1+c*pc; vc=1+c*c*pc*(1-pc);
 z=fzero(@(z)erfc(z/sqrt(2))/2-sqrt(mu/n)*z,[0 sqrt(L)+2]); d=sqrt(mu*n)*z;s=n*mu-d;
 rad=sqrt(110*n); kval=(max(0,floor(n*pc-rad)):min(n,ceil(n*pc+rad)))';
 w=exp(gammaln(n+1)-gammaln(kval+1)-gammaln(n-kval+1)-c*kval-n*log1p(exp(-c)));
 ts=n*mc+linspace(-3*sqrt(L),3*sqrt(L),nodes)*sqrt(n*vc); f=zeros(size(ts)); df=f;
 for left=1:4:numel(ts)
 ids=left:min(left+3,numel(ts)); g=ts(ids)-c*kval; good=g>0; pdf=exp((n-1)*log(max(g,realmin))-max(g,0)-gammaln(n));pdf(~good)=0;
 ds=zeros(size(g));ds(good)=(n-1)./g(good)-1;f(ids)=w'*pdf;df(ids)=w'*(pdf.*ds);
 end
 delta=max(ts-s+20*L,0); v=min(delta/n,log(2));
 for iter=1:12, v=max(0,min(log(2),v-(v+(n-1)*log1p(v)-delta)./(1+(n-1)./(1+v))));end
 coreCDF=w'*gammainc(max(s-20*L-c*kval,0),n,'lower');
 U=expm1(v); EU=trapz(ts,f.*U); EU2=trapz(ts,f.*U.^2); score=trapz(ts,U.*df);
 rr=struct('N',n,'offset',off,'nodes',nodes,'K',20,'cap',2,'core_CDF',coreCDF,'envelope_core_gap_over_L2N',(EU-coreCDF)*n/L^2,'EU_over_dN',EU/(d/n),'EU2_over_dN2',EU2/(d/n)^2,'E_U_location_score',score,'mass',trapz(ts,f),'tail_binomial_mass_bound',2*exp(-220),'inverse_log_equation_residual',max(abs(v(v<log(2))+(n-1)*log1p(v(v<log(2)))-delta(v<log(2)))),'finite',all(isfinite([EU EU2 score])));
 rows{ix}=rr; save63(fullfile(outdir,sprintf('case_%03d.mat',ix)),rr);
end
out=struct('status','complete','workers',p.NumWorkers,'version',version,'seconds',toc(start),'rows',[rows{:}],'scope','Trapezoidal moments of scalar upper envelope under Gamma-binomial kernels; numerical evidence only');fid=fopen(fullfile(outdir,'envelope63.json'),'w');fwrite(fid,jsonencode(out));fclose(fid);fprintf('Completed %d envelope cases in %.3f seconds\n',numel(rows),out.seconds);clear cleanup63
end
function save63(p,rr),save(p,'rr');end
