function out=sm_normalizer64(Ns,trials,workers,outdir)
% Conditional mixture importance sampling, not certified bounds.
if nargin<1,Ns=[2 256 4096 65536 1048576];end
if nargin<2,trials=[8192 4096 4096 4096 2048];end
if nargin<3,workers=12;end
if nargin<4,outdir=fullfile(tempdir,'strong-maximal-r064');end
assert(numel(Ns)==numel(trials));
if ~isfolder(outdir),mkdir(outdir);end
pool=parpool('threads',workers);cleanup=onCleanup(@()delete(pool));mu=fzero(@(x)x-1-exp(-x),[1 2]);pmu=1/(1+exp(mu));results=cell(size(Ns));[gx16,gw16]=gauss64(16);[gx8,gw8]=gauss64(8);start=tic;
for ci=1:numel(Ns)
 b=Ns(ci);L=log(b);
 if b==2,delta=1/20;s=3*mu;else,delta=sqrt(L/b);z=fzero(@(z)erfc(z/sqrt(2))/2-sqrt(mu/b)*z,[0 sqrt(L)+2]);s=b*mu-sqrt(mu*b)*z;end
 a=mu-delta;d=mu+delta;rad=sqrt(110*b);ks=(max(0,floor(b/(1+exp(d))-rad)):min(b,ceil(b/(1+exp(a))+rad)))';choose=gammaln(b+1)-gammaln(ks+1)-gammaln(b-ks+1);cgFine=linspace(a,d,257);lfFine=zeros(size(cgFine));
 for j=1:numel(cgFine)
 cc=cgFine(j);g=s-cc*ks;lp=choose-cc*ks-b*log1p(exp(-cc))+(b-1)*log(max(g,realmin))-max(g,0)-gammaln(b);lp(g<=0)=-Inf;mx=max(lp);lfFine(j)=mx+log(sum(exp(lp-mx)));
 end
 cg=cgFine(1:2:end);lf=lfFine(1:2:end);interpolationError=max(abs(interp1(cg,lf,cgFine,'pchip')-lfFine));
 nr=trials(ci);vals=zeros(nr,1);quadDiff=vals;sumResidual=vals;
 parfor it=1:nr
 stream=RandStream('mt19937ar','Seed',640064+ci*100000+it);cc=a+(d-a)*rand(stream);g=s-cc*ks;lp=choose-cc*ks+(b-1)*log(max(g,realmin))-max(g,0)-gammaln(b);lp(g<=0)=-Inf;ww=exp(lp-max(lp));cdf=cumsum(ww);kk=ks(find(cdf>=rand(stream)*cdf(end),1));
 E=-log(max(rand(stream,b,1),realmin));zz=(s-cc*kk)*E/sum(E);zz(1:kk)=zz(1:kk)+cc;inside=sort(zz(zz>a & zz<d));v=[a;inside;d];vals(it)=ratio64(v,b,pmu,mu,cg,lf,gx16,gw16);
 if it<=32,other=ratio64(v,b,pmu,mu,cg,lf,gx8,gw8);quadDiff(it)=abs(other/vals(it)-1);sumResidual(it)=abs(sum(zz)-s);end
 end
 r=struct('b',b,'trials',nr,'a',a,'d',d,'s',s,'estimate',mean(vals),'standard_error',std(vals)/sqrt(nr),'ratio_over_sqrt_b_logb',mean(vals)/sqrt(b*L),'max_8_16_node_relative_difference',max(quadDiff),'max_sample_sum_residual',max(sumResidual),'seed_rule','640064 + case_index*100000 + trial_index; mt19937ar','density_interpolation_nodes',129,'max_129_257_logdensity_interpolation_difference',interpolationError,'quadrature_nodes',16,'mode','Monte Carlo, not confidence-certified');
 r.exact_b2_normalizer=NaN;
 if b==2,t=(d-a)/(s-d);r.exact_b2_normalizer=1-t+2*log1p(t);end
 r.envelope_upper=1+b*log((1+exp(-a))/(1+exp(-d)))-(lf(end)-lf(1));results{ci}=r;save64(fullfile(outdir,sprintf('normalizer_%d.mat',b)),r);fprintf('b=%d: Z estimate %.8g +/- SE %.3g\n',b,r.estimate,r.standard_error);
end
out=struct('status','complete','version',version,'workers',pool.NumWorkers,'seconds',toc(start),'rows',[results{:}],'scope','Conditional likelihood-envelope normalizer from mixture importance sampling; not certified bounds');fid=fopen(fullfile(outdir,'normalizer64.json'),'w');fwrite(fid,jsonencode(out));fclose(fid);clear cleanup
end
function [x,w]=gauss64(n)
 k=(1:n-1)';beta=k./sqrt(4*k.^2-1);[v,d]=eig(diag(beta,1)+diag(beta,-1));[x,ix]=sort(diag(d));w=2*v(1,ix).^2;x=x';
end
function R=ratio64(v,b,pmu,mu,cg,lf,gx,gw)
 m=numel(v)-2;cand=[v(1);v(2:end-1);v(end)];j=[0;(0:m-1)';m];lc=-j*log(2)-b*log1p(pmu*expm1(-(cand-mu)))-interp1(cg,lf,cand,'pchip');mx=max(lc);mid=(v(1:end-1)+v(2:end))/2;half=(v(2:end)-v(1:end-1))/2;cq=mid+half.*gx;
 lq=-(0:m)'*log(2)-b*log1p(pmu*expm1(-(cq-mu)))-reshape(interp1(cg,lf,cq(:),'pchip'),size(cq));area=sum(half.*sum(exp(lq-mx).*gw,2));R=(v(end)-v(1))/area;
end
function save64(p,r),save(p,'r');end
