function sm_far_actions66(outdir)
if nargin<1,outdir=fullfile(tempdir,'sm66');end
if ~isfolder(outdir),mkdir(outdir);end
base=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'..','rounds','r065','matlab65.json')));rr=base.rows;rr=rr(ismember([rr.profile],[1 3 5]));mu=fzero(@(c)c-1-exp(-c),[1 2]);
Ns=unique([rr.N]);infos=cell(size(Ns));for j=1:numel(Ns),n=Ns(j);kk=(0:n)';infos{j}=struct('n',n,'choose',gammaln(n+1)-gammaln(kk+1)-gammaln(n-kk+1));end
pool=parpool('threads',12);cleanup=onCleanup(@()delete(pool));rows=cell(size(rr));start=tic;
parfor j=1:numel(rr)
 R=rr(j); checkpoint=fullfile(outdir,sprintf('case_%d_%d.json',R.N,R.profile)); if isfile(checkpoint),rows{j}=jsondecode(fileread(checkpoint));continue;end; I=infos{find(Ns==R.N,1)};rho=sqrt(log(R.N)/R.N);cg=unique([linspace(0,4,161),logspace(log10(4),4,30),max(0,mu+linspace(-8,8,129)*rho)]);vals=zeros(size(cg));mass=zeros(size(cg));
 for k=1:numel(cg),[vals(k),mass(k)]=mean66(cg(k),I,R.q,R.slope);end
 [v,ix]=max(vals);roots=[];for k=1:numel(cg)-1,if vals(k)*vals(k+1)<0,roots(end+1)=fzero(@(c)mean66(c,I,R.q,R.slope),[cg(k),cg(k+1)],optimset('TolX',1e-10));end,end
 far=zeros(1,3);for t=1:3,mask=abs(cg-mu)>2^(t-1)*rho;far(t)=max(vals(mask));end
 rows{j}=struct('N',R.N,'profile',R.profile,'slope',R.slope,'q',R.q,'grid_count',numel(cg),'c_range',[min(cg),max(cg)],'max_mean_minus_one',v,'c_grid_best',cg(ix),'grid_best_shift_over_rho',(cg(ix)-mu)/rho,'roots_shift_over_rho',(roots-mu)/rho,'max_mean_outside_radii_1_2_4_rho',far,'max_prior_mass_residual',max(abs(mass)),'all_finite',all(isfinite(vals))); write66(checkpoint,rows{j});
end
result=struct('status','complete','version',version,'workers',pool.NumWorkers,'seconds',toc(start),'rows',[rows{:}],'scope','Frozen bounded critical ramp family, finite full-action grid and sign-change roots only, not global optimization certificate');f=fopen(fullfile(outdir,'far_actions66.json'),'w');fwrite(f,jsonencode(result));fclose(f);clear cleanup;fprintf('Completed far-action scan: %d frozen profiles\n',numel(rows));
end
function [v,res]=mean66(c,I,q,alpha)
n=I.n;p=1/(1+exp(c));k=(max(0,floor(n*p-sqrt(110*n))):min(n,ceil(n*p+sqrt(110*n))))';w=exp(I.choose(k+1)-c*k-n*log1p(exp(-c)));raw=q-c*k;raw2=q+n/alpha-c*k;g=max(raw,0);g2=max(raw2,0);F=gammainc(g,n,'lower');Q=gammainc(g,n,'upper');Q2=gammainc(g2,n,'upper');pos=n*gammainc(g,n+1,'upper')-raw.*Q;pos2=n*gammainc(g2,n+1,'upper')-raw2.*Q2;v=sum(w.*(-F+alpha/n*(pos-pos2)));res=sum(w)-1;
end

function write66(path,row)
f=fopen(path,'w');fwrite(f,jsonencode(row));fclose(f);
end
