
function result=sm_upper_window62(Ns,offsets,workers,outdir,chunk,tailExponent,derivativeFormula)
% Reproducible finite Gamma--binomial window checks; floating evidence only.
% Batch exits automatically when invoked with MATLAB -batch.
if nargin<1, Ns=2.^(5:12); end
if nargin<2, offsets=-2:0.25:2; end
if nargin<3, workers=12; end
if nargin<4, outdir=fullfile(tempdir,'strong-maximal-r062'); end
if nargin<5, chunk=4; end
if nargin<6, tailExponent=0; end % zero retains full finite sum
assert(tailExponent>=0);
if nargin<7, derivativeFormula='centered'; end
assert(any(strcmp(derivativeFormula,{'centered','expanded'})));
assert(all(Ns>8 & Ns==fix(Ns)) && workers>=1 && chunk>=1);

if ~isfolder(outdir), mkdir(outdir); end
p=gcp('nocreate'); if ~isempty(p), delete(p); end
p=parpool('threads',workers);
cleanup62=onCleanup(@()delete(p));
mu=fzero(@(u)u-1-exp(-u),[1 2]);
shifts=-4:1:4;
zc=zeros(size(Ns));
for k=1:numel(Ns)
 zc(k)=fzero(@(z)erfc(z/sqrt(2))/2-sqrt(mu/Ns(k))*z,[0 sqrt(log(Ns(k)))+2]);
end
[ii,jj]=ndgrid(1:numel(Ns),1:numel(offsets));
tasks=[ii(:) jj(:)]; rows=cell(size(tasks,1),1); started=tic;
fprintf('Upper62: %d independent cases, %d thread workers\n',size(tasks,1),p.NumWorkers);
parfor aidx=1:size(tasks,1)
 ni=tasks(aidx,1); n=Ns(ni); L=log(n); rho=sqrt(L/n);
 c=mu+offsets(tasks(aidx,2))*rho;
 pc=1/(1+exp(c));
 if tailExponent>0
  radius=sqrt(tailExponent*n/2);
  kval=(max(0,floor(n*pc-radius)):min(n,ceil(n*pc+radius)))';
  omittedMassBound=min(1,2*exp(-tailExponent));
 else
  kval=(0:n)'; omittedMassBound=0;
 end
 logw=gammaln(n+1)-gammaln(kval+1)-gammaln(n-kval+1)-c*kval-n*log1p(exp(-c));
 w=exp(logw); pc=1/(1+exp(c)); mc=1+c*pc; vc=1+c*c*pc*(1-pc);
 tq=n*mu-sqrt(mu*n)*zc(ni)+shifts*L;
 tc=n*mc+linspace(-2*sqrt(L),2*sqrt(L),257)*sqrt(n*vc);
 ts=[tq tc]; dens=zeros(size(ts)); deriv=dens; cdfs=dens; cscores=dens; residuals=dens;
 for left=1:chunk:numel(ts)
  ids=left:min(left+chunk-1,numel(ts)); zz=ts(ids)-c*kval; mask=zz>0;
  pdf=exp((n-1)*log(max(zz,realmin))-max(zz,0)-gammaln(n)); pdf(~mask)=0;
  ds=zeros(size(zz)); ds(mask)=(n-1)./zz(mask)-1;
  dens(ids)=w'*pdf; deriv(ids)=w'*(pdf.*ds);
  cdfs(ids)=w'*gammainc(max(zz,0),n,'lower');
  if strcmp(derivativeFormula,'expanded')
  inv1=zeros(size(zz)); inv1(mask)=pdf(mask)./zz(mask);
  e1=(w.*kval)'*inv1; inv2=zeros(size(zz)); inv2(mask)=inv1(mask)./zz(mask); e2=(w.*(kval.^2))'*inv2;
  cscores(ids)=n*pc-(n-1)*e1./dens(ids);
  second=(n*pc*((n+1)*pc-1)*dens(ids)-2*n*pc*(n-1)*e1+(n-1)*(n-2)*e2)./dens(ids);
  else
   % Equivalent component score: cancel O(n) terms before squaring;
   % this avoids subtracting O(n^2) aggregate terms.
   kg=kval./max(zz,realmin); kg(~mask)=0;
   aa=n*pc-(n-1)*kg;
   cscores(ids)=(w'*(pdf.*aa))./dens(ids);
   second=(w'*(pdf.*(aa.^2-n*pc*(1-pc)-(n-1)*kg.^2)))./dens(ids);
  end
  msecond=-2*pc*(1-pc)+c*pc*(1-pc)*(1-2*pc);
  residuals(ids)=second+n*msecond*deriv(ids)./dens(ids);
 end
 fq=dens(1:numel(tq)); Fc=cdfs(1:numel(tq));
 fc=dens(numel(tq)+1:end); sc=deriv(numel(tq)+1:end)./fc;
 r=struct('N',n,'binomial_support_size',numel(kval),'binomial_omitted_mass_bound',omittedMassBound,'offset',offsets(tasks(aidx,2)),'c',c,'weight_mass_residual',sum(w)-1,...
  'root_residual',erfc(zc(ni)/sqrt(2))/2-sqrt(mu/n)*zc(ni),...
  'critical_shifts',shifts,'N_density_over_logN',n*fq/L,'critical_CDF',Fc,...
  'central_score_over_rho_max',max(abs(sc))/rho,...
  'central_parameter_score_over_logN_max',max(abs(cscores(numel(tq)+1:end)))/L,'central_residual_over_logN2_max',max(abs(residuals(numel(tq)+1:end)))/(L*L),'central_fisher_integral',trapz(tc,fc.*sc.^2),'fisher_upper',1/(n-2),...
  'finite_checks',all(isfinite([fq Fc sc cscores residuals])) && all(fq>=0) && all(Fc>=0) && all(Fc<=1+1e-7));
 rows{aidx}=r;
 saveCase62(fullfile(outdir,sprintf('case_%03d.mat',aidx)),r);
end
result=struct('status','complete','matlab',version,'workers',p.NumWorkers,'case_count',numel(rows),...
 'seconds',toc(started),'derivative_formula',derivativeFormula,'mu',mu,'scope','Floating Gamma-Binomial kernels; no certified upper-bound conclusion','rows',[rows{:}]);
fid=fopen(fullfile(outdir,'matlab_upper62.json'),'w'); fwrite(fid,jsonencode(result)); fclose(fid);
fprintf('Completed %d cases in %.2f seconds; all finite=%d\n',numel(rows),result.seconds,all(cellfun(@(r)r.finite_checks,rows)));
clear cleanup62
end

function saveCase62(path62,r)
 save(path62,'r');
end
