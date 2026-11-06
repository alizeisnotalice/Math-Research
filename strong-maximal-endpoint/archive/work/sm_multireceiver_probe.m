function result=sm_multireceiver_probe(n,N,spacing,spanFactor,gamma,out)
% Convex finite-grid dual for symmetric multiple receivers, NOT interval proof.
root='/Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal/strong-maximal-endpoint';old=jsondecode(fileread(fullfile(root,'rounds','r026','checks.json')));ref=jsondecode(fileread(fullfile(root,'rounds','r027','checks.json')));
baseline=old.rows([old.rows.n]==n & [old.rows.m]==1);c=baseline.c;b=n/2;mu=fzero(@(x)x-1-exp(-x),[1 2]);A=2;
rr=ref.rows;idx=find([rr.n]==n & reshape(arrayfun(@(r)numel(r.blocks)==2&&r.blocks(1)==r.blocks(2),rr),1,[]),1);tauRef=rr(idx).tau;
K=ceil(spanFactor*sqrt(n)/spacing);v=(0:K)*spacing;d=gamma*v.^2/(mu*n);u=[0 0;[d(2:end)'+v(2:end)' d(2:end)'-v(2:end)'];[d(2:end)'-v(2:end)' d(2:end)'+v(2:end)']];M=K+1;coeff=[1;2*ones(K,1)];tic;
p=1/(1+exp(c));mz=1+c*p;vz=1+c*c*p*(1-p);lo=min(u,[],'all');hi=max(u,[],'all')+b*mz+16*sqrt(b*vz);[t,qw]=quad(N,lo,hi);[t1,t2]=ndgrid(t,t);[w1,w2]=ndgrid(qw,qw);sumt=t1(:)+t2(:);logQ=log(w1(:))+log(w2(:));logRef=logQ+tauRef-sumt;
R=zeros(numel(sumt),M);P=zeros(size(R));
for j=1:M
 if j==1,points=[0 0];else,points=[d(j)+v(j) d(j)-v(j);d(j)-v(j) d(j)+v(j)];end
 for k=1:size(points,1)
  a=logR(b,c,t-points(k,1));bb=logR(b,c,t-points(k,2));lr=a+bb'+sum(points(k,:));lr=lr(:);
  R(:,j)=R(:,j)+exp(lr-tauRef);P(:,j)=P(:,j)+exp(logQ+lr-sumt);
 end
end
masses=sum(P,1)';assert(max(abs(masses-coeff))<1e-6);
wStart=zeros(M,1);wStart(1)=1;wStart(2:end)=.001;
options=optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true,'Display','off','MaxIterations',200,'MaxFunctionEvaluations',1500,'OptimalityTolerance',1e-8,'StepTolerance',1e-10);
[wOpt,~,flag,solver]=fmincon(@objective,wStart,[],[],[],[],zeros(M,1),100*ones(M,1),[],options);
[~,~,stats0]=objective(wOpt);loScale=0;hiScale=1;target=1+1e-5;
while min(objectiveStats(exp(hiScale)*wOpt).expectations./coeff)<target,hiScale=2*hiScale;assert(hiScale<64);end
if min(stats0.expectations./coeff)>=target,hiScale=0;else
 for it=1:40,mid=(loScale+hiScale)/2;s=objectiveStats(exp(mid)*wOpt);if min(s.expectations./coeff)<target,loScale=mid;else,hiScale=mid;end;end
end
wFeasible=exp(hiScale)*wOpt;st=objectiveStats(wFeasible);
[aa,ix]=sort(u(:,1));bb=cummin(u(ix,2));weights=exp(-bb);union=sum(weights(1:end-1).*exp(-aa(1:end-1)).*(-expm1(-diff(aa))))+exp(-aa(end)-bb(end));
result=struct('matlab',version,'status','finite positive quadrature dual/primal only; continuum average and cost not interval-certified','n',n,'nodes',numel(t),'spacing',spacing,'span_factor',spanFactor,'gamma',gamma,'receivers',size(u,1),'groups',M,'height_cap',A,'c',c,'tau_ref',tauRef,'solver_exitflag',flag,'solver_iterations',solver.iterations,'solver_funcCount',solver.funcCount,'pre_correction_min_average',min(stats0.expectations./coeff),'log_multiplier_correction',hiScale,'min_average',min(st.expectations./coeff),'primal_logJ',st.logJ,'dual_scaled',st.dual,'relative_duality_gap',(st.Jscaled-st.dual)/st.Jscaled,'union_volume',union,'union_log_score',log(union)-st.logJ,'same_c_scalar_log_score',baseline.log_score,'union_minus_scalar',log(union)-st.logJ-baseline.log_score,'weights',wFeasible,'v_nonnegative',v,'receiver_d',d,'averages_by_group',st.expectations./coeff,'probability_mass_error',max(abs(masses-coeff)),'seconds',toc);
fid=fopen(out,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fprintf('n%d N%d M%d spacing%g span%g gamma%g flag%d minAvg%.8g gap%.3g unionV%.8g score%.10f vsScalar%+.7g sec%.2f\n',n,numel(t),size(u,1),spacing,spanFactor,gamma,flag,min(st.expectations./coeff),result.relative_duality_gap,union,result.union_log_score,result.union_minus_scalar,toc);
function [f,g,st]=objective(w)
 arg=R*w;dd=log(arg);h=profile(dd,n,A);expect=P'*h;active=h>0;yy=log(h(active));logPhi=yy+(n-1)*log1p(max(yy,0));Jscaled=sum(exp(logRef(active)+logPhi));dual=coeff'*w-w'*expect+Jscaled;f=-dual;g=expect-coeff;
 if nargout>2,st=struct('expectations',expect,'Jscaled',Jscaled,'dual',dual,'logJ',-tauRef+log(Jscaled));end
end
function st=objectiveStats(w)
 [~,~,st]=objective(w);
end
end
function h=profile(dd,n,A)
h=zeros(size(dd));h(dd>=0)=1;soft=dd>log(n);ys=max(expm1((dd(soft)-log(n))/(n-1)),0);
for j=1:7,res=(n-2)*log1p(ys)+log(n+ys)-dd(soft);der=(n-2)./(1+ys)+1./(n+ys);ys=max(0,ys-res./der);end
h(soft)=exp(min(ys,log(A)));
end
function lr=logR(b,c,t)
k=0:b;lw=gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1);xx=t-c*k;lg=lw+(b-1)*log(max(xx,realmin));lg(xx<=0)=-Inf;pk=max(lg,[],2);lr=pk+log(sum(exp(lg-pk),2))-gammaln(b)-b*log1p(exp(-c));lr(t<=0)=-Inf;
end
function [t,w]=quad(N,lo,hi)
ii=(1:15)';bb=ii./sqrt(4*ii.^2-1);[V,D]=eig(diag(bb,1)+diag(bb,-1));[x,ix]=sort(diag(D));ww=2*V(1,ix)'.^2;panels=ceil(N/16);edges=linspace(lo,hi,panels+1);half=diff(edges)/2;mid=(edges(1:end-1)+edges(2:end))/2;t=x*half+mid;w=ww*half;t=t(:);w=w(:);
end
