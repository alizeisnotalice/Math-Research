function result=sm_receiver_validation30(seedPath,outPath,Ns,scan,innerOrder)
% Fixed-H floating quadrature validation for archived r029 theta=0 data.
% Roots of likelihood=1,n,Phi'(A) split every t2 integration panel.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(seedPath),seedPath=fullfile(root,'variable29_archived_theta0.json');end
if nargin<2 || isempty(outPath),outPath=fullfile(root,'sm_receiver_validation30.json');end
old=jsondecode(fileread(seedPath));
n=old.n; b=n/2; A=old.height_cap; c=old.c; tau=old.tau_ref; w=old.weights(:);
v=old.v_nonnegative(:); d=old.receiver_d(:); M=numel(v); coeff=[1;2*ones(M-1,1)];
assert(numel(w)==M && numel(d)==M);
u=[0 0;[d(2:end)+v(2:end),d(2:end)-v(2:end)];[d(2:end)-v(2:end),d(2:end)+v(2:end)]];
p=1/(1+exp(c)); mz=1+c*p; vz=1+c*c*p*(1-p);
lo=min(u,[],'all'); hi=max(u,[],'all')+b*mz+16*sqrt(b*vz);
levels=[1,n,phiPrime(A,n)];
if nargin<3 || isempty(Ns),Ns=[128 256 512];end
if nargin<4 || isempty(scan),scan=2049;end
if nargin<5 || isempty(innerOrder),innerOrder=16;end
runs=struct([]);
for q=1:numel(Ns)
    current=runOne(Ns(q),lo,hi,scan,levels,n,A,b,c,tau,u,w,coeff,innerOrder);
    if q==1,runs=current;else,runs(q)=current;end %#ok<AGROW>
    fprintf('N=%d nodes=%d cost=%.12g minAvg=%.12g roots=[%d,%d] sec=%.2f\n', ...
        Ns(q),runs(q).nodes,runs(q).scaled_cost,runs(q).min_average,runs(q).min_roots,runs(q).max_roots,runs(q).seconds);
end
reference=old.averages_by_group(:);
for q=1:numel(runs)
    runs(q).max_average_difference_from_r029=max(abs(runs(q).averages_by_group-reference));
    runs(q).scaled_cost_difference_from_previous=NaN;
    if q>1,runs(q).scaled_cost_difference_from_previous=runs(q).scaled_cost-runs(q-1).scaled_cost;end
end
result=struct('status','fixed archived H; segmented finite-domain Gauss-Legendre validation only, floating point and not an interval certificate', ...
 'source_file',seedPath,'n',n,'height_cap',A,'c',c,'tau_ref',tau, ...
 'weights',w,'receiver_d',d,'finite_domain',[lo hi],'thresholds_likelihood',levels, ...
 'root_scan_points',scan,'inner_quadrature_order',innerOrder,'runs',runs,'archived_r029_min_average',min(reference), ...
 'archived_r029_primal_logJ',old.primal_logJ);
fid=fopen(outPath,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end

function out=runOne(N,lo,hi,scan,levels,n,A,b,c,tau,u,w,coeff,innerOrder)
tic;[t,qw]=quad16(N,lo,hi); totalCost=0; totalAvg=zeros(numel(w),1); rawCost=0; rawAvg=zeros(numel(w),1); allRoots=zeros(numel(t),1);
for i=1:numel(t)
    [rawR,rawP]=kernels(t(i),t,b,c,tau,u); rawH=profileS(rawR*w,n,A); rawY=log(max(rawH,1));rawPhi=rawH.*(1+max(rawY,0)).^(n-1);
    rawMeasure=exp(tau-t(i)-t);rawCost=rawCost+qw(i)*sum(qw.*rawMeasure.*rawPhi);rawAvg=rawAvg+qw(i)*(rawP'*(qw.*rawH));
    roots=[]; grid=linspace(lo,hi,scan); [R,~]=kernels(t(i),grid,b,c,tau,u); s=R*w;
    for ell=levels
        z=s-ell; hit=find(z(1:end-1).*z(2:end)<0);
        hit=[hit; find(z(1:end-1)==0)]; %#ok<AGROW>
        for h=hit(:)'
            r=fzero(@(x) likelihood(t(i),x,b,c,tau,u,w)-ell,[grid(h),grid(h+1)],optimset('TolX',1e-11,'Display','off'));
            roots=[roots;r]; %#ok<AGROW>
        end
    end
    roots=unique(round(roots,11)); roots=roots(roots>lo & roots<hi); allRoots(i)=numel(roots);
    edge=[lo;roots;hi]; innerCost=0; innerAvg=zeros(numel(w),1);
    for h=1:numel(edge)-1
        [x,ww]=quadNodes(innerOrder,edge(h),edge(h+1)); [R,P]=kernels(t(i),x,b,c,tau,u); H=profileS(R*w,n,A);
        yy=log(max(H,1)); Phi=H.*(1+max(yy,0)).^(n-1);
        measure=exp(tau-t(i)-x(:)); innerCost=innerCost+sum(ww(:).*measure.*Phi);
        innerAvg=innerAvg+P'*(ww(:).*H);
    end
    totalCost=totalCost+qw(i)*innerCost; totalAvg=totalAvg+qw(i)*innerAvg;
end
out=struct('requested_N',N,'nodes',numel(t),'inner_quadrature_order',innerOrder, ...
 'scaled_cost',totalCost,'scaled_log_cost',log(totalCost),'actual_logJ',log(totalCost)-tau, ...
 'averages_by_group',totalAvg./coeff,'min_average',min(totalAvg./coeff), ...
 'unsplit_scaled_cost',rawCost,'unsplit_scaled_log_cost',log(rawCost),'unsplit_actual_logJ',log(rawCost)-tau,'unsplit_averages_by_group',rawAvg./coeff,'unsplit_min_average',min(rawAvg./coeff), ...
 'min_roots',min(allRoots),'max_roots',max(allRoots),'mean_roots',mean(allRoots),'seconds',toc);
end

function s=likelihood(t1,t2,b,c,tau,u,w)
[R,~]=kernels(t1,t2,b,c,tau,u);s=R*w;
end
function [R,P]=kernels(t1,t2,b,c,tau,u)
t2=t2(:); M=(size(u,1)+1)/2;R=zeros(numel(t2),M);P=zeros(numel(t2),M);
for j=1:M
    if j==1,points=[0 0];else,points=[u(j,1),u(j,2);u(j+M-1,1),u(j+M-1,2)];end
    for k=1:size(points,1)
        lr=logR(b,c,t1-points(k,1))+logR(b,c,t2-points(k,2))+sum(points(k,:));
        R(:,j)=R(:,j)+exp(lr-tau);P(:,j)=P(:,j)+exp(lr-t1-t2);
    end
end
end
function y=profileS(s,n,A)
y=zeros(size(s));y(s>=1)=1; soft=s>n; d=log(s(soft));z=max(expm1((d-log(n))/(n-1)),0);
for k=1:9,res=(n-2)*log1p(z)+log(n+z)-d;der=(n-2)./(1+z)+1./(n+z);z=max(0,z-res./der);end
y(soft)=min(exp(z),A);
end
function p=phiPrime(A,n),p=(1+log(A))^(n-2)*(n+log(A));end
function y=logR(b,c,t)
k=0:b;lw=gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1);xx=t-c*k;lg=lw+(b-1)*log(max(xx,realmin));lg(xx<=0)=-Inf;pk=max(lg,[],2);y=pk+log(sum(exp(lg-pk),2))-gammaln(b)-b*log1p(exp(-c));y(t<=0)=-Inf;
end
function [t,w]=quad16(N,lo,hi)
panels=ceil(N/16);edges=linspace(lo,hi,panels+1);[x,ww]=quadNodes(16,-1,1);half=diff(edges)/2;mid=(edges(1:end-1)+edges(2:end))/2;t=x*half+mid;w=ww*half;t=t(:);w=w(:);
end
function [t,w]=quadNodes(order,lo,hi)
ii=(1:order-1)';bb=ii./sqrt(4*ii.^2-1);[V,D]=eig(diag(bb,1)+diag(bb,-1));[x,ix]=sort(diag(D));ww=2*V(1,ix)'.^2;t=x*(hi-lo)/2+(lo+hi)/2;w=ww*(hi-lo)/2;t=t(:);w=w(:);
end
