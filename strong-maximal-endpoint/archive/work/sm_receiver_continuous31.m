function result=sm_receiver_continuous31(seedPath,outPath)
root=fileparts(mfilename('fullpath'));
if nargin<1||isempty(seedPath),seedPath=fullfile(root,'variable29_archived_theta0.json');end
if nargin<2||isempty(outPath),outPath=fullfile(root,'sm_receiver_continuous31.json');end
s=jsondecode(fileread(seedPath));n=s.n;b=n/2;A=s.height_cap;c=s.c;tau=s.tau_ref;w=s.weights(:);v=s.v_nonnegative(:);d=s.receiver_d(:);M=numel(w);coeff=[1;2*ones(M-1,1)];
u=[0 0;[d(2:end)+v(2:end),d(2:end)-v(2:end)];[d(2:end)-v(2:end),d(2:end)+v(2:end)]];
p=1/(1+exp(c));mz=1+c*p;vz=1+c*c*p*(1-p);lo=min(u,[],'all');hi=max(u,[],'all')+b*mz+16*sqrt(b*vz);levels=[1 n phiPrime(A,n)];
initial128=evaluate(1,128,32,levels,n,A,b,c,tau,u,w,coeff,lo,hi);initial256=evaluate(1,256,64,levels,n,A,b,c,tau,u,w,coeff,lo,hi);
target=1+1e-5;low=1;high=1;steps=0;
while initial128.min_average<target && high<4
    low=high;high=1.02*high;initial128=evaluate(high,128,32,levels,n,A,b,c,tau,u,w,coeff,lo,hi);steps=steps+1;
end
assert(high<4,'Scale bracket failed');
for k=1:18
    mid=(low+high)/2;test=evaluate(mid,128,32,levels,n,A,b,c,tau,u,w,coeff,lo,hi);
    if test.min_average<target,low=mid;else,high=mid;end
end
searchScale=high;safeScale=searchScale*1.0005;
final32=evaluate(safeScale,256,32,levels,n,A,b,c,tau,u,w,coeff,lo,hi);final=evaluate(safeScale,256,64,levels,n,A,b,c,tau,u,w,coeff,lo,hi);
tailActual=phi(A,n)*(2*exp(-lo-hi)-exp(-2*hi));tailAvg=tailAverageUpper(A,b,c,hi,u)./coeff;
result=struct('status','fixed-seed continuous quadrature with analytic constant-branch integration; floating numerical evidence only, not interval certified', ...
 'source_file',seedPath,'n',n,'height_cap',A,'c',c,'tau_ref',tau,'weights',w,'receiver_d',d,'finite_domain',[lo hi], ...
 'initial_N128_inner32',initial128,'initial_N256_inner64',initial256,'target_min_average',target,'scale_at_search_threshold',searchScale,'uniform_scale_with_extra_margin',safeScale,'scale_search_steps',steps, ...
 'final_N256_inner32',final32,'final_N256_inner64',final,'finite_domain_tail_actual_cost_upper',tailActual,'finite_domain_tail_scaled_cost_upper',exp(tau)*tailActual,'finite_domain_tail_average_upper_by_group',tailAvg);
fid=fopen(outPath,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fprintf('initial %.12g final scale %.12g minAvg %.12g actualLogJ %.12g\n',initial256.min_average,safeScale,final.min_average,final.actual_logJ);
end

function out=evaluate(scale,N,inner,levels,n,A,b,c,tau,u,w,coeff,lo,hi)
tic;[t,qw]=quad16(N,lo,hi);cost=0;avg=zeros(numel(w),1);roots=zeros(numel(t),3);
for i=1:numel(t)
    r=arrayfun(@(z) boundary(t(i),z,scale,b,c,tau,u,w,lo,hi),levels);roots(i,:)=r;
    [c1,a1]=constantPart(t(i),r(1),r(2),1,n,b,c,tau,u);[cA,aA]=constantPart(t(i),r(3),hi,A,n,b,c,tau,u);
    [cs,as]=softPart(t(i),r(2),r(3),inner,scale,n,A,b,c,tau,u,w);
    cost=cost+qw(i)*(c1+cA+cs);avg=avg+qw(i)*(a1+A*aA+as);
end
out=struct('outer_nodes',numel(t),'soft_inner_order',inner,'scaled_cost',cost,'scaled_log_cost',log(cost),'actual_logJ',log(cost)-tau, ...
 'averages_by_group',avg./coeff,'min_average',min(avg./coeff),'boundary_root_min',min(roots,[],1),'boundary_root_max',max(roots,[],1),'seconds',toc);
end

function r=boundary(t1,level,scale,b,c,tau,u,w,lo,hi)
f=@(x) likelihood(t1,x,scale,b,c,tau,u,w)-level;fl=f(lo);fh=f(hi);
if fl>=0,r=lo;return;end
if fh<=0,r=hi;return;end
left=lo;right=hi;
for k=1:70
    mid=(left+right)/2;
    if f(mid)<0,left=mid;else,right=mid;end
end
r=(left+right)/2;
end

function [cost,avg]=constantPart(t1,left,right,h,n,b,c,tau,u)
M=(size(u,1)+1)/2;avg=zeros(M,1);cost=0;if right<=left,return;end
cost=exp(tau-t1)*(exp(-left)-exp(-right))*phi(h,n);
for j=1:M
    pts=receiverPoints(j,u);
    for k=1:size(pts,1)
        z1=t1-pts(k,1);fac=exp(logPoly(b,c,z1)+pts(k,1)-t1);
        if isfinite(fac),avg(j)=avg(j)+fac*(cdfGB(b,c,right-pts(k,2))-cdfGB(b,c,left-pts(k,2)));end
    end
end
end

function [cost,avg]=softPart(t1,left,right,inner,scale,n,A,b,c,tau,u,w)
M=numel(w);avg=zeros(M,1);cost=0;if right<=left,return;end
[x,q]=quadNodes(inner,left,right);[R,P]=kernels(t1,x,b,c,tau,u);H=profileS(scale*(R*w),n,A);cost=sum(q.*exp(tau-t1-x).*phi(H,n));avg=P'*(q.*H);
end

function z=likelihood(t1,t2,scale,b,c,tau,u,w)
[R,~]=kernels(t1,t2,b,c,tau,u);z=scale*(R*w);
end
function [R,P]=kernels(t1,t2,b,c,tau,u)
t2=t2(:);M=(size(u,1)+1)/2;R=zeros(numel(t2),M);P=zeros(numel(t2),M);
for j=1:M
    pts=receiverPoints(j,u);
    for k=1:size(pts,1)
        lr=logPoly(b,c,t1-pts(k,1))+logPoly(b,c,t2-pts(k,2))+sum(pts(k,:));
        R(:,j)=R(:,j)+exp(lr-tau);P(:,j)=P(:,j)+exp(lr-t1-t2);
    end
end
end
function pts=receiverPoints(j,u)
if j==1,pts=[0 0];else,M=(size(u,1)+1)/2;pts=[u(j,1),u(j,2);u(j+M-1,1),u(j+M-1,2)];end
end
function y=profileS(s,n,A)
y=zeros(size(s));y(s>=1)=1;ix=s>n;d=log(s(ix));z=max(expm1((d-log(n))/(n-1)),0);
for k=1:10
    res=(n-2)*log1p(z)+log(n+z)-d;der=(n-2)./(1+z)+1./(n+z);z=max(0,z-res./der);
end
y(ix)=min(exp(z),A);
end
function z=phi(h,n),z=h.*(1+max(log(max(h,1)),0)).^(n-1);end
function z=phiPrime(A,n),z=(1+log(A))^(n-2)*(n+log(A));end
function y=logPoly(b,c,t)
k=0:b;lw=gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1);xx=t-c*k;lg=lw+(b-1)*log(max(xx,realmin));lg(xx<=0)=-Inf;pk=max(lg,[],2);y=pk+log(sum(exp(lg-pk),2))-gammaln(b)-b*log1p(exp(-c));y(t<=0)=-Inf;
end
function F=cdfGB(b,c,z)
k=0:b;pk=exp(gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1)-c*k-b*log1p(exp(-c)));F=sum(pk.*gammainc(max(z-c*k,0),b,'lower'));
end
function tail=tailAverageUpper(A,b,c,hi,u)
M=(size(u,1)+1)/2;tail=zeros(M,1);
for j=1:M
    pts=receiverPoints(j,u);
    for k=1:size(pts,1)
        t1=tailGB(b,c,hi-pts(k,1));t2=tailGB(b,c,hi-pts(k,2));tail(j)=tail(j)+A*(t1+(1-t1)*t2);
    end
end
end
function T=tailGB(b,c,z)
k=0:b;pk=exp(gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1)-c*k-b*log1p(exp(-c)));T=sum(pk.*gammainc(max(z-c*k,0),b,'upper'));
end
function [t,w]=quad16(N,lo,hi)
panels=ceil(N/16);e=linspace(lo,hi,panels+1);[x,ww]=quadNodes(16,-1,1);half=diff(e)/2;mid=(e(1:end-1)+e(2:end))/2;t=x*half+mid;w=ww*half;t=t(:);w=w(:);
end
function [t,w]=quadNodes(m,lo,hi)
i=(1:m-1)';bb=i./sqrt(4*i.^2-1);[V,D]=eig(diag(bb,1)+diag(bb,-1));[x,ix]=sort(diag(D));ww=2*V(1,ix)'.^2;t=x*(hi-lo)/2+(lo+hi)/2;w=ww*(hi-lo)/2;t=t(:);w=w(:);
end
