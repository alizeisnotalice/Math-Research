function result = sm_staircase(n,r,thresholds,margin)
% Global convex amplitude solve on a fixed threshold grid, via scalar dual.
% Probabilities are numerical, so the output is not an interval certificate.
if nargin<4, margin=1e-9; end
t=thresholds(:)';
assert(n>=2 && all(diff(t)>0) && all(t>=0));
anchor=t(1);
[alpha,~,failure]=arrayfun(@(q) sm_tail(n,r,q),t);
p=alpha-[alpha(2:end) 0];
% Lower-tail differences avoid cancellation when both upper tails round to 1.
useLower=find(failure(2:end)<0.5);
p(useLower)=failure(useLower+1)-failure(useLower);
w=exp(anchor-t)-[exp(anchor-t(2:end)) 0];
assert(all(p>0) && all(w>0),'Grid probabilities underflowed');
logRatio=log(p)-log(w);
assert(all(diff(logRatio)>-1e-10),'Likelihood ratio ordering failed');
target=1+margin;
lo=-max(logRatio)-50; hi=-min(logRatio)+50+log(n);
for k=1:100
    mid=(lo+hi)/2;
    h=amplitudes(n,mid+logRatio);
    if p*h'<target, lo=mid; else, hi=mid; end
end
logLambda=hi;
h=amplitudes(n,logLambda+logRatio);
% A discrete cell may lie at the derivative=1 jump: permit a fractional level.
hLo=amplitudes(n,lo+logRatio);
changed=find(h-hLo>0.5 & hLo<1e-8 & h<=1+1e-8);
if ~isempty(changed)
    j=changed(1); h(j)=0;
    h(j)=(target-p*h')/p(j);
end
expectation=p*h';
assert(abs(expectation-target)<2e-8,'Dual solve expectation residual');
assert(all(diff(h)>-1e-8),'Nonmonotone amplitudes');
phi=h.*(1+max(log(max(h,realmin)),0)).^(n-1);
cost=w*phi';
logScore=anchor-log(cost);
der=ones(size(h));
above=h>1;
der(above)=(1+log(h(above))).^(n-2).*(n+log(h(above)));
dualSlope=exp(logLambda+logRatio);
errors=zeros(size(h));
errors(h==0)=max(dualSlope(h==0)-1,0);
errors(h==1)=max(max(1-dualSlope(h==1),dualSlope(h==1)-n),0);
free=h>0 & h~=1;
errors(free)=abs(der(free)-dualSlope(free))./max(1,dualSlope(free));
result=struct('n',n,'r',r,'thresholds',t,'heights',h, ...
    'expectation',expectation,'margin',margin,'log_lower_candidate',logScore, ...
    'log_lambda_scaled',logLambda,'kkt_residual',max(errors), ...
    'precision','IEEE double; not directed bounds');
end

function h=amplitudes(n,z)
h=ones(size(z));
h(z<0)=0;
idx=z>log(n);
v=max((z(idx)-log(n))/(n-1),0);
for j=1:30
    f=(n-2)*log1p(v)+log(n+v)-z(idx);
    d=(n-2)./(1+v)+1./(n+v);
    v=max(v-f./d,0);
end
h(idx)=exp(v);
end
