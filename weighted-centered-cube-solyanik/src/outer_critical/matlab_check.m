% Independent numerical finite-atom halo and A1 comparison controls.
tauValues=[0.25,0.5,0.75,0.9,1]; cases=0;
for tau=tauValues
    for n=[2,4,6,8]
        x=linspace(-1,1,n); m=ones(1,n)*4/(3*(n-1));
        H=atomicHalo(x,m);
        haloMass=sum((sign(H(:,2)).*abs(H(:,2)).^tau ...
            -sign(H(:,1)).*abs(H(:,1)).^tau)/tau);
        value=haloMass/sum(m.*abs(x).^(tau-1));
        L=1+m(1)/2;
        expected=(2*L^tau/tau)/sum(m.*abs(x).^(tau-1));
        assert(abs(value-expected)<1e-11);
        if n==2
            assert(abs(value-3/(4*tau)*(5/3)^tau)<1e-11);
        end
        cases=cases+1;
    end
end
fprintf('MATLAB %s: %d atomic chain controls PASS.\n',version,cases);
for tau=tauValues(1:4)
    u=fzero(@(v) v^(1-tau)+(1-tau)*v-tau,[0,1]);
    assert(abs(u^(1-tau)+(1-tau)*u-tau)<1e-10);
    upper=2*u^(tau-1);
    fprintf('tau=%.2f pair lower=%.15g weighted comparison upper=%.15g\n', ...
        tau,3/(4*tau)*(5/3)^tau,upper);
    if tau==0.5
        assert(abs(upper-2*(1+sqrt(2)))<1e-10);
    end
end

function H=atomicHalo(x,m)
% Enumerate subsets, independently of the searcher's contiguous-block code.
intervals=[];
for mask=1:2^numel(x)-1
    keep=logical(bitget(mask,1:numel(x)));
    mass=sum(m(keep));a=max(x(keep))-mass/2;b=min(x(keep))+mass/2;
    if a<b,intervals(end+1,:)=[a,b];end %#ok<AGROW>
end
intervals=sortrows(intervals,1);H=intervals(1,:);
for j=2:size(intervals,1)
    if intervals(j,1)<=H(end,2)+1e-14
        H(end,2)=max(H(end,2),intervals(j,2));
    else
        H(end+1,:)=intervals(j,:); %#ok<AGROW>
    end
end
end
