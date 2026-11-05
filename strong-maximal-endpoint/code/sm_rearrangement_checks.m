function result=sm_rearrangement_checks
% Finite-cell diagnostics of the anchored rearrangement inequality.
root=fileparts(fileparts(mfilename('fullpath')));
old=rng; cleanup=onCleanup(@()rng(old)); rng(20261106,'twister');
B=fzero(@(b)log(b)-1-1/b,[3 4]); m=80; edges=linspace(0,1,m+1);
profiles=2*rand(30,m); sorted=sort(profiles,2,'descend');
maxExcess=-Inf; maxWeightRise=-Inf; maxCostDifference=0; comparisons=0;
for n=[2 4 8]
    phi=@(a)a.*(1+max(log(max(a,realmin)),0)).^(n-1);
    maxCostDifference=max(maxCostDifference,max(abs(sum(phi(profiles),2)-sum(phi(sorted),2)))/m);
    for r=[0 1/B 1]
        for x=[0.25 0.5 1 2 4]
            v=edges/x; F=zeros(size(v)); F(v>=1)=1;
            active=v>0 & v<1;
            F(active)=arrayfun(@(q)sm_tail(n,r,q),-log(v(active)));
            weights=diff(F);
            maxWeightRise=max(maxWeightRise,max(diff(weights)));
            excess=profiles*weights'-sorted*weights';
            maxExcess=max(maxExcess,max(excess)); comparisons=comparisons+numel(excess);
            assert(all(excess<1e-12),'Rearrangement diagnostic failed');
        end
    end
end
result=struct('matlab_version',version,'seed',20261106,'cells',m, ...
    'comparisons',comparisons,'max_pointwise_excess',maxExcess, ...
    'max_weight_increase',maxWeightRise,'max_modulus_difference',maxCostDifference, ...
    'status','finite-cell numerical diagnostics, not proof of rearrangement theorem');
fid=fopen(fullfile(root,'rounds','r009','rearrangement_checks.json'),'w');
assert(fid>=0); fileCleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fprintf('rearrangement comparisons=%d maxExcess=%.3g maxCostDiff=%.3g\n', ...
    comparisons,maxExcess,maxCostDifference);
end
