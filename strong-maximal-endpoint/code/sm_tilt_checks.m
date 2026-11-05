function result=sm_tilt_checks
% Conditional entropy-cost calibration; analytic uniform bound is in round 12.
root=fileparts(fileparts(mfilename('fullpath')));
B=fzero(@(b)log(b)-1-1/b,[3 4]); mu=log(B); rstar=1/B;
rows=cell(0,1);
for delta=[0.05 0.02 0.01 0.005 0.002 0.001]
    target=mu-delta;
    low=fzero(@(r)meanOriginal(r)-target,[0 rstar]);
    high=fzero(@(r)meanOriginal(r)-target,[rstar 1]);
    grid=linspace(low,high,501);
    costs=arrayfun(@(r)costAt(r,target,delta),grid);
    [~,j]=max(costs); left=grid(max(1,j-2)); right=grid(min(501,j+2));
    [r,negative]=fminbnd(@(r)-costAt(r,target,delta),left,right,optimset('TolX',1e-12));
    [cost,lambda]=costAt(r,target,delta);
    c=-log(r); variance=@(t)1./(1+t).^2+c^2*exp(-c*(1+t))./(1+exp(-c*(1+t))).^2;
    integralCost=integral(@(t)t.*variance(t),0,lambda,'AbsTol',1e-20,'RelTol',1e-12);
    assert(lambda<=2*delta && abs(cost-integralCost)<1e-13);
    row=struct('delta',delta,'sampled_best_r',r,'cost_over_delta_squared',cost/delta^2, ...
        'lambda_over_delta',lambda/delta,'integral_cost_difference',abs(cost-integralCost), ...
        'asymptotic_cost_coefficient',1/(2*mu),'grid_size',501,'solver_negative_cost',negative);
    rows{end+1}=row; %#ok<AGROW>
    fprintf('tilt delta=%.4g cost/delta2=%.12g limit=%.12g lambda/delta=%.9g\n', ...
        delta,row.cost_over_delta_squared,row.asymptotic_cost_coefficient,row.lambda_over_delta);
end
result=struct('matlab_version',version,'status','deterministic scalar calibration, not interval certification','rows',{rows});
fid=fopen(fullfile(root,'rounds','r012','tilt_checks.json'),'w');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end

function m=meanOriginal(r)
if r==0, m=1; else, m=1-r*log(r)/(1+r); end
end

function [cost,lambda]=costAt(r,target,delta)
if meanOriginal(r)<=target, cost=0; lambda=0; return; end
c=-log(r); p=r/(1+r);
meanTilt=@(t)1/(1+t)+c*exp(-c*(1+t))/(1+exp(-c*(1+t)));
lambda=fzero(@(t)meanTilt(t)-target,[0 2*delta],optimset('TolX',1e-14));
cost=-lambda*target+log1p(lambda)-log1p(p*expm1(lambda*log(r)));
end
