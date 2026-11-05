function rows = sm_baseline(outputPath)
% Reproducible equal-r search: numerical candidates, not a proof of optimal H_n.
if nargin == 0
    outputPath = fullfile(fileparts(fileparts(mfilename('fullpath'))),'rounds','r001','baseline.json');
end
ns = [2 4 8 16 32 64 128 256];
B = fzero(@(b) log(b)-1-1/b,[3 4]);
mu = log(B); sigma = sqrt(mu);
opts = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
    'OptimalityTolerance',1e-9,'StepTolerance',1e-11, ...
    'MaxFunctionEvaluations',1500,'MaxIterations',400);
rows = struct([]);
for n = ns
    limit = n*mu+20*sqrt(n)+50;
    q0 = max(0.2,n*mu-sigma*sqrt(n*log(n)));
    bestValue = -Inf; bestX = []; bestFlag = NaN;
    for r0 = [0.08 1/B 0.65]
        for qStart = [q0 max(0.5,0.7*n*mu)]
            [x,value,flag] = fmincon(@(x) -score(n,x),[r0 qStart], ...
                [],[],[],[],[1e-7 1e-7],[1 limit],[],opts);
            if -value > bestValue
                bestValue = -value; bestX = x; bestFlag = flag;
            end
        end
    end
    % Explicit endpoint families and q=0 are additional controls.
    for r = [0 1]
        [q,value] = fminbnd(@(q) -score(n,[r q]),0,limit, ...
            optimset('TolX',1e-10,'Display','off'));
        if -value > bestValue
            bestValue = -value; bestX = [r q]; bestFlag = 1;
        end
    end
    [alpha,logAlpha,failure] = sm_tail(n,bestX(1),bestX(2));
    qLo = 0; qHi = limit;
    for j = 1:80
        mid = (qLo+qHi)/2;
        [~,logMid] = sm_tail(n,bestX(1),mid);
        if logMid > logAlpha, qLo = mid; else, qHi = mid; end
    end
    qRecovered = (qLo+qHi)/2;
    qGrid = linspace(0,limit,101);
    tailGrid = arrayfun(@(q) sm_tail(n,bestX(1),q),qGrid);
    assert(all(diff(tailGrid) <= 2e-14),'Tail monotonicity failed');
    asymptotic = n*mu-sigma*sqrt(n*log(n)) ...
        +sigma*(log(sigma*sqrt(2*pi))-1)*sqrt(n/log(n));
    row = struct('n',n,'r',bestX(1),'q',bestX(2),'alpha',alpha, ...
        'failure',failure,'log_lower_candidate',bestValue, ...
        'scaled_gap',(n*mu-bestValue)/sqrt(n*log(n)), ...
        'three_term_prediction',asymptotic,'exitflag',bestFlag, ...
        'quantile_residual',abs(qRecovered-bestX(2)), ...
        'q_search_cap',limit,'r_search_floor',1e-7, ...
        'small_threshold_control_log',n*log(2)-gammaln(n));
    if isempty(rows)
        rows = row;
    else
        rows(end+1) = row; %#ok<AGROW>
    end
    fprintf('n=%3d r=%.10g q=%.10g alpha=%.10g logL=%.12g gap=%.9g\n', ...
        n,row.r,row.q,row.alpha,bestValue,row.scaled_gap);
end
metadata = struct('matlab_version',version,'precision','IEEE double', ...
    'method','Gamma/binomial sum; six SQP starts plus r=0,1 controls', ...
    'status','numerical fixed-r witness candidates; no global optimality claim', ...
    'B',B,'mu',mu,'rows',rows);
fid = fopen(outputPath,'w');
assert(fid>=0,'Cannot open output');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(metadata,PrettyPrint=true));
end

function value = score(n,x)
[~,logAlpha] = sm_tail(n,x(1),x(2));
value = logAlpha+x(2)-(n-1)*log1p(-logAlpha);
if ~isfinite(value), value = -1e100; end
end
