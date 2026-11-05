function result = sm_nearflat_ceiling_checks(outputPath)
% Deterministic diagnostics for the proved near-flat full-cube ceiling.
% No optimization/random sampling. Double arithmetic, not interval certificates.
if nargin < 1
    outputPath = fullfile(fileparts(mfilename('fullpath')), '..', 'rounds', 'r024', 'checks.json');
end
ns = [2 4 8 16 32 64 128 256 1e4 1e8 1e16 1e64 1e256];
rows = zeros(numel(ns), 6);
for j = 1:numel(ns)
    n = ns(j); ep = 1/n; delta = 1/(n*sqrt(log(n)));
    logR = log1p(-delta);
    logKappa = log1p(-delta*(1-log(delta)));
    shell = -expm1(n*logR);
    logPhi = log1p(ep)+(n-1)*log1p(log1p(ep));
    denom = 1+shell*expm1(logPhi);
    lower = 1; upper = 1;
    for k = 1:min(n,40)
        logChoose = 0;
        for ii = 0:k-1
            logChoose = logChoose + log(n-ii)-log(ii+1);
        end
        meanG = -expm1(k*logR+(n-k)*logKappa);
        ellLow = log1p((1-ep^2)*ep*meanG);
        ellUp = log1p(ep);
        lower = lower + exp(logChoose+k*log(2)+logH(k,ellLow));
        upper = upper + exp(logChoose+k*log(2)+logH(k,ellUp));
    end
    rows(j,:) = [n shell denom lower/denom upper upper/exp(logPhi)];
end
mu = fzero(@(x) x-1-exp(-x),[1 2]);
baseNs = [2 4 8 16 32 64 128 256];
capRows = zeros(3*numel(baseNs),6); jj=0;
for n = baseNs
    epsilons = [1/n 1/sqrt(n*log(n)) sqrt(log(n)/n)];
    for scale = 1:3
        jj=jj+1; ep=epsilons(scale); A=1+ep;
        q=n*mu-sqrt(2*mu*n*log(2*A/ep));
        logPhi=log1p(ep)+(n-1)*log1p(log1p(ep));
        % I0(2 sqrt(x)) <= exp(2 sqrt(x)): analytic halo upper bound.
        logCubeBound=log(A)+2*sqrt(2*n*log1p(ep));
        capRows(jj,:)=[n scale ep q q-logPhi logCubeBound];
    end
end
result=struct('matlab',version,'arithmetic','double; no random seed; no optimizer', ...
    'columns',{{'n','shell_mass','denominator','Jensen_witness_ratio','cube_halo_partial','constant_endpoint_ratio'}}, ...
    'rows',rows,'series_cutoff',40,'omitted_halo_tail_bound',2*(2^41/exp(2*gammaln(42)))/(1-2/42^2), ...
    'S',besseli(0,2*sqrt(2)),'S_over_e',besseli(0,2*sqrt(2))/exp(1), ...
    'cap_columns',{{'n','scale_id','epsilon','q','hyperbolic_log_lower','cube_log_upper'}},'cap_rows',capRows);
fid=fopen(outputPath,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
disp(rows); disp([result.S result.S_over_e]);
end
function value=logH(k,ell)
% H_k(t)=ell^k/k! sum_{j>=0} ell^j/j! * k/(k+j).
term=1; series=1;
for j=1:20
    term=term*ell/j;
    series=series+term*k/(k+j);
end
value=k*log(ell)-gammaln(k+1)+log(series);
end
