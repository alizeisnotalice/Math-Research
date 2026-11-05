function [alpha, logAlpha, failure] = sm_tail(n, r, q)
% Exact equal-r Gamma/binomial mixture, evaluated in double precision.
% alpha=P(sum Z_r >= q); failure is evaluated separately to retain small tails.
arguments
    n (1,1) double {mustBeInteger,mustBePositive}
    r (1,1) double {mustBeGreaterThanOrEqual(r,0),mustBeLessThanOrEqual(r,1)}
    q (1,1) double {mustBeFinite}
end
if r == 0
    x = max(q,0);
    alpha = gammainc(x,n,'upper');
    failure = gammainc(x,n,'lower');
    logAlpha = log(alpha);
    return
end
k = 0:n;
p = r/(1+r);
logWeights = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1) ...
    + k*log(p)+(n-k)*log1p(-p);
weights = exp(logWeights);
weights = weights/sum(weights);
x = max(q+k*log(r),0);
failure = sum(weights.*gammainc(x,n,'lower'));
upper = gammainc(x,n,'upper');
if failure < 0.5
    logAlpha = log1p(-failure);
    alpha = 1-failure;
else
    terms = logWeights+log(upper);
    peak = max(terms);
    if isinf(peak)
        logAlpha = -Inf;
    else
        logAlpha = peak+log(sum(exp(terms-peak)))-log(sum(exp(logWeights)));
    end
    alpha = exp(logAlpha);
end
end
