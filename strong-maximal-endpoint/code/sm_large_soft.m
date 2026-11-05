function rows=sm_large_soft
% Large-dimensional calibration with central binomial truncation, not a certificate.
root=fileparts(fileparts(mfilename('fullpath')));
B=fzero(@(b)log(b)-1-1/b,[3 4]); mu=log(B); sigma=sqrt(mu); r=1/B;
rows=cell(0,1);
for n=[1e4 1e5 1e6 1e7 1e8]
    theta=1/n; baseP=r/(1+r);
    logMGF=n*(-log1p(-theta)+log1p(baseP*expm1(-theta*log(r))));
    p=r^(1-theta)/(1+r^(1-theta)); mode=floor((n+1)*p);
    sd=sqrt(n*p*(1-p)); lower=max(0,floor(mode-12*sd)); upper=min(n,ceil(mode+12*sd));
    k=lower:upper; logWeights=zeros(size(k)); mid=mode-lower+1;
    ku=mode:upper-1;
    logWeights(mid+1:end)=cumsum(log(((n-ku)*p)./((ku+1)*(1-p))));
    kd=mode:-1:lower+1;
    logWeights(mid-1:-1:1)=cumsum(log((kd*(1-p))./((n-kd+1)*p)));
    weights=exp(logWeights); weights=weights/sum(weights);
    omitted=binocdf(lower-1,n,p)+binocdf(upper,n,p,'upper');
    lo=0; hi=n*mu;
    for j=1:70
        q=(lo+hi)/2;
        failure=sum(weights.*gammainc(max((1-theta)*(q+k*log(r)),0),n,'lower'));
        logE=logMGF-theta*q+log1p(-failure);
        if logE>1e-10, lo=q; else, hi=q; end
    end
    q=lo;
    K=min(n-1,ceil(12*sqrt(n)));
    terms=[0 cumsum(log1p(-(0:K-1)/(n-1)))];
    integralSum=sum(exp(terms));
    ratio=max(0,1-K/(n-1));
    remainder=exp(terms(end))*ratio/(1-ratio);
    logI=log(n/(n-1))+log(integralSum+remainder);
    score=q-logI; thirdScale=sqrt(n/log(n));
    residual=(score-n*mu+sigma*sqrt(n*log(n)))/thirdScale-sigma*log(log(n));
    row=struct('n',n,'q',q,'binomial_terms',numel(k), ...
        'binomial_omitted_probability_numeric',omitted, ...
        'integral_tail_upper_estimate_numeric',remainder, ...
        'log_lower_candidate',score,'log_integral_factor',logI, ...
        'third_constant_residual',residual, ...
        'predicted_third_constant',sigma*log(sigma*sqrt(2*pi)));
    rows{end+1}=row; %#ok<AGROW>
    fprintf('n=%.0e terms=%d logL=%.12g thirdConstant=%.9g target=%.9g omitted=%.2g\n', ...
        n,numel(k),score,residual,row.predicted_third_constant,omitted);
end
fid=fopen(fullfile(root,'rounds','r008','large_soft.json'),'w');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(struct('matlab_version',version,'B',B,'mu',mu,'r',r,'status', ...
    'double-precision truncated-mixture calibration; not directed certification','rows',{rows}),PrettyPrint=true));
end
