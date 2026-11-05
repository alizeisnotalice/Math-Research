function rows=sm_exponential
% Infinite soft profile, with exact Gamma/binomial exponential tilting.
% C_n >= exp(q)/I_n(a) if E exp(a(T-q)/n)1_(T>=q)>1.
root=fileparts(fileparts(mfilename('fullpath')));
B=fzero(@(b)log(b)-1-1/b,[3 4]); mu=log(B); r=1/B;
rows=cell(0,1);
for n=[4 8 16 32 64 128 256 512 1024 2048 4096]
    best=-Inf;
    for a=[0.25 0.5 0.75 0.9 1]
        theta=a/n;
        p=r^(1-theta)/(1+r^(1-theta));
        k=0:n;
        lw=gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1) ...
            +k*log(p)+(n-k)*log1p(-p);
        w=exp(lw); w=w/sum(w);
        logMGF=n*(-log1p(-theta)+log1p(r^(1-theta))-log1p(r));
        lo=0; hi=n*mu+10*sqrt(n);
        for iter=1:70
            q=(lo+hi)/2;
            failure=sum(w.*gammainc(max((1-theta)*(q+k*log(r)),0),n,'lower'));
            logE=logMGF-theta*q+log1p(-failure);
            if logE>1e-10, lo=q; else, hi=q; end
        end
        q=lo;
        j=0:n-1;
        terms=gammaln(n)-gammaln(n-j)+j*log(a/n)-(j+1)*log1p(-a/n);
        peak=max(terms); logI=peak+log(sum(exp(terms-peak)));
        score=q-logI;
        if score>best
            best=score;
            row=struct('n',n,'r',r,'a',a,'q',q,'log_integral_factor',logI, ...
                'log_lower_candidate',score,'expectation_log_target',1e-10, ...
                'scaled_gap',(n*mu-score)/sqrt(n*log(n)), ...
                'gap_after_leading',(score-n*mu+sqrt(mu*n*log(n)))/sqrt(n/log(n)));
        end
    end
    rows{end+1}=row; %#ok<AGROW>
    fprintf('n=%4d a=%.2g q=%.9g logL=%.12g scaledGap=%.8g thirdScale=%.8g\n', ...
        n,row.a,row.q,row.log_lower_candidate,row.scaled_gap,row.gap_after_leading);
end
fid=fopen(fullfile(root,'rounds','r003','exponential.json'),'w');
assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(struct('matlab_version',version,'rows',{rows}),PrettyPrint=true));
end
