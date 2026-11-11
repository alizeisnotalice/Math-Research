% Read-only formula sanity check executed for work/r086/audit86.md.
% This file records the exact successful MATLAB command; do not rerun it
% as part of this archival step.
format long g; mu=fzero(@(u) u-1-exp(-u),1.3); r=exp(-mu); p=r/(1+r); m=@(x) 1-x.*log(x)./(1+x); fprintf('mu=%.15g r=%.15g p=%.15g H=%.15g m(r)=%.15g\n',mu,r,p,-p*log(p)-(1-p)*log(1-p),m(r)); rr=[1e-8 0.05 r 0.7 1]; fprintf('max m grid=%.15g, m(1)=%.15g\n',max(m(rr)),m(1)); b=.2; tailSecond=4/(1-b)^3; fprintf('tail moment bound 4/(1-b)^3=%.15g\n',tailSecond);
