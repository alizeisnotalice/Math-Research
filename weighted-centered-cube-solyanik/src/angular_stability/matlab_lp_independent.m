% Independent finite-vector floating diagnostic for the stated Lp bounds.
% No maximal/halo quantity appears here.  No clearvars or path reset is used.
r22_root=fileparts(mfilename('fullpath'));
r22_p=[1.1 1.25 1.5 1.9 2 2.1 3 4 8 16];
r22_den=[10 4 2 10 1 10 1 1 1 1];
r22_total=0; r22_failure=0; r22_maxViolation=-Inf; r22_records=struct([]);
for r22_seed=1:20
    for r22_ip=1:numel(r22_p)
        p=r22_p(r22_ip); den=r22_den(r22_ip);
        r22_family=mod(r22_seed-1,4);
        if r22_family==0
            delta=2^(-16); r=[.5 .5]; f=[1-delta 1+delta]; label='balanced';
        elseif r22_family==1
            q=2^(-9); r=[q 1-q]; f=[q 1+q]; label='hole';
        elseif r22_family==2
            M=2^(den*10); eps=M^(-(p+1)); c=(1-eps*M)/(1-eps);
            r=[eps 1-eps]; f=[M c]; label='spike';
        else
            rw=mod(r22_seed*7+(1:4)*3,11)+1; r=rw/sum(rw);
            f0=(mod(r22_seed*13+(1:4)*5,29)+1)./(mod(r22_seed+(1:4),7)+2);
            f=f0/(r*f0.'); label='random';
        end
        u=sum(r.*r22_bregman(f,p));
        dp=sum(r.*abs(f-1).^p); d2=sum(r.*(f-1).^2);
        if p<=2
            rhs=expm1((2/p)*log1p(u))/(p-1); lhs=dp^(2/p); violation=lhs-rhs;
            valid=violation<=2e-11*max([1 abs(lhs) abs(rhs)]);
        else
            lpv=dp-u; varv=d2-u/(p-1); violation=max(lpv,varv);
            valid=violation<=2e-11*max([1 abs(u) abs(d2)]);
        end
        r22_total=r22_total+1; r22_failure=r22_failure+(~valid); r22_maxViolation=max(r22_maxViolation,violation);
        r22_records(r22_total).seed=r22_seed; %#ok<SAGROW>
        r22_records(r22_total).p=p;
        r22_records(r22_total).family=label;
        r22_records(r22_total).bregman_u=u;
        r22_records(r22_total).dp=dp;
        r22_records(r22_total).d2=d2;
        r22_records(r22_total).maxViolation=violation;
        r22_records(r22_total).valid=valid;
    end
end
r22_summary=struct('version','matlab-lp-angular-v1-20261105','examples',r22_total, ...
    'failures',r22_failure,'max_violation',r22_maxViolation,'matlab_version',version,'records',{r22_records});
r22_fid=fopen(fullfile(r22_root,'matlab_lp_summary.json'),'w');
assert(r22_fid>=0,'Cannot open MATLAB Lp summary.');
fprintf(r22_fid,'%s\n',jsonencode(r22_summary)); fclose(r22_fid);
assert(r22_failure==0,'A finite-vector diagnostic failed.');
fprintf('MATLAB Lp independent: %d examples, failures=%d, max violation %.3g\n',r22_total,r22_failure,r22_maxViolation);

function y=r22_bregman(f,p)
d=f-1;
y=expm1(p*log1p(d))-p*d;
end
