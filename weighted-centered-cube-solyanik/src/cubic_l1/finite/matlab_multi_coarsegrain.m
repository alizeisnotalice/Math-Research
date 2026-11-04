% Independent float coarse-graining diagnostic for finite p=3 vectors.
% It does not clear the user's workspace or reset MATLAB paths.
r23_root=fileparts(mfilename('fullpath')); r23_dims=[3 4 8 16 64];
r23_C=(11+5*sqrt(5))/8; r23_exponents=[8 16 32]; r23_total=0; r23_fail=0; r23_nonzero=0; r23_constant=0; r23_minPositiveRatio=Inf;
r23_maxJensenViolation=-Inf; r23_records=struct([]);
for r23_n=r23_dims
    for r23_seed=1:30
        r23_kind=mod(r23_seed-1,5);
        if r23_kind==0
            q=0.381966+2e-6*mod(r23_seed,3); a=2^(-r23_exponents(mod(r23_seed-1,3)+1));
            [r,f]=r23_embedded(r23_n,q,a); label='embedded';
        elseif r23_kind==1
            q=0.381966+2e-6*mod(r23_seed,3); a=2^(-r23_exponents(mod(r23_seed-1,3)+1));
            [r,f]=r23_spread(r23_n,q,a); label='spread';
        elseif r23_kind==2
            rw=mod(r23_seed*11+(1:r23_n)*7,97)+1; r=rw/sum(rw);
            x=(mod(r23_seed*17+(1:r23_n)*13,89)+1)./(mod(r23_seed+(1:r23_n)*3,23)+2);
            f=x/(r*x.'); label='random';
        elseif r23_kind==3
            r=ones(1,r23_n)/r23_n; f=ones(1,r23_n); f(1)=.5; f(2)=1.5; label='equal_atom';
        else
            r=ones(1,r23_n)/r23_n; f=ones(1,r23_n); label='constant';
        end
        T=sum(r.*f.^3)-1; d=sum(r.*abs(f-1));
        if d==0
            Tbin=0; dbin=0; ratio=0; jv=0; valid=abs(T)<=5e-12;
        else
            low=f<1; high=~low; q=sum(r(low)); al=sum(r(low).*f(low))/q; b=sum(r(high).*f(high))/(1-q);
            Tbin=q*al^3+(1-q)*b^3-1; dbin=q*(1-al)+(1-q)*(b-1);
            ratio=T/d^2; jv=Tbin-T; valid=(T+5e-10*max(1,abs(T))>=r23_C*d^2) && (jv<=5e-10) && (abs(d-dbin)<=5e-10);
        end
        r23_total=r23_total+1; r23_fail=r23_fail+(~valid);
        if d>0, r23_nonzero=r23_nonzero+1; r23_minPositiveRatio=min(r23_minPositiveRatio,ratio); else, r23_constant=r23_constant+1; end
        r23_maxJensenViolation=max(r23_maxJensenViolation,jv);
        r23_records(r23_total).n=r23_n; %#ok<SAGROW>
        r23_records(r23_total).seed=r23_seed; r23_records(r23_total).kind=label;
        r23_records(r23_total).ratio=ratio; r23_records(r23_total).jensenGap=T-Tbin;
        r23_records(r23_total).dDifference=d-dbin; r23_records(r23_total).valid=valid;
    end
end
r23_summary=struct('version','matlab-multi-p3-v2-boundary-20261105','examples',r23_total,'nonzero_examples',r23_nonzero,'constant_examples',r23_constant,'failures',r23_fail, ...
    'min_positive_distance_ratio',r23_minPositiveRatio,'max_jensen_violation',r23_maxJensenViolation,'matlab_version',version,'records',{r23_records});
r23_fid=fopen(fullfile(r23_root,'matlab_multi_summary.json'),'w');assert(r23_fid>=0,'Cannot write MATLAB multi summary.');
fprintf(r23_fid,'%s\n',jsonencode(r23_summary));fclose(r23_fid);assert(r23_fail==0,'Finite coarse-graining diagnostic failed.');
fprintf('MATLAB multi coarsegrain: %d rows (%d nonzero, %d constant), failures=%d, min positive ratio %.15g\n',r23_total,r23_nonzero,r23_constant,r23_fail,r23_minPositiveRatio);

function [r,f]=r23_embedded(n,q,a)
lo=max(1,floor(n/2));hi=n-lo;r=[repmat(q/lo,1,lo),repmat((1-q)/hi,1,hi)];b=(1-q*a)/(1-q);f=[repmat(a,1,lo),repmat(b,1,hi)];
end
function [r,f]=r23_spread(n,q,a)
lo=max(1,floor(n/2));hi=n-lo;r=[repmat(q/lo,1,lo),repmat((1-q)/hi,1,hi)];b=(1-q*a)/(1-q);
x=a*ones(1,lo);for i=1:lo-1,x(i)=a*(1+((-1)^i)/8);end;x(end)=lo*a-sum(x(1:end-1));
y=b*ones(1,hi);for i=1:hi-1,y(i)=b*(1+((-1)^i)/9);end;y(end)=hi*b-sum(y(1:end-1));f=[x y];
end
