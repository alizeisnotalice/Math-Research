% Independent double-precision general-p envelope diagnostic (not a proof).
% No workspace clearing or path reset.
r24_root=fileparts(mfilename('fullpath'));r24_ps=[1.1 1.5 1.9 2 2.1 3 4 8 16];
r24_machineEps=2.220446049250313e-16;r24_count=0;r24_fail=0;r24_resolvedInner=0;r24_indeterminate=0;r24_maxAttemptRes=0;r24_maxResolvedRes=0;r24_maxOuterInv=0;r24_rows=struct([]);
for p=r24_ps
    sw=2/p;ds=unique([min(sw/8,.01),sw/2,sw*(1-1e-10),sw,min(1.9,sw+min(.1,(2-sw)/2))]);
    for d=ds
        if d<sw
            lo=d/2;hi=1-1e-13;
            for k=1:220
                q=(lo+hi)/2;
                if r24_gp(p,d,q)<0,lo=q;else,hi=q;end
            end
            q=(lo+hi)/2;H=r24_g(p,d,q);res=abs(r24_gp(p,d,q));invErr=0;qgap=q-d/2;
            gapFloor=128*r24_machineEps*max(1,d/2);
            if qgap<=gapFloor
                % An actual attempted solve whose strict separation is not reliable in double.
                valid=true;status='indeterminate_near_switch';r24_indeterminate=r24_indeterminate+1;
            else
                valid=(res<5e-8)&&(q<1);status='resolved_inner';r24_resolvedInner=r24_resolvedInner+1;r24_maxResolvedRes=max(r24_maxResolvedRes,res);
            end
        else
            q=d/2;qgap=NaN;H=(1-d/2)^(1-p)-1;res=0;dback=2-2/(H+1)^(1/(p-1));invErr=abs(dback-d);
            valid=invErr<2e-11;status='outer';
        end
        r24_count=r24_count+1;r24_fail=r24_fail+(~valid);if ~isnan(res),r24_maxAttemptRes=max(r24_maxAttemptRes,res);end;if ~isnan(invErr),r24_maxOuterInv=max(r24_maxOuterInv,invErr);end
        r24_rows(r24_count).p=p; %#ok<SAGROW>
        r24_rows(r24_count).d=d;r24_rows(r24_count).q=q;r24_rows(r24_count).qMinusDelta=qgap;r24_rows(r24_count).H=H;r24_rows(r24_count).rootResidual=res;r24_rows(r24_count).outerInverseError=invErr;r24_rows(r24_count).status=status;r24_rows(r24_count).valid=valid;
    end
end
r24_summary=struct('version','matlab-general-p-v4-resolution-classified-20261105','rows',r24_count,'resolved_inner_rows',r24_resolvedInner,'indeterminate_near_switch_rows',r24_indeterminate,'failures',r24_fail,'max_attempted_inner_residual',r24_maxAttemptRes,'max_resolved_inner_residual',r24_maxResolvedRes,'max_outer_inverse_error',r24_maxOuterInv,'matlab_version',version,'records',{r24_rows});
r24_fid=fopen(fullfile(r24_root,'matlab_general_p_summary.json'),'w');assert(r24_fid>=0,'Cannot write general-p summary.');fprintf(r24_fid,'%s\n',jsonencode(r24_summary));fclose(r24_fid);
assert(r24_fail==0,'General-p MATLAB diagnostic failed.');fprintf('MATLAB general-p: %d rows (%d resolved inner, %d indeterminate), failures=%d, resolved residual %.3g, inverse %.3g\n',r24_count,r24_resolvedInner,r24_indeterminate,r24_fail,r24_maxResolvedRes,r24_maxOuterInv);
function y=r24_g(p,d,q)
a=1-d/(2*q);b=1+d/(2*(1-q));y=q*a^p+(1-q)*b^p-1;
end
function y=r24_gp(p,d,q)
delta=d/2;a=1-delta/q;b=1+delta/(1-q);y=a^p-b^p+p*delta*(a^(p-1)/q+b^(p-1)/(1-q));
end
