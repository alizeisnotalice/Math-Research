% Independent MATLAB numeric root/inversion diagnostic for round 24.
% No clearvars or path reset: safe for the connected user MATLAB session.
r24_root=fileparts(mfilename('fullpath')); r24_ds=unique([linspace(1e-5,1.999,91),2/3-2.^(-[8 16 32]),2/3,2/3+2.^(-[8 16 32])]);
r24_rows=struct([]);r24_maxRootResidual=0;r24_maxInverseError=0;r24_count=0;r24_fail=0;
for r24_d=r24_ds
    if r24_d<2/3
        lo=1/3;hi=1/2;
        for k=1:180
            q=(lo+hi)/2; dq=3*q*(1-q)*(1-2*q)/(1-3*q+3*q^2);
            if dq>r24_d,lo=q;else,hi=q;end
        end
        q=(lo+hi)/2;
        H=(r24_d^2/4)*(3*q*(1-q)+(r24_d/2)*(2*q-1))/(q^2*(1-q)^2);
        residual=6*q^3-3*(3+r24_d)*q^2+3*(1+r24_d)*q-r24_d;
        % Independently invert D=H through the parametric H(q) bisection.
        lo2=1/3;hi2=1/2;
        for k=1:180
            qi=(lo2+hi2)/2;ri=1-3*qi+3*qi^2;Hi=27*qi*(1-qi)*(1-2*qi)^2*(1-2*qi+2*qi^2)/(8*ri^3);
            if Hi>H,lo2=qi;else,hi2=qi;end
        end
        qi=(lo2+hi2)/2;dback=3*qi*(1-qi)*(1-2*qi)/(1-3*qi+3*qi^2);
    else
        q=r24_d/2;H=r24_d*(4-r24_d)/(2-r24_d)^2;residual=0;dback=2-2/sqrt(H+1);
    end
    invErr=abs(dback-r24_d);valid=abs(residual)<=2e-11 && invErr<=2e-11;
    r24_count=r24_count+1;r24_fail=r24_fail+(~valid);r24_maxRootResidual=max(r24_maxRootResidual,abs(residual));r24_maxInverseError=max(r24_maxInverseError,invErr);
    r24_rows(r24_count).d=r24_d; %#ok<SAGROW>
    r24_rows(r24_count).q=q;r24_rows(r24_count).H=H;r24_rows(r24_count).rootResidual=residual;r24_rows(r24_count).inverseError=invErr;r24_rows(r24_count).valid=valid;
end
r24_summary=struct('version','matlab-envelope-independent-v1-20261105','rows',r24_count,'failures',r24_fail,'max_root_residual',r24_maxRootResidual,'max_inverse_error',r24_maxInverseError,'matlab_version',version,'records',{r24_rows});
r24_fid=fopen(fullfile(r24_root,'matlab_envelope_summary.json'),'w');assert(r24_fid>=0,'Cannot write MATLAB envelope summary.');fprintf(r24_fid,'%s\n',jsonencode(r24_summary));fclose(r24_fid);
assert(r24_fail==0,'MATLAB envelope diagnostic failed.');
fprintf('MATLAB envelope: %d rows, failures=%d, root residual %.3g, inverse error %.3g\n',r24_count,r24_fail,r24_maxRootResidual,r24_maxInverseError);
