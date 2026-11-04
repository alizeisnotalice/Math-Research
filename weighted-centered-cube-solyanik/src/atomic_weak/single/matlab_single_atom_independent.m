% Independent single-atom halo/integral and optimization diagnostic.
% No clearvars or path reset.
r25_root=fileparts(mfilename('fullpath'));r25_taus=[1/512 1/256 1/128 1/64 1/32 1/16 1/12 1/10 1/8 1/7 1/6 1/5 1/4 .3 1/3];
r25_rows=struct([]);r25_count=0;r25_fail=0;r25_maxGeomErr=0;r25_maxCrit=0;
for t=r25_taus
    lo=1+1e-12;hi=4;
    assert(r25_critical(t,lo)>0 && r25_critical(t,hi)<0);
    for k=1:220
        mid=(lo+hi)/2;
        if r25_critical(t,mid)>0,lo=mid;else,hi=mid;end
    end
    r=(lo+hi)/2;J=r25_J(t,r);crit=r25_critical(t,r);
    mass=(r-1)^t/t+(r+1)^t/t;direct=mass/(2*r);geomErr=abs(direct-J);
    valid=(r>1)&&(abs(crit)<2e-8)&&(geomErr<2e-12);r25_count=r25_count+1;r25_fail=r25_fail+(~valid);r25_maxGeomErr=max(r25_maxGeomErr,geomErr);r25_maxCrit=max(r25_maxCrit,abs(crit));
    r25_rows(r25_count).tau=t; %#ok<SAGROW>
    r25_rows(r25_count).r=r;r25_rows(r25_count).J=J;r25_rows(r25_count).criticalResidual=crit;r25_rows(r25_count).geometryError=geomErr;r25_rows(r25_count).valid=valid;
end
t=1/8;r=17/16;Jfixed=r25_J(t,r);B=6*(5/3)^(1/8);assert(Jfixed>B);
r25_summary=struct('version','matlab-single-atom-v1-20261105','rows',r25_count,'failures',r25_fail,'max_geometry_error',r25_maxGeomErr,'max_critical_residual',r25_maxCrit,'fixed_J',Jfixed,'old_B',B,'matlab_version',version,'records',{r25_rows});
r25_fid=fopen(fullfile(r25_root,'matlab_single_atom_summary.json'),'w');assert(r25_fid>=0,'Cannot write MATLAB single-atom summary.');fprintf(r25_fid,'%s\n',jsonencode(r25_summary));fclose(r25_fid);assert(r25_fail==0,'MATLAB single-atom diagnostic failed.');
fprintf('MATLAB single atom: %d rows, failures=%d, geometry %.3g, critical %.3g\n',r25_count,r25_fail,r25_maxGeomErr,r25_maxCrit);
function y=r25_J(t,r)
y=((r-1)^t+(r+1)^t)/(2*t*r);
end
function y=r25_critical(t,r)
N=(r-1)^t+(r+1)^t;Np=t*((r-1)^(t-1)+(r+1)^(t-1));y=r*Np-N;
end
