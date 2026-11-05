function rows=sm_box
% Exact single-box halo formula evaluated by deterministic quadrature.
root=fileparts(fileparts(mfilename('fullpath')));
rows=cell(0,1);
for n=[8 32 128 512]
    A=sqrt(log(n)/n); k=1:n;
    integrals=arrayfun(@(j)integral(@(u)exp(A*u).*u.^(j-1),0,1, ...
        'AbsTol',1e-13,'RelTol',1e-11),k);
    logs=gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1) ...
        +k*log(2*A)-gammaln(k)+log(integrals);
    peak=max([0 logs]); logVolume=peak+log(exp(-peak)+sum(exp(logs-peak)));
    logCost=A+(n-1)*log1p(A);
    bound=2*sqrt(2*exp(A)*n*A);
    assert(logVolume<=bound+1e-10);
    row=struct('n',n,'A',A,'log_halo_volume',logVolume, ...
        'volume_log_upper_bound',bound,'log_modulus',logCost, ...
        'log_normalized_ratio',logVolume-logCost);
    rows{end+1}=row; %#ok<AGROW>
    fprintf('box n=%d logV=%.10g logCost=%.10g logRatio=%.10g\n', ...
        n,logVolume,logCost,row.log_normalized_ratio);
end
fid=fopen(fullfile(root,'rounds','r005','box.json'),'w');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(struct('matlab_version',version,'rows',{rows}),PrettyPrint=true));
end
