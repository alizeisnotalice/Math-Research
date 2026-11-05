function result=sm_controlled_gaussian_checks
% Compare stored Bellman diagnostics with the limiting Gaussian tail.
root=fileparts(fileparts(mfilename('fullpath')));
coarse=jsondecode(fileread(fullfile(root,'rounds','r013','bellman_1024.json')));
fine=jsondecode(fileread(fullfile(root,'rounds','r011','bellman_00125.json')));
rows=cell(0,1);
for idx=1:numel(coarse.rows)
    row=coarse.rows(idx);
    if row.n<=256
        match=find([fine.rows.n]==row.n,1); row=fine.rows(match);
    end
    normal=erfc(sqrt(log(row.n)/2))/2;
    out=struct('n',row.n,'state_step',row.step,'gaussian_tail',normal, ...
        'bellman_failure',row.adaptive_failure, ...
        'bellman_over_gaussian',row.adaptive_failure/normal, ...
        'fixed_over_gaussian',row.fixed_failure/normal, ...
        'status','state-discretized diagnostic; no finite-n theorem certificate');
    rows{end+1}=out; %#ok<AGROW>
    fprintf('n=%4d adaptive/Gaussian=%.9g fixed/Gaussian=%.9g\n', ...
        out.n,out.bellman_over_gaussian,out.fixed_over_gaussian);
end
result=struct('matlab_version',version,'rows',{rows});
fid=fopen(fullfile(root,'rounds','r013','gaussian_checks.json'),'w');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
