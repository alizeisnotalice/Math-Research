function rows=sm_dense
% Dense monotone profiles: numerical grid calibration, not continuum optimality.
root=fileparts(fileparts(mfilename('fullpath')));
B=fzero(@(b)log(b)-1-1/b,[3 4]); mu=log(B);
rows=cell(0,1);
for n=[16 64 256 1024]
    for resolution=[0.15 0.075]
        sd=sqrt(mu*n); step=resolution*sqrt(n/log(n));
        grid=max(0,n*mu-6*sd):step:n*mu+8*sd;
        candidate=sm_staircase(n,1/B,grid);
        row=struct('n',n,'resolution',resolution,'cells',numel(grid), ...
            'log_lower_candidate',candidate.log_lower_candidate, ...
            'scaled_gap',(n*mu-candidate.log_lower_candidate)/sqrt(n*log(n)), ...
            'expectation',candidate.expectation,'kkt_residual',candidate.kkt_residual);
        rows{end+1}=row; %#ok<AGROW>
        fprintf('dense n=%d cells=%d logL=%.12g scaledGap=%.8g\n', ...
            n,row.cells,row.log_lower_candidate,row.scaled_gap);
    end
end
fid=fopen(fullfile(root,'rounds','r004','dense.json'),'w');
assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(struct('matlab_version',version,'rows',{rows}),PrettyPrint=true));
end
