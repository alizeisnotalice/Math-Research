function rows=sm_multilayer
root=fileparts(fileparts(mfilename('fullpath')));
baseline=jsondecode(fileread(fullfile(root,'rounds','r001','baseline.json')));
rows=cell(0,1);
for b=baseline.rows'
    n=b.n;
    if n==2, continue; end % Boundary single-layer control has no finite optimum here.
    best=b.log_lower_candidate;
    winner=[];
    for levels=[3 5]
        for width=[0.25 0.5 1 2 4]
            for offset=[-2 -1 0 1 2]*sqrt(n/log(n))
                t=b.q+offset+linspace(-1,1,levels)*width*sqrt(n/log(n));
                t=unique(max(t,0));
                if numel(t)<2, continue; end
                candidate=sm_staircase(n,b.r,t);
                if candidate.log_lower_candidate>best
                    best=candidate.log_lower_candidate;
                    winner=candidate;
                end
            end
        end
    end
    assert(~isempty(winner),'No improvement found');
    winner.baseline_log=b.log_lower_candidate;
    winner.log_gain=best-b.log_lower_candidate;
    winner.scaled_gap=(n*baseline.mu-best)/sqrt(n*log(n));
    rows{end+1}=winner; %#ok<AGROW>
    fprintf('n=%d steps=%d logL=%.12g gain=%.8g gap=%.8g maxh=%.8g KKT=%.3g\n', ...
        n,numel(winner.heights),best,winner.log_gain,winner.scaled_gap, ...
        max(winner.heights),winner.kkt_residual);
end
out=struct('matlab_version',version,'status','fixed-grid numerical convex amplitude optima; grid search is not global', ...
    'rows',{rows});
fid=fopen(fullfile(root,'rounds','r002','multilayer.json'),'w');
assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(out,PrettyPrint=true));
end
