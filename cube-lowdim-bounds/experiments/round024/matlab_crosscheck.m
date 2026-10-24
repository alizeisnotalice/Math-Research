% Round 24: one finite LP cross-check; no new kernel evaluations.
% MATLAB and Python both use HiGHS; this is not an independent solver core.
r24here = fileparts(mfilename('fullpath'));
r24out = fullfile(r24here,'matlab_result.json');
assert(~isfile(r24out),'Existing result: no automatic repeat.');
r24m = load(fullfile(r24here,'last_model.mat'));
assert(size(r24m.A_ub,1)<=60000 && size(r24m.A_ub,2)<=288);
assert(all(isfinite(nonzeros(r24m.A_ub))) && all(isfinite(r24m.b_ub)));
r24options = optimoptions('linprog','Algorithm','dual-simplex-highs', ...
    'Display','off','MaxTime',30,'ConstraintTolerance',1e-9, ...
    'OptimalityTolerance',1e-9);
r24clock = tic;
r24record = struct('status','started','rows',size(r24m.A_ub,1), ...
    'columns',size(r24m.A_ub,2),'max_time_soft_seconds',30, ...
    'shared_solver_core','HiGHS','kernel_evaluations',0);
r24fid=fopen(r24out,'w'); fprintf(r24fid,'%s\n',jsonencode(r24record)); fclose(r24fid);
try
    [r24x,r24f,r24flag,r24output] = linprog(r24m.c(:),r24m.A_ub, ...
        r24m.b_ub(:),[],[],r24m.lb(:),r24m.ub(:),r24options);
    r24record.status='returned';
    r24record.exitflag=r24flag;
    r24record.output=r24output;
    r24record.x=r24x;
    r24record.objective=r24f;
    r24record.residual_entries=0;
    if ~isempty(r24x) && all(isfinite(r24x))
        r24residual=r24m.b_ub(:)-r24m.A_ub*r24x;
        r24record.residual_entries=numel(r24residual);
        r24record.minimum_row_slack=min(r24residual);
        r24record.minimum_bound_slack=min([r24x-r24m.lb(:);r24m.ub(:)-r24x]);
        r24record.full_cost=double(r24m.base_cost)+r24m.c(:)'*r24x;
        if ~isempty(r24m.python_objective)
            r24record.objective_difference=r24f-double(r24m.python_objective);
        end
    end
catch r24error
    r24record.status='error';
    r24record.error_identifier=r24error.identifier;
    r24record.error_message=r24error.message;
end
r24record.elapsed_seconds=toc(r24clock);
r24record.matlab_version=version;
r24fid=fopen(r24out,'w'); fprintf(r24fid,'%s\n',jsonencode(r24record)); fclose(r24fid);
disp(r24record);
