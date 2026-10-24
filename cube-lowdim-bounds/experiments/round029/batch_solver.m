% One dedicated -batch, three fixed-beta models; never the shared MATLAB GUI.
r29here=fileparts(mfilename('fullpath'));
r29ipc=fullfile(r29here,'ipc.local');
r29control=jsondecode(fileread(fullfile(r29ipc,'control.json')));
assert(r29control.maximum_calls==3 && r29control.maximum_soft_seconds==20);
assert(r29control.rows==51574 && r29control.columns==128);
r29attempted=0; r29completed=0;
atomicJson(fullfile(r29ipc,'worker_started.json'),struct('status','started', ...
    'pid',matlabProcessID,'version',version,'maximum_calls',3, ...
    'basis_warm_start',false,'primitive_formula_calls',0));
try
    for r29iteration=1:3
        r29requestfile=fullfile(r29ipc,sprintf('request_%02d.ready.json',r29iteration));
        while ~isfile(r29requestfile)
            if isfile(fullfile(r29ipc,'stop.json')) || unixNow()>=r29control.deadline_unix
                atomicJson(fullfile(r29ipc,'worker_done.json'),struct('status','stopped_waiting', ...
                    'LP_started',r29attempted,'LP_completed',r29completed));
                return
            end
            pause(0.05);
        end
        if isfile(fullfile(r29ipc,'stop.json')) || unixNow()>=r29control.deadline_unix
            atomicJson(fullfile(r29ipc,'worker_done.json'),struct('status','stopped_before_request', ...
                'LP_started',r29attempted,'LP_completed',r29completed));
            return
        end
        r29request=jsondecode(fileread(r29requestfile));
        assert(r29request.iteration==r29iteration);
        r29modelname=sprintf('model_%02d.mat',r29iteration);
        assert(strcmp(r29request.model_file,r29modelname));
        r29modelpath=fullfile(r29ipc,r29modelname);
        assert(strcmp(shaFile(r29modelpath),r29request.input_sha256));
        r29model=load(r29modelpath);
        assert(size(r29model.A_ub,1)==51574 && size(r29model.A_ub,2)==128);
        assert(numel(r29model.b_ub)==51574 && numel(r29model.c)==128);
        assert(numel(r29model.lb)==128 && numel(r29model.ub)==128);
        assert(all(isfinite(nonzeros(r29model.A_ub))) && all(isfinite(r29model.b_ub(:))));
        assert(all(isfinite(r29model.c(:))) && all(isfinite(r29model.lb(:))) && all(isfinite(r29model.ub(:))));
        assert(all(r29model.lb(:)<=r29model.ub(:)));
        r29soft=min([20,double(r29request.soft_seconds),r29control.deadline_unix-unixNow()]);
        assert(r29soft>0);
        r29options=optimoptions('linprog','Algorithm','dual-simplex-highs', ...
            'Display','off','MaxTime',r29soft,'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);
        r29record=struct('iteration',r29iteration,'input_sha256',r29request.input_sha256, ...
            'model_file',r29modelname,'status','started','rows',size(r29model.A_ub,1), ...
            'columns',size(r29model.A_ub,2),'soft_seconds',r29soft, ...
            'started_unix',unixNow(),'basis_warm_start',false,'shared_solver_core','HiGHS');
        r29attempted=r29attempted+1;
        atomicJson(fullfile(r29ipc,sprintf('started_%02d.json',r29iteration)),r29record);
        r29clock=tic;
        try
            [r29x,r29f,r29flag,r29output]=linprog(r29model.c(:),r29model.A_ub, ...
                r29model.b_ub(:),[],[],r29model.lb(:),r29model.ub(:),r29options);
            r29record.linprog_seconds=toc(r29clock);
            r29record.status='returned';r29record.exitflag=r29flag;
            r29record.output=r29output;r29record.x=r29x;r29record.objective=r29f;
            r29record.finite_primal=~isempty(r29x) && all(isfinite(r29x));
            if r29record.finite_primal
                r29record.minimum_row_slack=min(r29model.b_ub(:)-r29model.A_ub*r29x);
                r29record.minimum_bound_slack=min([r29x-r29model.lb(:);r29model.ub(:)-r29x]);
                r29record.full_cost=double(r29model.base_cost)+r29model.c(:)'*r29x;
                r29record.existing_matrix_residual_entries=size(r29model.A_ub,1);
            end
        catch r29error
            r29record.status='solver_error';r29record.error_identifier=r29error.identifier;
            r29record.error_message=r29error.message;r29record.finite_primal=false;
        end
        r29record.elapsed_seconds=toc(r29clock);r29record.completed_unix=unixNow();
        r29completed=r29completed+1;
        r29resultname=sprintf('result_%02d.json',r29iteration);
        atomicJson(fullfile(r29ipc,r29resultname),r29record);
        atomicJson(fullfile(r29ipc,sprintf('result_%02d.ready.json',r29iteration)), ...
            struct('iteration',r29iteration,'result_file',r29resultname, ...
            'result_sha256',shaFile(fullfile(r29ipc,r29resultname)), ...
            'input_sha256',r29request.input_sha256));
        if ~strcmp(r29record.status,'returned')
            break
        end
        % A returned nonoptimal status rejects this beta, without a retry.
        % The next beta is a different predeclared model, not a warm start.
    end
    atomicJson(fullfile(r29ipc,'worker_done.json'),struct('status','returned', ...
        'LP_started',r29attempted,'LP_completed',r29completed));
catch r29error
    atomicJson(fullfile(r29ipc,'worker_error.json'),struct('status','worker_error', ...
        'identifier',r29error.identifier,'message',r29error.message, ...
        'LP_started',r29attempted,'LP_completed',r29completed));
end

function t=unixNow()
    t=posixtime(datetime('now','TimeZone','UTC'));
end

function atomicJson(path,value)
    assert(~isfile(path),'Immutable IPC artifact already exists.');
    temporary=[path '.tmp'];
    fid=fopen(temporary,'w');assert(fid>=0);
    cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',jsonencode(value));clear cleanup;
    [ok,message]=movefile(temporary,path);assert(ok,message);
end

function value=shaFile(path)
    fid=fopen(path,'rb');assert(fid>=0);
    cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
    md=java.security.MessageDigest.getInstance('SHA-256');
    md.update(typecast(bytes(:),'int8'));
    raw=typecast(md.digest(),'uint8');
    value=lower(reshape(dec2hex(raw,2).',1,[]));
end
