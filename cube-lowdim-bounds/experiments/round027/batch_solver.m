% One owned -batch process, one linprog call; never the shared MATLAB GUI.
r27here=fileparts(mfilename('fullpath'));
r27ipc=fullfile(r27here,'ipc.local');
r27control=jsondecode(fileread(fullfile(r27ipc,'control.json')));
assert(r27control.maximum_calls==1 && r27control.maximum_soft_seconds==20);
r27attempted=0; r27completed=0;
atomicJson(fullfile(r27ipc,'worker_started.json'),struct('status','started', ...
    'pid',matlabProcessID,'version',version,'maximum_calls',1, ...
    'basis_warm_start',false,'primitive_formula_calls',0));
try
    for r27iteration=1:1
        r27requestfile=fullfile(r27ipc,sprintf('request_%02d.ready.json',r27iteration));
        while ~isfile(r27requestfile)
            if isfile(fullfile(r27ipc,'stop.json')) || unixNow()>=r27control.deadline_unix
                atomicJson(fullfile(r27ipc,'worker_done.json'),struct('status','stopped_waiting', ...
                    'LP_started',r27attempted,'LP_completed',r27completed));
                return
            end
            pause(0.05);
        end
        r27request=jsondecode(fileread(r27requestfile));
        assert(r27request.iteration==r27iteration);
        r27modelname=sprintf('model_%02d.mat',r27iteration);
        assert(strcmp(r27request.model_file,r27modelname));
        r27modelpath=fullfile(r27ipc,r27modelname);
        assert(strcmp(shaFile(r27modelpath),r27request.input_sha256));
        r27model=load(r27modelpath);
        assert(size(r27model.A_ub,1)==51574 && size(r27model.A_ub,2)==288);
        assert(all(isfinite(nonzeros(r27model.A_ub))) && all(isfinite(r27model.b_ub(:))));
        assert(all(isfinite(r27model.c(:))) && all(isfinite(r27model.lb(:))) && all(isfinite(r27model.ub(:))));
        r27soft=min([20,double(r27request.soft_seconds),r27control.deadline_unix-unixNow()]);
        assert(r27soft>0);
        r27options=optimoptions('linprog','Algorithm','dual-simplex-highs', ...
            'Display','off','MaxTime',r27soft,'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);
        r27record=struct('iteration',r27iteration,'input_sha256',r27request.input_sha256, ...
            'model_file',r27modelname,'status','started','rows',size(r27model.A_ub,1), ...
            'columns',size(r27model.A_ub,2),'soft_seconds',r27soft, ...
            'started_unix',unixNow(),'basis_warm_start',false,'shared_solver_core','HiGHS');
        r27attempted=r27attempted+1;
        atomicJson(fullfile(r27ipc,sprintf('started_%02d.json',r27iteration)),r27record);
        r27clock=tic;
        try
            [r27x,r27f,r27flag,r27output]=linprog(r27model.c(:),r27model.A_ub, ...
                r27model.b_ub(:),[],[],r27model.lb(:),r27model.ub(:),r27options);
            r27record.linprog_seconds=toc(r27clock);
            r27record.status='returned';r27record.exitflag=r27flag;
            r27record.output=r27output;r27record.x=r27x;r27record.objective=r27f;
            r27record.finite_primal=~isempty(r27x) && all(isfinite(r27x));
            if r27record.finite_primal
                r27record.minimum_row_slack=min(r27model.b_ub(:)-r27model.A_ub*r27x);
                r27record.minimum_bound_slack=min([r27x-r27model.lb(:);r27model.ub(:)-r27x]);
                r27record.full_cost=double(r27model.base_cost)+r27model.c(:)'*r27x;
                r27record.existing_matrix_residual_entries=size(r27model.A_ub,1);
            end
        catch r27error
            r27record.status='solver_error';r27record.error_identifier=r27error.identifier;
            r27record.error_message=r27error.message;r27record.finite_primal=false;
        end
        r27record.elapsed_seconds=toc(r27clock);r27record.completed_unix=unixNow();
        r27completed=r27completed+1;
        r27resultname=sprintf('result_%02d.json',r27iteration);
        atomicJson(fullfile(r27ipc,r27resultname),r27record);
        atomicJson(fullfile(r27ipc,sprintf('result_%02d.ready.json',r27iteration)), ...
            struct('iteration',r27iteration,'result_file',r27resultname, ...
            'result_sha256',shaFile(fullfile(r27ipc,r27resultname)), ...
            'input_sha256',r27request.input_sha256));
        if ~strcmp(r27record.status,'returned') || r27record.exitflag~=1 || ~r27record.finite_primal
            break
        end
    end
    atomicJson(fullfile(r27ipc,'worker_done.json'),struct('status','returned', ...
        'LP_started',r27attempted,'LP_completed',r27completed));
catch r27error
    atomicJson(fullfile(r27ipc,'worker_error.json'),struct('status','worker_error', ...
        'identifier',r27error.identifier,'message',r27error.message, ...
        'LP_started',r27attempted,'LP_completed',r27completed));
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
