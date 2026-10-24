% One owned -batch process, <=6 linprog calls; never the shared MATLAB GUI.
r26here=fileparts(mfilename('fullpath'));
r26ipc=fullfile(r26here,'ipc.local');
r26control=jsondecode(fileread(fullfile(r26ipc,'control.json')));
assert(r26control.maximum_calls==6 && r26control.maximum_soft_seconds==20);
r26attempted=0; r26completed=0;
atomicJson(fullfile(r26ipc,'worker_started.json'),struct('status','started', ...
    'pid',matlabProcessID,'version',version,'maximum_calls',6, ...
    'basis_warm_start',false,'primitive_formula_calls',0));
try
    for r26iteration=1:6
        r26requestfile=fullfile(r26ipc,sprintf('request_%02d.ready.json',r26iteration));
        while ~isfile(r26requestfile)
            if isfile(fullfile(r26ipc,'stop.json')) || unixNow()>=r26control.deadline_unix
                atomicJson(fullfile(r26ipc,'worker_done.json'),struct('status','stopped_waiting', ...
                    'LP_started',r26attempted,'LP_completed',r26completed));
                return
            end
            pause(0.05);
        end
        r26request=jsondecode(fileread(r26requestfile));
        assert(r26request.iteration==r26iteration);
        r26modelname=sprintf('model_%02d.mat',r26iteration);
        assert(strcmp(r26request.model_file,r26modelname));
        r26modelpath=fullfile(r26ipc,r26modelname);
        assert(strcmp(shaFile(r26modelpath),r26request.input_sha256));
        r26model=load(r26modelpath);
        assert(size(r26model.A_ub,1)<=60000 && size(r26model.A_ub,2)<=288);
        assert(all(isfinite(nonzeros(r26model.A_ub))) && all(isfinite(r26model.b_ub(:))));
        assert(all(isfinite(r26model.c(:))) && all(isfinite(r26model.lb(:))) && all(isfinite(r26model.ub(:))));
        r26soft=min([20,double(r26request.soft_seconds),r26control.deadline_unix-unixNow()]);
        assert(r26soft>0);
        r26options=optimoptions('linprog','Algorithm','dual-simplex-highs', ...
            'Display','off','MaxTime',r26soft,'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);
        r26record=struct('iteration',r26iteration,'input_sha256',r26request.input_sha256, ...
            'model_file',r26modelname,'status','started','rows',size(r26model.A_ub,1), ...
            'columns',size(r26model.A_ub,2),'soft_seconds',r26soft, ...
            'started_unix',unixNow(),'basis_warm_start',false,'shared_solver_core','HiGHS');
        r26attempted=r26attempted+1;
        atomicJson(fullfile(r26ipc,sprintf('started_%02d.json',r26iteration)),r26record);
        r26clock=tic;
        try
            [r26x,r26f,r26flag,r26output]=linprog(r26model.c(:),r26model.A_ub, ...
                r26model.b_ub(:),[],[],r26model.lb(:),r26model.ub(:),r26options);
            r26record.linprog_seconds=toc(r26clock);
            r26record.status='returned';r26record.exitflag=r26flag;
            r26record.output=r26output;r26record.x=r26x;r26record.objective=r26f;
            r26record.finite_primal=~isempty(r26x) && all(isfinite(r26x));
            if r26record.finite_primal
                r26record.minimum_row_slack=min(r26model.b_ub(:)-r26model.A_ub*r26x);
                r26record.minimum_bound_slack=min([r26x-r26model.lb(:);r26model.ub(:)-r26x]);
                r26record.full_cost=double(r26model.base_cost)+r26model.c(:)'*r26x;
                r26record.existing_matrix_residual_entries=size(r26model.A_ub,1);
            end
        catch r26error
            r26record.status='solver_error';r26record.error_identifier=r26error.identifier;
            r26record.error_message=r26error.message;r26record.finite_primal=false;
        end
        r26record.elapsed_seconds=toc(r26clock);r26record.completed_unix=unixNow();
        r26completed=r26completed+1;
        r26resultname=sprintf('result_%02d.json',r26iteration);
        atomicJson(fullfile(r26ipc,r26resultname),r26record);
        atomicJson(fullfile(r26ipc,sprintf('result_%02d.ready.json',r26iteration)), ...
            struct('iteration',r26iteration,'result_file',r26resultname, ...
            'result_sha256',shaFile(fullfile(r26ipc,r26resultname)), ...
            'input_sha256',r26request.input_sha256));
        if ~strcmp(r26record.status,'returned') || r26record.exitflag~=1 || ~r26record.finite_primal
            break
        end
    end
    atomicJson(fullfile(r26ipc,'worker_done.json'),struct('status','returned', ...
        'LP_started',r26attempted,'LP_completed',r26completed));
catch r26error
    atomicJson(fullfile(r26ipc,'worker_error.json'),struct('status','worker_error', ...
        'identifier',r26error.identifier,'message',r26error.message, ...
        'LP_started',r26attempted,'LP_completed',r26completed));
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
