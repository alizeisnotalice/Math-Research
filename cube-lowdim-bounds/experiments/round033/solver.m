% One fixed discovery LP, one dedicated MATLAB batch. No GUI changes.
r33here=fileparts(mfilename('fullpath'));
r33=load(fullfile(r33here,'model.local.mat'));
assert(size(r33.A,1)==67870 && size(r33.A,2)==26);
assert(all(isfinite(nonzeros(r33.A))) && all(isfinite(r33.b)));
assert(all(isfinite(r33.f)) && all(r33.lb<=r33.ub));
r33out=struct('status','started','calls_started',0,'calls_completed',0,'version',version,'shared_solver_core','HiGHS');
r33out.input_sha256=shaFile(fullfile(r33here,'model.local.mat'));
writeJson(fullfile(r33here,'solver_result.json'),r33out);
r33clock=tic;
try
 options=optimoptions('linprog','Algorithm','dual-simplex-highs','Display','off','MaxTime',20,'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);
 r33out.calls_started=1;writeJson(fullfile(r33here,'solver_result.json'),r33out);
 [x,f,flag,output]=linprog(r33.f,r33.A,r33.b,[],[],r33.lb,r33.ub,options);
 r33out.status='returned';r33out.calls_completed=1;r33out.x=x;r33out.objective=f;r33out.exitflag=flag;r33out.output=output;
catch ex
 r33out.status='solver_error';r33out.message=ex.message;
end
r33out.elapsed_seconds=toc(r33clock);
writeJson(fullfile(r33here,'solver_result.json'),r33out);
function writeJson(path,obj)
 fid=fopen([path '.tmp'],'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(obj));fclose(fid);movefile([path '.tmp'],path);
end

function value=shaFile(path)
 fid=fopen(path,'rb');assert(fid>=0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
 md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes(:),'int8'));
 raw=typecast(md.digest(),'uint8');value=lower(reshape(dec2hex(raw,2).',1,[]));
end
