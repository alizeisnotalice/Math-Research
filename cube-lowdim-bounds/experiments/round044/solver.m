% One fixed discovery LP, one dedicated MATLAB batch. No GUI changes.
r44here=fileparts(mfilename('fullpath'));
r44=load(fullfile(r44here,'model.local.mat'));
assert(size(r44.A,1)==88319 && size(r44.A,2)==26);
assert(all(isfinite(nonzeros(r44.A))) && all(isfinite(r44.b)));
assert(all(isfinite(r44.f)) && all(r44.lb<=r44.ub));
r44out=struct('status','started','calls_started',0,'calls_completed',0,'version',version,'shared_solver_core','HiGHS');
r44out.input_sha256=shaFile(fullfile(r44here,'model.local.mat'));
writeJson(fullfile(r44here,'solver_result.json'),r44out);
r44clock=tic;
try
 options=optimoptions('linprog','Algorithm','dual-simplex-highs','Display','off','MaxTime',20,'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);
 r44out.calls_started=1;writeJson(fullfile(r44here,'solver_result.json'),r44out);
 [x,f,flag,output]=linprog(r44.f,r44.A,r44.b,[],[],r44.lb,r44.ub,options);
 r44out.status='returned';r44out.calls_completed=1;r44out.x=x;r44out.objective=f;r44out.exitflag=flag;r44out.output=output;
catch ex
 r44out.status='solver_error';r44out.message=ex.message;
end
r44out.elapsed_seconds=toc(r44clock);
writeJson(fullfile(r44here,'solver_result.json'),r44out);
function writeJson(path,obj)
 fid=fopen([path '.tmp'],'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(obj));fclose(fid);movefile([path '.tmp'],path);
end

function value=shaFile(path)
 fid=fopen(path,'rb');assert(fid>=0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
 md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes(:),'int8'));
 raw=typecast(md.digest(),'uint8');value=lower(reshape(dec2hex(raw,2).',1,[]));
end
