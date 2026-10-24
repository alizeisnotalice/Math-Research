% Exactly one finite necessary-point LP, returning the fifth output lambda.
r47here=fileparts(mfilename('fullpath'));
r47=load(fullfile(r47here,'model.local.mat'));
assert(isequal(size(r47.A),[51601,26]));
assert(all(isfinite(nonzeros(r47.A))) && all(isfinite(r47.b)));
assert(all(isfinite(r47.f)) && all(r47.lb<=r47.ub));
r47out=struct('status','started','calls_started',0,'calls_completed',0,'version',version);
r47out.input_sha256=shaFile(fullfile(r47here,'model.local.mat'));
writeJson(fullfile(r47here,'solver_result.json'),r47out);
r47clock=tic;
try
 options=optimoptions('linprog','Algorithm','dual-simplex-highs','Display','off','MaxTime',20,'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);
 r47out.calls_started=1;writeJson(fullfile(r47here,'solver_result.json'),r47out);
 [x,f,flag,output,lambda]=linprog(r47.f,r47.A,r47.b,[],[],r47.lb,r47.ub,options);
 r47out.status='returned';r47out.calls_completed=1;r47out.x=x;r47out.objective=f;
 r47out.exitflag=flag;r47out.output=output;r47out.lambda=lambda;
catch ex
 r47out.status='solver_error';r47out.message=ex.message;
end
r47out.elapsed_seconds=toc(r47clock);
writeJson(fullfile(r47here,'solver_result.json'),r47out);
function writeJson(path,obj)
 fid=fopen([path '.tmp'],'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(obj));fclose(fid);movefile([path '.tmp'],path);
end
function value=shaFile(path)
 fid=fopen(path,'rb');assert(fid>=0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
 md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes(:),'int8'));
 raw=typecast(md.digest(),'uint8');value=lower(reshape(dec2hex(raw,2).',1,[]));
end
