function result=sm_root_checks(outputPath)
% Gaussian root only: no finite-n certification of C_n.
if nargin==0, outputPath=fullfile(pwd,'root_checks.json'); end
B=fzero(@(b)log(b)-1-1/b,[3 4]); mu=log(B); sigma=sqrt(mu);
logns=[log([2 4 8 16 32 64 128 256]) [4 8 16 32 64 128 256]*log(10)];
rows=cell(1,numel(logns));
for idx=1:numel(logns)
 L=logns(idx); fun=@(z)log(erfcx(z/sqrt(2))/2)-z.^2/2-log(sigma)-log(z)+L/2;
 lo=realmin; hi=sqrt(L)+2;
 for j=1:70
  mid=(lo+hi)/2; if fun(mid)>0, lo=mid; else, hi=mid; end
 end
 z=(lo+hi)/2; H=log(L)+log(sigma*sqrt(2*pi));
 z3=sqrt(L)-H/sqrt(L); z5=z3+(2*H-1-H^2/2)/L^(3/2);
 rows{idx}=struct('log_n',L,'z',z,'log_equation_residual',fun(z), ...
  'three_term_z_error',z-z3,'next_term_z_error',z-z5);
 fprintf('logn=%10.5g z=%.15g residual=%+.3g next_error=%+.6g\n',L,z,fun(z),z-z5);
end
result=struct('matlab_version',version,'precision','double', ...
 'method','70 deterministic bisections; log(erfcx) tail', ...
 'mu',mu,'sigma',sigma,'rows',{rows});
fid=fopen(outputPath,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
