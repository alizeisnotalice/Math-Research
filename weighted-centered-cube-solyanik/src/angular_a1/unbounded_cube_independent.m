function unbounded_cube_independent()
% Independent double-precision formula evaluation; no workspace/path reset.
ns=[2,3]; gs=[.25,.5,.75,.9]; ts=[.1,.01,.0001,1e-6]; rows=[];
for n=ns
 for g=gs
  ea=1/n+(n-1)/(n-1+g)*(1/(1-g)-1/n);
  for t=ts
   R=localratio(n,g,t); target=(n-g)/(1-g);
   a1=2*3^abs(n-t-g)/(1-g);
   rows=[rows;n,g,t,ea,R,t*R,target,t*R-target,a1]; %#ok<AGROW>
  end
 end
end
payload=struct('version','round28-cube-matlab-v2-20261105','matlab_version',version, ...
 'columns',{{'n','gamma','tau','Ea','atom_R','tau_atom_R','target','error','local_A1_bound'}},'rows',rows);
root=fileparts(mfilename('fullpath'));fid=fopen(fullfile(root,'unbounded_cube_matlab_summary.json'),'w');fwrite(fid,jsonencode(payload),'char');fclose(fid);
assert(all(rows(:,4)>1) && all(rows(:,5)>0) && all(rows(:,9)>1),'Nonpositive angular/atom/A1 diagnostic');
fprintf('unbounded cube MATLAB: %d rows, including q=0 logarithmic branch\n',size(rows,1));
end
function v=localratio(n,g,t)
r=1+t;
v=(localprim(2+t,n,g,t,r)+localprim(t,n,g,t,r))/(2^n*r^n);
end
function v=localprim(L,n,g,t,r)
if L<=r, v=localsmall(L,n,g,t,r); return; end
v=localsmall(r,n,g,t,r)+(2*r)^(n-1)*(L^(t-n+1)-r^(t-n+1))/(t-n+1);
end
function v=localsmall(a,n,g,t,r)
d=n-1;q=t+g-1;z=1-g;
if abs(q)<32*eps
 v=2^d*a^z*((1+d*log(r/a))/z+d/z^2);
else
 v=2^d*((1-d/q)*a^t/t+d*r^q*a^(1-g)/(q*(1-g)));
end
end
