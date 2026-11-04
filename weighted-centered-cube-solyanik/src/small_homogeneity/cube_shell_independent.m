function cube_shell_independent()
% Independent MATLAB adaptive quadrature of the given cube shell formula.
% This is numerical evidence only; no workspace/path reset is performed.
taus=[.1,.01,.0001,1e-6]; ns=[1,2,3]; bs=[1,2,5]; rows=[];
for n=ns
 for b=bs
  for t=taus
   r=1+t; breaks=[0,r-1,r,r+1];
   dA=@(s) localdA(s,n,b,r);
   % The first shell has an integrable s^(t-1) endpoint singularity; integrate
   % it explicitly.  On the other shells use s=exp(u), avoiding endpoint loss.
   mu0=(1+b)*2^(n-1)*n/t*(r-1)^t;
   logintegrand=@(u) exp((t-n+1)*u).*dA(exp(u));
   mu1=integral(logintegrand,log(breaks(2)),log(breaks(3)), ...
       'RelTol',1e-11,'AbsTol',1e-11,'ArrayValued',true);
   mu2=integral(logintegrand,log(breaks(3)),log(breaks(4)), ...
       'RelTol',1e-11,'AbsTol',1e-11,'ArrayValued',true);
   mu=mu0+mu1+mu2;
   R=mu/(2^n*r^n); target=n*(1+b)/2; ea=(1+b)/2;
   ballLower=n*ea/t*t^t/(1+t)^n;
   eta=t^(1/(n+1)); triUpper=2^n*b*3^(n-t)+n*2^n*ea/(n-t)*(1+1/eta)^(n-t)+n*ea/t*eta^t*(1+eta)^(n-t);
   rows=[rows; n,b,t,mu,R,t*R,target,t*R-target,ballLower,triUpper]; %#ok<AGROW>
  end
 end
end
payload=struct('version','round27-cube-shell-matlab-v1-20261105','matlab_version',version, ...
 'columns',{{'n','b','tau','shell_mass','R','tauR','target','error','ball_lower','three_scale_upper'}},'rows',rows);
root=fileparts(mfilename('fullpath'));fid=fopen(fullfile(root,'cube_shell_matlab_summary.json'),'w');fwrite(fid,jsonencode(payload),'char');fclose(fid);
assert(all(rows(:,4)>0) && all(rows(:,5)>0) && all(rows(:,5)>=rows(:,9)) && all(rows(:,5)<=rows(:,10)), ...
    'Shell comparison with lower/upper control failed');
fprintf('cube shell MATLAB: %d rows, all positive\n',size(rows,1));
end
function y=localdA(s,n,b,r)
y=zeros(size(s));
i=s<=r-1; y(i)=(1+b)*2^(n-1)*n*s(i).^(n-1);
i=s>r-1 & s<=r; y(i)=2^(n-1)*(n*s(i).^(n-1)+b*(r-1)*(n-1)*s(i).^(n-2));
i=s>r & s<=r+1; y(i)=(2*r)^(n-1);
end
