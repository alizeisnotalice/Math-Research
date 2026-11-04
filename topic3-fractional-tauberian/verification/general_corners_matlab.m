% Independent adaptive integration of all exact rational rectangle regressions.
base=fullfile(fileparts(mfilename('fullpath')),'..','results');
source=jsondecode(fileread(fullfile(base,'general_corners_python.json')));
maxError=0; checked=0;
for j=1:numel(source.independent_rectangle_regressions)
 row=source.independent_rectangle_regressions(j);k=zeros(4,2);
 for c=1:4
  for d=1:2
   k(c,d)=ratval(row.corners{c}{d});
  end
 end
 result=0;
 for c=0:3
  result=result+integral(@(x)cornerValue(x,k,c),0,1, ...
   'Waypoints',linspace(0,1,33),'AbsTol',1e-9,'RelTol',1e-9);
 end
 err=abs(result-ratval(row.B));assert(err<2e-7);
 maxError=max(maxError,err);checked=checked+1;
end
% Independent random-density defect check, including nonmonotone pixel holes.
rng(20261105,'twister');maxResidual=0;densityChecks=10000;
for j=1:densityChecks
 x=rand;m=min(x,1-x);u=(2*rand-1)*m;s=x-u;q=1-x-u;
 f=rand(4);M=sum(f,'all')/4;Tx=0;Ty=0;Txy=0;lower=0;
 for a=1:4
  for b=1:4
   loX=(a-1)/2;hiX=a/2;loY=(b-1)/2;hiY=b/2;
   dx=max(hiX-max(loX,s),0);dy=max(hiY-max(loY,q),0);
   Tx=Tx+f(a,b)*dx/2;Ty=Ty+f(a,b)*dy/2;Txy=Txy+f(a,b)*dx*dy;
   lower=lower+f(a,b)*max(min(hiX,s)-loX,0)*max(min(hiY,q)-loY,0);
  end
 end
 phi=u^2+u+x*(1-x)+M-Tx-Ty-Txy;
 defect=2*(m^2-u^2)+s*q-lower+2*Txy;
 maxResidual=max(maxResidual,abs(2*m-phi-defect));
 assert(defect>=-1e-12 && abs(2*m-phi-defect)<1e-12);
end
out=struct('matlab',version,'seed',20261105,'exact_rectangle_cases',checked, ...
 'adaptive_integral_max_error',maxError,'density_identity_samples',densityChecks, ...
 'max_identity_residual',maxResidual, ...
 'scope','independent double regressions; universal inequalities and remainder are analytic');
fid=fopen(fullfile(base,'general_corners_matlab.json'),'w');fprintf(fid,'%s\n',jsonencode(out));fclose(fid);disp(out);
function v=ratval(s)
 parts=split(string(s),'/');v=str2double(parts(1));if numel(parts)==2,v=v/str2double(parts(2));end
end
function val=cornerValue(x,k,c)
 m=min(x,1-x);ix=bitxor(c,1)+1;iy=bitxor(c,2)+1;ij=bitxor(c,3)+1;
 U=[-m;m;x-k(ix,1);1-x-k(iy,2);x-k(ij,1);1-x-k(ij,2)];
 valid=U>=-m-1e-14 & U<=m+1e-14;U=min(max(U,-m),m);s=x-U;q=1-x-U;
 mass=k(:,1).*k(:,2);
 G=U.^2+U+x.*(1-x)+sum(mass)/2-mass(c+1) ...
  -max(k(ix,1)-s,0)*k(ix,2)-k(iy,1)*max(k(iy,2)-q,0) ...
  -max(k(ij,1)-s,0).*max(k(ij,2)-q,0);
 G(~valid)=-inf;val=max(G,[],1);
end
