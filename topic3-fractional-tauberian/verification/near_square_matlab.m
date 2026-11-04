% Regenerate the simplex coefficients through the general intersection formula.
base=fullfile(fileparts(mfilename('fullpath')),'..','results');
source=jsondecode(fileread(fullfile(base,'near_square_python.json')));
checked=0;maxError=0;
for row=1:numel(source.certificate_rows)
 info=source.certificate_rows(row);order=info.order+1;V=zeros(5,4);
 for k=0:4
  V(k+1,order)=[-ones(1,k),ones(1,4-k)];
 end
 assert(ismember(string(info.x_endpoint),["0","1/2"]));endpoint=0;
 if string(info.x_endpoint)=="1/2",endpoint=.5;end
 index=0;
 for i=1:5
  for j=i:5
   if i==j
    coefficient=qvalue(V(i,:),endpoint);
   else
    coefficient=2*qvalue((V(i,:)+V(j,:))/2,endpoint) ...
     -(qvalue(V(i,:),endpoint)+qvalue(V(j,:),endpoint))/2;
   end
   index=index+1;expected=str2double(info.bernstein_upper_triangle{index});
   maxError=max(maxError,abs(coefficient-expected));
   assert(coefficient>=-1e-12 && abs(coefficient-expected)<1e-12);
   checked=checked+1;
  end
 end
end
rng(20261105,'twister');geometrySamples=10000;maxRankResidual=0;maxExpansionResidual=0;
for j=1:geometrySamples
 x=rand;r=1+.01*(2*rand(1,4)-1);u=2*rand(1,4)-1;t=.001;
 [H,R,W]=geometry(x,r,u);assert(H<=2*sum(u)-sum(r)+1e-12);
 exact=2+sum((1+t*W(1:4,1)).*(1+t*W(1:4,2))) ...
  -2*(1+t*W(5,1))*(1+t*W(5,2));
 maxExpansionResidual=max(maxExpansionResidual,abs(exact-(4+t*H+t^2*R)));
 for axis=1:2
  coord=x;if axis==2,coord=1-x;end
  signs=1-2*bitget(0:3,axis);right=signs.*r*coord+u;left=signs.*r*coord-u;
  right=sort(right);left=sort(left);
  rankValue=sum(right(1:2))-sum(left(3:4));
  identity=sum(W(1:4,axis))-2*W(5,axis);
  maxRankResidual=max(maxRankResidual,abs(rankValue-identity));
 end
end
assert(maxRankResidual<1e-12 && maxExpansionResidual<1e-12);
out=struct('matlab',version,'seed',20261105,'bernstein_coefficients',checked, ...
 'coefficient_max_error',maxError,'general_geometry_samples',geometrySamples, ...
 'max_rank_identity_residual',maxRankResidual,'max_expansion_residual',maxExpansionResidual, ...
 'scope','independent double regeneration; exact certificate and analytic transfer establish the theorem');
fid=fopen(fullfile(base,'near_square_matlab.json'),'w');fprintf(fid,'%s\n',jsonencode(out));fclose(fid);disp(out);
function q=qvalue(v,endpoint)
 if endpoint==0
  q=2*loss(.25,v)-loss(.5,v);
 else
  q=loss(.5,v);
 end
end
function q=loss(x,v)
 [~,R]=geometry(x,ones(1,4),x*v);q=(8*x-2*x*sum(v)-R)/x;
end
function [H,R,W]=geometry(x,r,u)
 sx=[1,-1,1,-1];sy=[1,1,-1,-1];W=zeros(5,2);
 for omit=1:5
  members=1:4;if omit<5,members(omit)=[];end
  W(omit,1)=min(sx(members).*r(members)*x+u(members)) ...
   -max(sx(members).*r(members)*x-u(members));
  W(omit,2)=min(sy(members).*r(members)*(1-x)+u(members)) ...
   -max(sy(members).*r(members)*(1-x)-u(members));
 end
 H=sum(W(1:4,:),'all')-2*sum(W(5,:));
 R=sum(W(1:4,1).*W(1:4,2))-2*prod(W(5,:));
end
