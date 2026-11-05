% Independent common-interval budget through general intersection min/max.
base=fullfile(fileparts(mfilename('fullpath')),'..','results');
source=jsondecode(fileread(fullfile(base,'near_rectangle_python.json')));
checked=0;maxCoefficientError=0;
for row=1:numel(source.certificate_rows)
 info=source.certificate_rows(row);order=info.order+1;V=zeros(5,4);
 for k=0:4,V(k+1,order)=[-ones(1,k),ones(1,4-k)];end
 index=0;
 for i=1:5
  for j=i:5
   if i==j
    coefficient=8-Tvalue(V(i,:));
   else
    coefficient=2*(8-Tvalue((V(i,:)+V(j,:))/2)) ...
     -((8-Tvalue(V(i,:)))+(8-Tvalue(V(j,:))))/2;
   end
   index=index+1;expected=str2double(info.bernstein_upper_triangle{index});
   maxCoefficientError=max(maxCoefficientError,abs(coefficient-expected));
   assert(coefficient>=-1e-12 && abs(coefficient-expected)<1e-12);checked=checked+1;
  end
 end
end
rng(20261105,'twister');maxIdentityError=0;maxSharpError=0;samples=10000;
for j=1:samples
 x=rand;p=rand;L=max(p/2-x,-p/2-1+x);U=min(p/2+x,-p/2+1-x);
 h=(U-L)/2;mu=(U+L)/2;v=2*rand(1,4)-1;u=mu+h*v;
 actual=capacityR(x,u)+2*sum(u);
 identity=8*x*(1-x)+h^2*Tvalue(v)+4*mu*h*sum(v)+8*mu^2;
 bound=8*x*(1-x)+8*max(L^2,U^2);
 maxIdentityError=max(maxIdentityError,abs(actual-identity));assert(actual<=bound+1e-12);
 endpoint=L;if U^2>L^2,endpoint=U;end
 sharp=capacityR(x,repmat(endpoint,1,4))+8*endpoint;
 maxSharpError=max(maxSharpError,abs(sharp-bound));
end
maxIntegralError=0;
for p=linspace(0,1,65)
 a=(1-p)/2;b=(1+p)/2;
 exact=integral(@(x)envelope(x,p),0,1,'Waypoints',[a,.5,b],'AbsTol',1e-12,'RelTol',1e-12);
 maxIntegralError=max(maxIntegralError,abs(exact-(2+2*p*(1-p))));
end
x=linspace(0,1,10001);L=max(.25-x,-1.25+x);U=min(.25+x,.75-x);
Uc=[L;U;x-1.25;.5-x];valid=Uc>=L-1e-14 & Uc<=U+1e-14;u=min(max(Uc,L),U);
s=x-u+.25;q=1-x-u-.25;
G=u.^2+u+x.*(1-x)+.5*(x-.5)+1/16+3/8 ...
 -max(1.5-s,0)/4-1.5*max(.25-q,0)-max(1.5-s,0).*max(.25-q,0);
G(~valid)=-inf;G=max(G,[],1);expected=exampleG(x);maxExampleError=max(abs(G-expected));
exampleCubic=4*integral(@exampleG,0,1,'Waypoints',[1/8,7/8],'AbsTol',1e-12);
assert(maxIdentityError<1e-12 && maxSharpError<1e-12 && maxIntegralError<1e-12 ...
 && maxExampleError<1e-12 && abs(exampleCubic-29/16)<1e-12);
out=struct('matlab',version,'seed',20261105,'bernstein_coefficients',checked, ...
 'max_coefficient_error',maxCoefficientError,'common_interval_samples',samples, ...
 'max_shift_identity_error',maxIdentityError,'max_sharp_envelope_error',maxSharpError, ...
 'envelope_integral_parameters',65,'max_integral_error',maxIntegralError, ...
 'example_samples',numel(x),'max_example_error',maxExampleError,'example_cubic',exampleCubic, ...
 'scope','independent double regressions; analytic transfer and exact certificates establish continuous claims');
fid=fopen(fullfile(base,'near_rectangle_matlab.json'),'w');fprintf(fid,'%s\n',jsonencode(out));fclose(fid);disp(out);
function T=Tvalue(v)
 T=4*(capacityR(.5,v/2)+sum(v)-2);
end
function R=capacityR(x,u)
 sx=[1,-1,1,-1];sy=[1,1,-1,-1];W=zeros(5,2);
 for omit=1:5
  members=1:4;if omit<5,members(omit)=[];end
  W(omit,1)=min(sx(members)*x+u(members))-max(sx(members)*x-u(members));
  W(omit,2)=min(sy(members)*(1-x)+u(members))-max(sy(members)*(1-x)-u(members));
 end
 R=sum(W(1:4,1).*W(1:4,2))-2*prod(W(5,:));
end
function e=envelope(x,p)
 L=max(p/2-x,-p/2-1+x);U=min(p/2+x,-p/2+1-x);
 e=8*x.*(1-x)+8*max(L.^2,U.^2);
end
function g=exampleG(x)
 g=.5+zeros(size(x));g(x<1/8)=3*x(x<1/8)+1/8;g(x>7/8)=25/8-3*x(x>7/8);
end
