run[m_,d_,fac_]:=Module[{a=1+d,p,l,t,e,price},p=a/2;l=d^2/(2a);t=p-l;e=fac*d^2/a^2;
price=MinValue[{x+y+m(z+w),0<=x<=t&&0<=y<=1-t&&0<=z<=l&&0<=w<=d-l&&x+y+z+w>=p&&t-x+y+l-z+w<=e a},{x,y,z,w},Reals];price==p+(m-1)Max[0,l-e a/2]];
tests=Flatten[Table[run[m,d,f],{m,{2,8,32}},{d,{1/4,1/32}},{f,{1/2,1,2}}]];Print[{Length[tests],And@@tests}];If[!TrueQ[And@@tests],Exit[1]];
