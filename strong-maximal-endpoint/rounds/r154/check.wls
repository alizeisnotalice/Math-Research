cases=Flatten[Table[{m,t,e},{m,{2,4,8}},{t,{1/4,3/4}},{e,{0,1/8,1/4,1/2,3/4}}],2];
test[{m_,t_,e_}]:=Module[{v},v=MinValue[{x+y+m(z+w),0<=x<=t&&0<=y<=1-t&&0<=z<=1-t&&0<=w<=t&&x+y+z+w>=1&&t-x+y+1-t-z+w<=2e},{x,y,z,w},Reals];v==m-(m-1)Min[1,t+e]];
tests=test/@cases;Print[{Length[tests],And@@tests}];If[!TrueQ[And@@tests],Exit[1]];
