ClearAll[a,u,j,t,p,x];c=2a^2-1;w=1/a-a;
int=Integrate[a-c/(a-u),{u,0,w},Assumptions->1/Sqrt[2]<a<1];upper=1-a^2+c Log[c/a^2];
co=FullSimplify[2(2^(j-1)-1)/(j-1)-(2^j-1)/j-(2^j-j-1)/(j(j-1)),Assumptions->Element[j,Integers]&&j>=2];
Print[<|"quadrant_integral_residual"->FullSimplify[int-upper,Assumptions->1/Sqrt[2]<a<1],"coefficient_identity"->co,"upper_series"->ToString[Series[upper/.a->1-t,{t,0,8}],InputForm],"lower_series"->ToString[Series[((1/a-a)^2/2)/.a->1-t,{t,0,6}],InputForm],"moment_upper"->ToString[Series[4/((p+1)(p+2))+20/((p+1)(p+2)(p+3)),{p,Infinity,3}],InputForm],"moment_lower"->ToString[Series[4/((p+1)(p+2))+12/((p+1)(p+2)(p+3)),{p,Infinity,3}],InputForm],"version"->$Version|>];
