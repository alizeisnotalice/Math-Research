ClearAll[t,s,r,q]; r=1+2s*t+q*t^2;v1=r(1-s);v2=r*s;
end=((1-t*r)(1-2t*r*s)-t^2(2-2r*s))/((1-2t*r*s)Sqrt[1-4t^2]);
a=2v1-2;b1=v1-a/2;b2=v2-a/2;u1=1;u2=2-v2+a/2;
corner=((1-t*b1)(1-t*b2)-t^2*u1*u2)/((1+t*a)Sqrt[1-4t^2]);
e=FullSimplify[Coefficient[Normal[Series[end-(1-t),{t,0,3}]],t,3]];
c=FullSimplify[Coefficient[Normal[Series[corner-(1-t),{t,0,3}]],t,3]];
qe=q/.First[Solve[e==0,q]];qc=q/.First[Solve[c==0,q]];
Print[<|"end_residual"->e,"corner_residual"->c,"end_radial_q"->qe,"corner_radial_q"->qc,"difference"->FullSimplify[qe-qc],"quartic"->FullSimplify[8Integrate[qe,{s,0,1/2}]+4Integrate[(2s)^2,{s,0,1/2}]+8],"version"->$Version|>];
