Clear[t,x,m,g,t0,bb,aa];
f=Exp[-x^2/(2m t)]/Sqrt[2Pi m t];
z2=x^2/(m t);
h=g*((1-t0/t)*(z2-1)+Log[t/t0]);
u=-h*f;
res=FullSimplify[D[u,t]-m/2 D[u,{x,2}]+g*z2/t*f,Assumptions->{m>0,t>t0>0,Element[x,Reals]}];
initial=FullSimplify[u/.t->t0,Assumptions->{m>0,t0>0,Element[x,Reals]}];
hsimple=bb*z2+aa*Log[t/t0];
firstResidual=FullSimplify[(D[-hsimple*f,t]-m/2 D[-hsimple*f,{x,2}])/f,Assumptions->{m>0,t>0,t0>0,Element[x,Reals]}];
mu=1+ProductLog[Exp[-1]];
nums=Table[With[{zz=z/.FindRoot[Erfc[z/Sqrt[2]]/2==Sqrt[mu]*z/Sqrt[n],{z,Sqrt[Log[n]]},WorkingPrecision->70]}, {n,N[zz,30],N[(mu-1)/8*((1-1/n)*(zz^2-1)+Log[n]),30]}],{n,{64,256,1024,2048}}];
angular=FullSimplify[Limit[ll*((((ll-qq)*alpha+xi+zz)/(ll-qq+tt))-alpha),ll->Infinity]-(xi+zz-alpha*tt)];
(* Exactly positive angular crossing probability in n=2, theta1=1/4. *)
c=mu;pr=1/(1+Exp[c]);
laplaceTail[v_]:=If[v>=0,Exp[-v]/2,1-Exp[v]/2];
angleRows=Table[With[{vv=ll/2}, {ll,N[Sum[(If[b1==1,pr,1-pr]) (If[b2==1,pr,1-pr]) laplaceTail[vv-c*(b1-b2)],{b1,0,1},{b2,0,1}],50]}],{ll,{10,20,40,80}}];
r=<|"angular_expansion_residual"->ToString[angular,InputForm],"angular_crossing_rows"->Map[ToString[#,InputForm]&,angleRows,{2}],"version"->$Version,"linear_heat_residual"->ToString[res,InputForm],"initial_residual"->ToString[initial,InputForm],"simple_shift_heat_residual_over_f"->ToString[firstResidual,InputForm],"precision"->70,"linear_model_predictions"->Map[ToString[#,InputForm]&,nums,{2}]|>;
base=DirectoryName[DirectoryName[$InputFileName]];Export[FileNameJoin[{base,"rounds","r025","symbolic_checks.json"}],r,"RawJSON"];Print[r];
