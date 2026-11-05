root=DirectoryName[DirectoryName[$InputFileName]];
Clear[z];
v1[z_]:=2/(1+Exp[z]);
g1[z_]:=Exp[-z](1+2Log[(1+Exp[z])/2]);
identity=FullSimplify[D[g1[z],z]-v1[z]+g1[z],Assumptions->z>0];
files={"bellman_005.json","bellman_0025.json","bellman_00125.json"};
checks=Map[Function[file,Module[{data,row,q,lo=0,hi,c,r,res,value},
 data=Import[FileNameJoin[{root,"rounds","r011",file}],"RawJSON"];
 row=First[data["rows"]]; q=N[Rationalize[row["probe_q"],0],80]; hi=q;
 Do[c=(lo+hi)/2;r=Exp[-c];
 res=g1[q]+r g1[q-c]-(1+r)v1[q-c];
 If[res>0,lo=c,hi=c],{200}];
 c=(lo+hi)/2;r=Exp[-c];value=(g1[q]+r g1[q-c])/(1+r);
 <|"step"->data["step"],"two_step_hit80"->ToString[N[value,30],InputForm],
 "first_control_r80"->ToString[N[r,30],InputForm],
 "grid_hit_error"->N[Abs[1-row["adaptive_failure"]-value]]|>
]],files];
mu=1+ProductLog[1/E]; aa=Log[2];
varianceChecks=Table[<|"n"->n,"variance"->N[(n-1)(1+mu)/2+1-aa^2+((n-1)(mu-1)+2aa)^2/4,20]|>,{n,{16,64,256,1024}}];
result=<|"wolfram_version"->$Version,"convolution_ode_residual"->ToString[identity,InputForm],
 "checks"->checks,"branch_policy_variances"->varianceChecks|>;
Export[FileNameJoin[{root,"rounds","r011","crosscheck.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
