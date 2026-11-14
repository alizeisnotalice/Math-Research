folder=DirectoryName[$InputFileName];
verify[k_]:=Module[{n=k^2,th,b,a,ua,ub,u,vA,vars,c,x,y,z,w,l,cost,avA,avB,true,resp,target,g,diff},
 th=2^(-n);b=(3k-1)/2;a=1+2/k;ua=(k+1)^2;ub=6(k+1);u=ua+ub-th;vA=2+2/k;vars={c,x,y,z,w,l};cost=c+x+y+n z+b w;avA=(c+x+z)/vA;avB=(c+y+w)/12;
 true=First[Minimize[{cost,0<=c<=th&&0<=x<=1-th&&0<=y<=8-th&&0<=z<=a&&0<=w<=4&&avA>=1/2&&avB>=1/2},{c,x,y,z,w}]];
 resp=First[Minimize[{cost,0<=c<=th&&0<=x<=1-th&&0<=y<=8-th&&0<=z<=a&&0<=w<=4&&0<=l<=1&&l<=avA&&l<=avB&&(ua-th)avA+(ub-th)avB+th l>=u/2},vars]];
 target=8-th+th/(6(k+1)-th);g=u/(k+7-th);diff=g-(k+1);
 <|"k"->k,"n"->n,"true_source"->ToString[true,InputForm],"responsible_source"->ToString[resp,InputForm],"true_pass"->(true==k+7-th),"responsible_pass"->(resp==target),"strict_score_difference_positive"->(diff>0),"score_identity"->(diff==k th/(k+7-th)),"core_mass"->ToString[th,InputForm]|>];
ks={3,5,7,9,11,13,15,17,19,25,49,101};LaunchKernels[12];DistributeDefinitions[verify];res=ParallelMap[Function[k,TimeConstrained[MemoryConstrained[verify[k],268435456,$Failed],120,$Aborted]],ks];CloseKernels[];
Export[FileNameJoin[{folder,"wolfram.json"}],<|"version"->$Version,"results"->res,"scope"->"Exact five source classes from membership-atom permutation symmetry; no ambient coordinate grid enumerated"|>,"RawJSON"];Print[<|"cases"->Length[res],"passed"->Count[res,r_Association/;TrueQ[r["true_pass"]&&r["responsible_pass"]&&r["score_identity"]&&r["strict_score_difference_positive"]]]|>];Quit[];
