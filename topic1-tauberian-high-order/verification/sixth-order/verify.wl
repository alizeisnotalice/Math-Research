(* Algebra only; the all-input v^6 theorem is a separate analytic proof. *)
Print["VERSION: ",$Version," MASTER: ",$ProcessID];
Clear[dim,v,x,q,ee,aa,av,ss,mi,loss,dup,rawCheck,physicalCheck];
LaunchKernels[15];Print["KERNELS: ",Length[Kernels[]]," PIDS: ",InputForm[ParallelEvaluate[$ProcessID]]];
rawCheck[dim_]:=Module[{poly,got,want},
 poly=Normal[Series[v^dim/((1+v)^dim-v^dim),{v,0,dim+5}]];
 got=Table[Coefficient[poly,v,dim+j],{j,0,5}];
 want=Table[Sum[(-1)^(j-(kk-1)dim)Binomial[dim+j-1,j-(kk-1)dim],{kk,1,1+Floor[j/dim]}],{j,0,5}];
 {dim,got,got===want,$ProcessID,$KernelID}];
physicalCheck[dim_]:=Module[{poly,sigma,root,ser,co,raw,expected},
 poly=1+(1-q)x^dim/2+(1-q)(2-q)x^(2dim)/6;
 sigma=x^dim poly;root=2x poly^(1/dim);
 ser=Normal[Series[2^dim sigma/((1+root)^dim-2^dim sigma),{x,0,dim+5}]];
 co=Factor/@Table[Coefficient[ser,x,dim+j],{j,0,5}];
 raw=Table[Sum[(-1)^(j-(kk-1)dim)Binomial[dim+j-1,j-(kk-1)dim],{kk,1,1+Floor[j/dim]}],{j,0,5}];
 expected=Switch[dim,2,66-2q,3,-572-4q,4,9224-8q,5,-127984-16q];
 {dim,co,FullSimplify[(co/.q->1)==Table[2^(dim+j)raw[[j+1]],{j,0,5}]],
  FullSimplify[Coefficient[ser,x,2dim]==expected],$ProcessID,$KernelID}];
DistributeDefinitions[rawCheck,physicalCheck];
raw=ParallelTable[rawCheck[dim],{dim,2,16},Method->"CoarsestGrained"];
physical=ParallelTable[physicalCheck[dim],{dim,2,5},Method->"CoarsestGrained"];
Scan[Print["RAW: ",InputForm[#]]&,raw];Scan[Print["PHYSICAL: ",InputForm[#]]&,physical];
f=(1-(1-ee)^q)/q;
flags={FullSimplify[Normal[Series[f,{ee,0,3}]]==ee+(1-q)ee^2/2+(1-q)(2-q)ee^3/6,0<q<=1],
 FullSimplify[Limit[f,q->0,Direction->"FromAbove"]==-Log[1-ee],0<ee<1],
 FullSimplify[D[f,{ee,4}]==(1-q)(2-q)(3-q)(1-ee)^(q-4),0<q<=1&&0<=ee<1],
 Reduce[0<=q<=1&&((1-q)(2-q)(3-q)<0||(1-q)(2-q)(3-q)>6),q,Reals]===False,
 Reduce[0<=q<=1&&((1-q)/2<0||(1-q)(2-q)/6<0),q,Reals]===False,
 Reduce[nn>=2&&jj>=0&&3nn+jj-6<0,{nn,jj},Reals]===False,
 FullSimplify[(1-aa)ss+mi-av mi/aa== (1-aa)(ss-mi/aa)-(av-1)mi/aa,aa>0],
 FullSimplify[(av/aa-loss/aa-1)dup==(av-loss-aa)dup/aa,aa>0],
 FullSimplify[(1+ee)/(1-ee)-1==2ee/(1-ee),0<ee<1],
 FullSimplify[D[Log[(1+ee)/(1-ee)],ee]<4,0<ee<1/2],
 raw[[2,2,6]]===0,
 FullSimplify[physical[[2,2,6]]==160(1-q)]};
Print["SCALAR FLAGS: ",InputForm[flags]];
failures=Count[Join[raw[[All,3]],physical[[All,3]],physical[[All,4]],flags],Except[True]];
Print["FAILURES: ",failures];CloseKernels[];Exit[If[failures==0,0,1]];
