(* Actual one-dimensional catalogue: [0,3/2] and [-1/2,1]. *)
ker={{0,2/3,1/3},{1/3,2/3,0}};
ref={1/4,1/2,1/4};
pi={1/2,1/2};
nu={{3/8,1/8},{1/8,3/8}};
ent[p_]:=-Total[(If[#==0,0,# Log[#]]&)/@p];
kl[p_,q_]:=Total[MapThread[If[#1==0,0,#1 Log[#1/#2]]&,{p,q}]];
nuA=Total/@nu;
outA=nuA.ker;
outAB=Transpose[ker].nu.ker;
condB=Table[(Transpose[ker].nu)[[z]]/outA[[z]],{z,3}];
fAB=kl[Flatten[outAB],Flatten[Outer[Times,ref,ref]]]-kl[Flatten[nu],Flatten[Outer[Times,pi,pi]]];
fA=kl[outA,ref]-kl[nuA,pi];
fBgivenA=Sum[outA[[z]](kl[condB[[z]].ker,ref]-kl[condB[[z]],pi]),{z,3}];
cmi=Sum[outA[[z]]ent[condB[[z]]],{z,3}]-Sum[nuA[[i]]ent[nu[[i]]/nuA[[i]]],{i,2}];
lambda=Log[4/(3 2^(1/3))];
Print[FullSimplify[fAB-fA-fBgivenA+cmi]];
Print[N[{fAB,2lambda,cmi,2lambda-fAB},30]];
