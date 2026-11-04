Print["Version: ", $Version];
js={0,1,3,4,5,6,8,9}; g1={0,1,8,9}; g2={3,4,5,6};
(* Exact cubic primitive, with each linear positive part resolved over 0<e<1. *)
Clear[pos]; pos[z_] := pos[z] = Module[{v=Together[z], roots, sign}, roots=Reduce[0<e<1 && v==0,e,Reals]; If[roots=!=False && v=!=0, Print["INTERIOR ROOT: ",InputForm[v]," ",InputForm[roots]]; Abort[]]; sign=Sign[v/.e->1/2]; If[sign>0,v,0]];
pr[z_,rr_]:=(pos[z+rr]^3-2pos[z]^3+pos[z-rr]^3)/(6rr);
pair[i_,j_,rr_,k_]:= Module[{a=i/4,b=(i+e)/4,c=j/4,d=(j+e)/4}, pr[b-c-3k,rr]-pr[a-c-3k,rr]-pr[b-d-3k,rr]+pr[a-d-3k,rr]];
e1=Expand[Sum[pair[i,j,1,k],{i,g1},{j,js},{k,-1,1}]]; e2=Expand[Sum[pair[i,j,1/4,k],{i,g2},{j,js},{k,-1,1}]];
Print["E1=",InputForm[e1],"; E1/4=",InputForm[e2]]; Print["exact checks: ",InputForm[{FullSimplify[e1==19e^2/32,0<e<1],FullSimplify[e2==e^2/4-e^3/48,0<e<1]}]];
(* Independent direct Integrate for representative self-pair and adjacent-pair. *)
Print["direct self 1/4=",InputForm[Integrate[1-4 Abs[x-y],{x,0,e/4},{y,0,e/4},Assumptions->0<e<1,GenerateConditions->False]]];
Print["direct adjacent 1/4=",InputForm[Integrate[4 (x-w),{x,0,e/4},{w,0,x},Assumptions->0<e<1,GenerateConditions->False]]];
