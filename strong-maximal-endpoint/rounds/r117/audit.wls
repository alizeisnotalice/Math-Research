(* Independent exact dot products at the two worst cells, without polynomial convolution. *)
m=1024; l=512; tt=2048; u=512; q=10^6; acc=10^15; jmax=12000;
pow=NestList[Ceiling[1000 #/1001]&, acc, jmax];
fL=Quotient[q (acc-pow),acc+pow]; fE=Quotient[q (acc-pow),acc];
weights[f_]:=Differences[Prepend[f,0]];
idx[i_]:=Ceiling[Range[0,jmax](i+1)/m];
(* At a given y-lattice index, integrate / sum only the remaining X CDF. *)
i=1013; ell=255; z=Range[-u,u-2,2]; ks=Floor[1279 (m+i)(l u + z (ell+Boole[#<0]& /@ z))/(m l u)];
small=Total[Function[k,weights[fL].fL[[1+Clip[k-idx[i],{0,jmax}]]]] /@ ks]-u q^2/2;
i=1004;t=1024;
left=Ceiling[1280(m+i+1)(tt-t)/(m tt)];right=Floor[(m+i)(1000tt+1279t)/(m tt)];
rup=Ceiling[1280(m+i+1)(t+1)/(m tt)];
prim=Prepend[Accumulate[fE],0]; ids=idx[i];
large=weights[fE].(prim[[1+Clip[right-ids,{0,jmax}]]]-prim[[1+Clip[left-ids,{0,jmax}]]])-rup q^2;
(* Rational bracket for the implicit root and logarithmic mesh constant. *)
a=12784/10000;b=12785/10000;
expUpper=Sum[a^k/k!,{k,0,12}]+a^13/13!/(1-a/14);
expLower=Sum[b^k/k!,{k,0,12}];
checks={small==1090807670824,large==276296984034,
 (a-1)expUpper<1,(b-1)expLower>1,
 1279(1/1000-1/(2 10^6)+1/(3 10^9))<a,
 b<1280(1/1000-1/(2 10^6))};
Print[InputForm[{small,large,checks}]];If[!And@@checks,Exit[1]];
