(* Independent exact regeneration of the quadratic simplex certificate. *)
Clear[x,t,a,b,c,d,A,v];
rr[xx_,u_]:=Module[{yy=1-xx,wx,wy},
 wx=-2xx+Min[u[[2]],u[[4]]]+Min[u[[1]],u[[3]]];
 wy=-2yy+Min[u[[3]],u[[4]]]+Min[u[[1]],u[[2]]];
 Sum[(-2xx+u[[BitXor[j,2]+1]]+Min[u[[BitXor[j,1]+1]],u[[BitXor[j,3]+1]]])
 (-2yy+u[[BitXor[j,1]+1]]+Min[u[[BitXor[j,2]+1]],u[[BitXor[j,3]+1]]]),{j,0,3}]-2wx wy];
loss[xx_,vv_]:=(8xx-2Total[xx vv]-rr[xx,xx vv])/xx;
coeff={};affinity={};
Do[
 vertices=Table[ReplacePart[ConstantArray[0,4],Thread[order->Join[ConstantArray[-1,k],ConstantArray[1,4-k]]]],{k,0,4}];
 Do[
  q[vv_]:=If[endpoint==0,2loss[1/4,vv]-loss[1/2,vv],loss[1/2,vv]];
  Do[AppendTo[coeff,If[i==j,q[vertices[[i]]],2q[(vertices[[i]]+vertices[[j]])/2]-(q[vertices[[i]]]+q[vertices[[j]]])/2]],{i,1,5},{j,i,5}],
  {endpoint,{0,1/2}}];
 (* On this order chamber each Min is the first ordered variable; verify x-affinity. *)
 vv=Array[v,4];pos=Ordering[order];
 mn[i_,j_]:=If[pos[[i+1]]<pos[[j+1]],vv[[i+1]],vv[[j+1]]];
 wx=-2x+x(mn[1,3]+mn[0,2]);wy=-2(1-x)+x(mn[2,3]+mn[0,1]);
 poly=Sum[(-2x+x(vv[[BitXor[j,2]+1]]+mn[BitXor[j,1],BitXor[j,3]]))(-2(1-x)+x(vv[[BitXor[j,1]+1]]+mn[BitXor[j,2],BitXor[j,3]])),{j,0,3}]-2wx wy;
 AppendTo[affinity,Expand[D[Cancel[(8x-2x Total[vv]-poly)/x],{x,2}]]===0],
 {order,Permutations[Range[4]]}];
rank=FullSimplify[Total[Min@@#& /@Subsets[{a,b,c,d},{3}]]-2Min[a,b,c,d]==a+b&&Total[Max@@#& /@Subsets[{a,b,c,d},{3}]]-2Max[a,b,c,d]==c+d,Assumptions->a<=b<=c<=d];
coverage=And@@Table[Max[k-2,0]==Binomial[k,3]-2Binomial[k,4],{k,0,4}];
alpha=1-t;w=Sqrt[1-A t^2];s=(w+1/w)/2;
strip=FullSimplify[Normal[Series[2(s-alpha)(1/alpha-s),{t,0,6}]],Assumptions->A>=0];
stripCheck=Expand[strip-(2t^2+2t^3+2t^4+2t^5+(2+A^2/4)t^6)]===0;
If[Length[coeff]!=720||Min[coeff]<0||!And@@affinity||!TrueQ[rank]||!coverage||!stripCheck,Print[{Length[coeff],Min[coeff],affinity,rank,coverage,strip}];Exit[1]];
out=<|"kernel"->$Version,"bernstein_coefficients"->Length[coeff],"minimum"->Min[coeff],"histogram"->Association[Rule@@@SortBy[Tally[coeff],First]],"order_chamber_affinity_checks"->Length[affinity],"rank_identity"->rank,"coverage_identity"->coverage,"thin_strip_series"->ToString[strip,InputForm],"thin_strip_series_check"->stripCheck|>;
(* JSON requires string histogram keys. *)
out["histogram"]=Association[(ToString[#[[1]]]->#[[2]])&/@SortBy[Tally[coeff],First]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"..","results","near_square_wolfram.json"}],out,"RawJSON"];Print[out];
