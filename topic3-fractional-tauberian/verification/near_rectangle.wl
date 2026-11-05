Clear[x,p,mu,h,u,v,vx,vy,rho,t];
tfun[vv_]:=Module[{aa,bb},aa=Min[vv[[2]],vv[[4]]]+Min[vv[[1]],vv[[3]]];bb=Min[vv[[3]],vv[[4]]]+Min[vv[[1]],vv[[2]]];Sum[(vv[[BitXor[j,2]+1]]+Min[vv[[BitXor[j,1]+1]],vv[[BitXor[j,3]+1]]])(vv[[BitXor[j,1]+1]]+Min[vv[[BitXor[j,2]+1]],vv[[BitXor[j,3]+1]]]),{j,0,3}]-2aa bb];
coeff={};shiftChecks={};
Do[
 V=Table[ReplacePart[ConstantArray[0,4],Thread[order->Join[ConstantArray[-1,k],ConstantArray[1,4-k]]]],{k,0,4}];
 q[vv_]:=8-tfun[vv];
 Do[AppendTo[coeff,If[i==j,q[V[[i]]],2q[(V[[i]]+V[[j]])/2]-(q[V[[i]]]+q[V[[j]]])/2]],{i,1,5},{j,i,5}];
 vv=Array[v,4];pos=Ordering[order];mn[i_,j_]:=If[pos[[i+1]]<pos[[j+1]],vv[[i+1]],vv[[j+1]]];
 aa=mn[1,3]+mn[0,2];bb=mn[2,3]+mn[0,1];
 tt=Sum[(vv[[BitXor[j,2]+1]]+mn[BitXor[j,1],BitXor[j,3]])(vv[[BitXor[j,1]+1]]+mn[BitXor[j,2],BitXor[j,3]]),{j,0,3}]-2aa bb;
 wx=-2x+2mu+h aa;wy=-2(1-x)+2mu+h bb;
 rpoly=Sum[(-2x+2mu+h(vv[[BitXor[j,2]+1]]+mn[BitXor[j,1],BitXor[j,3]]))(-2(1-x)+2mu+h(vv[[BitXor[j,1]+1]]+mn[BitXor[j,2],BitXor[j,3]])),{j,0,3}]-2wx wy;
 AppendTo[shiftChecks,Expand[rpoly+8mu+2h Total[vv]-(8x(1-x)+h^2 tt+4mu h Total[vv]+8mu^2)]===0],
 {order,Permutations[Range[4]]}];
a=(1-p)/2;
envelope=FullSimplify[2(Integrate[8x(1-x)+8(p/2+x)^2,{x,0,a}]+Integrate[8x(1-x)+8(1-p/2-x)^2,{x,a,1/2}]),Assumptions->0<=p<=1];
maxCheck=Resolve[ForAll[p,0<=p<=1,2+2p(1-p)<=5/2],Reals];
dmin=TimeConstrained[Resolve[ForAll[{vx,vy,rho},vx>=0&&vy>=0,Exists[u,Max[vx,Abs[u-rho/2]]+Max[vy,Abs[u+rho/2]]==Max[vx+vy,Abs[rho]]]],Reals],20,"timeout"];
alpha=1-t;w=Exp[-rho t];s=(w+1/w)/2;reference=FullSimplify[Normal[Series[2(s-alpha)(1/alpha-s),{t,0,5}]]];
sixthGap=Normal[Series[(1/alpha-alpha)^2/2-2(Cosh[t]-alpha)(1/alpha-Cosh[t]),{t,0,8}]];
sixthCheck=Expand[sixthGap-(t^6/2+11t^7/12+409t^8/288)]===0;
referenceCheck=Expand[reference-(2t^2+2t^3+(2+rho^2-rho^4/2)t^4+(2+rho^2)t^5)]===0;
If[Length[coeff]!=360||Min[coeff]<0||!And@@shiftChecks||envelope!=2+2p-2p^2||!TrueQ[maxCheck]||!referenceCheck||!sixthCheck||!TrueQ[dmin],Print[{Length[coeff],Min[coeff],shiftChecks,envelope,maxCheck,referenceCheck,dmin}];Exit[1]];
out=<|"kernel"->$Version,"bernstein_coefficients"->360,"coefficient_histogram"->Association[(ToString[#[[1]]]->#[[2]])&/@SortBy[Tally[coeff],First]],"order_chamber_shift_identity_checks"->Length[shiftChecks],"envelope_integral"->ToString[envelope,InputForm],"envelope_maximum_check"->maxCheck,"first_defect_minimum_quantifier_check"->dmin,"reference_series"->ToString[reference,InputForm],"rho_one_gap_series"->ToString[sixthGap,InputForm],"critical_rho_series"->ToString[Series[-Log[1-t]/t,{t,0,4}],InputForm]|>;
Export[FileNameJoin[{DirectoryName[$InputFileName],"..","results","near_rectangle_wolfram.json"}],out,"RawJSON"];Print[out];
