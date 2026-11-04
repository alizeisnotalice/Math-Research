(* Identities in the analytic proof, not a numerical proof of the operator bound. *)
Clear[nn,gg,tt,qq,zz,kk];
Print["Wolfram version: ",$Version];
ea=Integrate[zz^(nn-1)+(nn-1)/(nn-1+gg)(zz^(-gg)-zz^(nn-1)),{zz,0,1},Assumptions->nn>1&&0<gg<1];
dfun=qq^nn/(nn Log[qq]);
checks={
 TrueQ[FullSimplify[ea==(nn-gg)/(nn(1-gg)),nn>1&&0<gg<1]],
 TrueQ[1+Sum[2^kk 2^(-3kk-4)/8,{kk,1,Infinity}]==385/384],
 TrueQ[(1+Sqrt[1/4])/((1/2)(1+1/4))==12/5>2],
 TrueQ[FullSimplify[D[dfun,qq]==qq^nn(nn Log[qq]-1)/(nn qq Log[qq]^2),nn>0&&qq>1]],
 TrueQ[FullSimplify[dfun/.qq->Exp[1/nn],nn>0]==E],
 TrueQ[FullSimplify[Limit[(tt^tt(1+tt)^(nn-tt)-1)/(tt Log[1/tt]),tt->0,Direction->"FromAbove"]==-1,nn>0]],
 TrueQ[FullSimplify[Limit[(tt^tt(nn-tt)^(nn-tt)/nn^nn-1)/(tt Log[1/tt]),tt->0,Direction->"FromAbove"]==-1,nn>0]]
};
finite=Flatten[Table[
 TrueQ[Sum[(q0+2q0^(-j))^n0,{j,0,count-1}]==q0^n0 count+Sum[Binomial[n0,k]q0^(n0-k)2^k(1-q0^(-k count))/(1-q0^(-k)),{k,1,n0}]],
 {n0,1,5},{q0,{3/2,2}},{count,{1,2,7}}]];
Print["analytic guards: ",checks];Print["exact finite sum controls: ",Length[finite],"; all=",And@@finite];
If[!And@@Join[checks,finite],Exit[1]];Print["ANGULAR/SHELL GUARDS PASS"];
