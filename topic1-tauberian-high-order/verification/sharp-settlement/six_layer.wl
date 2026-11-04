Print["VERSION ",$Version," PID ",$ProcessID];
tau=b v/cc;lhs=Sum[tau^(6-h)(cc f[h]-b v f[h-1]),{h,1,6}];rhs=cc f[6]-b v tau^5 f[0];
Print["weighted telescope=",InputForm[FullSimplify[lhs==rhs,cc>0]]];
W=Sum[tau^k,{k,0,5}];Print["weight exact=",InputForm[Expand[W]]];Print["coarse weights divided v=",InputForm[Factor[(W-1)/v]]];
Print["endpoint identity=",InputForm[FullSimplify[b v tau^5 A0 v M/a==A0 b^6 v^7 M/(cc^5 a),cc>0&&a>0]]];
(* A represents A_v, k=K+epsilon; x=M/E, y=S/E. *)
h1=a(1+k v^6 y)-(A-a)x;h2=1-(1-a)y;
gap=A/(A-a)+a k v^6/((A-a)(1-a))-(1+x);
Print["ratio slack identity=",InputForm[FullSimplify[Together[gap-(h1/(A-a)+a k v^6 h2/((A-a)(1-a)))]]]];
Print["positive denominator violation=",InputForm[Reduce[0<a<1&&A>a&&(A-a)(1-a)<=0,{a,A},Reals]]];
Print["ratio slack violation=",InputForm[Reduce[0<a<1&&A>a&&k>=0&&v>0&&h1>=0&&h2>=0&&gap<0,{a,A,k,v,x,y},Reals]]];
Print["prefactor limit=",InputForm[Limit[1/(((1+v)^n-v^n)(1-v^n)),v->0,Assumptions->Element[n,Integers]&&n>=2]]];
