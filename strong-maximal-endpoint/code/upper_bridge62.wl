ClearAll[c,s,a,y];
den=Integrate[2,{y,0,c}]+Integrate[4,{y,c,s-c}]+Integrate[2,{y,s-c,s}];
num=Integrate[2,{y,0,a}];
ratio=FullSimplify[num/den,Assumptions->0<a<c<s/2];
componentIdentity=FullSimplify[(n p-(n-1)k/g)^2-n p(1-p)-(n-1)k^2/g^2-(n p((n+1)p-1)-2n p(n-1)k/g+(n-1)(n-2)k^2/g^2)];
checks={componentIdentity==0,FullSimplify[den==4(s-c)],FullSimplify[ratio==a/(2(s-c))],FullSimplify[D[ratio,c]==a/(2(s-c)^2)]};
den2=Integrate[2,{y,0,1}]+Integrate[1,{y,1,2}]+Integrate[2,{y,2,3}];
Print[<|"version"->$Version,"symbolic_checks"->checks,"c1"->(ratio/.{a->1/2,s->3,c->1}),"c2"->1/den2,"scope"->"Exact conditional radialization counterexample; not endpoint-efficient"|>];
If[!And@@checks,Exit[1]];
