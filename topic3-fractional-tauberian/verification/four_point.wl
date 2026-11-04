(* All inputs exact; no numerical minimizer used as proof. *)
ClearAll[a,r,q,z,g];
g[z_] := (z-2r)^2-2a z+1;
ex={a->9/10,r->9/160,q->901/1000};
vals={g[a],g[q],g[1/a],(a+2r)^2-4(2q-a)(2q-2a)-1}/.ex;
If[!And@@Thread[Sign[vals]=={1,-1,-1,1}],Exit[1]];
rs=(a-Sqrt[2a^2-1])/2;
identities=FullSimplify[{g[a]/.r->rs,(D[g[z],z]/.z->a)/.r->rs,
 (g[1/a]-g[a]-(1/a-a)*(1/a-a-4r))/.r->rs},1/Sqrt[2]<a<1];
Print[$Version];Print["exact margins = ",vals];Print["endpoint identities = ",identities];
Print["union excess sign = ",FullSimplify[2a-Sqrt[2a^2-1]>1,1/Sqrt[2]<a<1]];
eps=1/100000;
robust=({g[a]-4eps(a-2r)-4eps^2,-g[q]-4eps(q-2r)-4eps^2,
 -g[1/a]-4eps(1/a-2r)-4eps^2,(a+2r)^2-4(2q-a)(2q-2a)-1-16a eps}/.ex);
If[!And@@Thread[robust>0],Exit[2]];Print["robust margins = ",robust];
