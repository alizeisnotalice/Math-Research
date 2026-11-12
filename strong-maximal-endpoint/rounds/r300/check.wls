Clear[t,s,q,b,u,v,k];
z[t_]:=(v[t]-b)/u;
h[t_]:=Erfc[z[t]/Sqrt[2]]/2;
checks={FullSimplify[h''[t]/h'[t]==-z[t] v'[t]/u+v''[t]/v'[t],Assumptions->u>0&&v'[t]<0],
FullSimplify[D[k[s Log[q]],{q,2}]==s^2/q^2 (k''[s Log[q]]-k'[s Log[q]]/s),Assumptions->q>0&&s>0]};
Print[checks];If[!And@@checks,Exit[1]];
