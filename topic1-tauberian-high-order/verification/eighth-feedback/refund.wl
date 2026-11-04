ClearAll[d,q,x,n,e,f];
Print["Wolfram ",$Version];
(* 0<d<1/2, theta=d/(1-d), intersection of unit-square two arms. *)
th=d/(1-d); mlo=q-q^2/(2th); mhi=FullSimplify[(d-d^2/(2th))+th((q-q^2/2)-(d-d^2/2))];
Print["M0 low, high: ",{mlo,mhi}];
I0=th(5+3th)/(3(1+th)^2);
Print["I0 exact simplification: ",FullSimplify[I0,0<d<1/2]];
weighted[q_,m_]:=Integrate[1+q-2x,{x,d(q-m),d q}];
diff=FullSimplify[2 Integrate[weighted[q,mlo],{q,0,d}]+2 Integrate[weighted[q,mhi],{q,d,1}],0<d<1/2];
Print["I0-I1 exact: ",Factor[diff]];
Print["I0-I1 series: ",Series[diff,{d,0,3}]];
Print["I1 series: ",Series[I0-diff,{d,0,3}]];
Print["triple integral: ",Integrate[(e-x)(e-y)/f,{x,0,e},{y,0,e},Assumptions->e>0&&f>0]];
