ClearAll[t,z,a,r,q,ep,M,g];
a=1-t;r=t/2+t^2/2+15t^3/16+9t^4/4;q=a+t^2/4;ep=t^3/16;M=2q-a;
g[z_]:=(z-2r)^2-2a z+1;
margins={g[a]-4ep(a-2r)-4ep^2,-g[q]-4ep(q-2r)-4ep^2,
 -g[1/a]-4ep(1/a-2r)-4ep^2,(a+2r)^2-4M(M-a)-1-16a ep};
polys=Cancel/@(margins/{t^5,t^4,t^3/a^2,t^3});
lower=Table[Coefficient[p,t,0]+Sum[Min[Coefficient[p,t,k],0]/100^k,{k,1,Exponent[p,t]}],{p,polys}];
If[!And@@Thread[lower>0],Exit[1]];
If[!TrueQ[FullSimplify[And[0<2r+2ep<a,a<q<(a+1/a)/2],0<t<=1/100]],Exit[2]];
Print[$Version];Print["normalized polynomials = ",polys];Print["whole-interval positive lower bounds = ",lower];
Print["Fourier reciprocal sum = ",Sum[1/((2m+1)^2-1),{m,1,Infinity}]];
Print["gap coefficient = ",FullSimplify[ep^2/(3t^6)]];
(* General template: r=r*-d t^3, q=a+e t^2, epsilon=h t^3. *)
rs=(a-Sqrt[2a^2-1])/2;rr=rs-d t^3;qq=a+e t^2;ee=h t^3;
gg[z_]:=(z-2rr)^2-2a z+1;
lead=FullSimplify[Limit[#/t^3,t->0,Direction->"FromAbove"]]&/@{gg[a]-4ee(a-2rr)-4ee^2,-gg[qq]-4ee(qq-2rr)-4ee^2,-gg[1/a]-4ee(1/a-2rr)-4ee^2};
Print["template leading margins = ",lead];Print["template mass leading = ",Limit[((a+2rr)^2-4(2qq-a)(2qq-2a)-1-16a ee)/t^2,t->0,Direction->"FromAbove"]];
Print["moment beta coefficient = ",FullSimplify[Factorial[6]/768]];
