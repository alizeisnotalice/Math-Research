Print["VERSION ",$Version," PID ",$ProcessID];
(* Capacity for 0<=s<=S, kink t=h+m/C: derivatives suffice for the global maximum. *)
cap=aa((R+s)^n-R^n);Print["capacity derivative=",InputForm[FullSimplify[D[cap,s],Element[n,Integers]&&n>=2&&aa>0&&R>0&&s>=0]]];
Print["pre-kink derivative positive=",InputForm[FullSimplify[D[cap,s]>0,Element[n,Integers]&&n>=2&&aa>0&&R>0&&s>=0]]];
Print["post-kink derivative negative=",InputForm[FullSimplify[D[cap,s]-C<0,Element[n,Integers]&&n>=2&&aa>0&&R>0&&S>=s>=0&&C>n aa(R+S)^(n-1)]]];
Print["power monotonic derivative=",InputForm[FullSimplify[D[(R+s)^(n-1),s]>=0,Element[n,Integers]&&n>=2&&R>0&&s>=0]]];
Print["order transitivity violation=",InputForm[Reduce[0<=lo<=hi&&hi<C&&lo>=C,{lo,hi,C},Reals]]];
(* Linear prototype, interior kink. *)
tt=h+m/C;linear0=L s-Max[0,C(s-h)-m];linearT=L tt;
Print["linear before kink=",InputForm[FullSimplify[linearT-(L s),C>0&&L>0&&0<=s<=tt]]];
Print["linear after kink=",InputForm[FullSimplify[linearT-(L s-C(s-h)+m),C>L>0&&s>=tt]]];
Print["linear pre-kink violation=",InputForm[Reduce[C>L>0&&h>=0&&m>=0&&0<=s<=h+m/C&&L(h+m/C)<L s,{C,L,h,m,s},Reals]]];
Print["linear post-kink violation=",InputForm[Reduce[C>L>0&&h>=0&&m>=0&&s>=h+m/C&&L(h+m/C)<L s-C(s-h)+m,{C,L,h,m,s},Reals]]];
(* Two-arm lower estimate: subtract two 8h² errors and at most h² overlap. *)
arm1=Integrate[h(1-y/e)-8h^2/e,{y,0,e(1-z)}];arm2=Integrate[h(1-y/e)-8h^2/e,{y,0,e}];
expected1=e h((1-z)-(1-z)^2/2)-8h^2(1-z);Print["first arm identity=",InputForm[FullSimplify[arm1==expected1,e>0]]];
Print["second arm identity=",InputForm[FullSimplify[arm2==e h/2-8h^2,e>0]]];
Print["two arm slack=",InputForm[FullSimplify[arm1+arm2-h^2-e h(1-z^2/2-17h/e),e>0]]];
Print["positive arm bound violation=",InputForm[Reduce[0<del<1/35&&0<=z<=1-2del&&1-z^2/2-17del<=0,{del,z},Reals]]];
Z=1-2del;integ=Integrate[1-z^2/2-17del,{z,0,Z}];Print["recovery integral=",InputForm[Factor[integ]]];Print["zero-delta integral=",InputForm[integ/.del->0]];
bar=8n/3;Kstar=B n^2/27;beta=5B bar^4/(96n^3);reccoef=4B (1/4)bar^4/(16n^2) (5/6);
Print["recovery coefficient=",InputForm[FullSimplify[reccoef==5B bar^4/(96n^2),n>=2]]];
Print["coefficient improvement=",InputForm[FullSimplify[2n beta/2^7==(10/9)Kstar,n>=2]]];
Print["resultant shift=",InputForm[FullSimplify[-(n+4)Kstar+(10/9)Kstar==-(n+26/9)Kstar,n>=2]]];
