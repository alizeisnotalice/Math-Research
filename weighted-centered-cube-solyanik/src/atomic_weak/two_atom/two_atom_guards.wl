(* Exact independent guards for tau=1/10 rational two-atom witness. *)
Print["wolfram_version=",$Version];
Clear[x];
radNeg=1/8;radPos=15/16;radBoth=17/16;
blocks={{-1-radNeg,-1+radNeg},{1-radBoth,-1+radBoth},{1-radPos,1+radPos}};
a=5059/5000;b=2467/2500;c=3789/5000;d=10683/10000;
twoLower=(8/17)(a-b+c+d);
A0=535887/500000;B0=53589/1000000;d0=173/100-A0;c0=173/100-B0;s0=4361/100000;
guards={
 TrueQ[blocks=={{-9/8,-7/8},{-1/16,1/16},{1/16,31/16}}],
 TrueQ[a^10<9/8&&b^10>7/8&&c^10<1/16&&d^10<31/16],
 TrueQ[twoLower==18511/21250&&twoLower>87/100],
 TrueQ[A0^10>2&&(1/(10B0))^10<2^9&&s0^9>(1/(10*c0))^10&&s0<(10*d0/9)^10],
 TrueQ[twoLower>173/200]
 };
Print["guard_vector=",guards];If[!And@@guards,Exit[1]];Print["all_guards=True"];
