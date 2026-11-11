mu = x /. FindRoot[x == 1 + Exp[-x], {x, 13/10}, WorkingPrecision -> 80];
excdf[x_, a_, s_] := NIntegrate[Exp[-u] CDF[NormalDistribution[0,s],x-a u],{u,0,Infinity},WorkingPrecision->60,AccuracyGoal->35,PrecisionGoal->35];
hypo[x_,a_,b_,s_] := Which[s==0&&x<=0,0,s==0&&b==0,1-Exp[-x/a],s==0&&a==b,1-Exp[-x/a](1+x/a),s==0,1-(a Exp[-x/a]-b Exp[-x/b])/(a-b),b==0,excdf[x,a,s],a==b,NIntegrate[u Exp[-u] CDF[NormalDistribution[0,s],x-a u],{u,0,Infinity},WorkingPrecision->60,AccuracyGoal->35,PrecisionGoal->35],True,(a excdf[x,a,s]-b excdf[x,b,s])/(a-b)];
prob[r_,c_,cp_,s_] := Total[Flatten[Table[(If[i==0,1-1/(1+Exp[c]),1/(1+Exp[c])]) (If[j==0,1-1/(1+Exp[cp]),1/(1+Exp[cp])]) hypo[Chop[mu(1+r)-(If[i==1,c,0]+r If[j==1,cp,0]),10^-40],1,r,s],{i,0,1},{j,0,1}]]];
Print["mu=",N[mu,45]];
pts={{1,mu,mu,0},{1,mu,mu,1},{1/100,mu,mu,4},{1/100,mu,mu,0},{1/1000000,2,3,1},{0,mu,0,0},{1,0,0,0},{1,30,30,34},{1/1000,mu,mu,34},{7777737242868376/10^20,13553876772791646/10^16,11952418680764778/10^15,19385004174467444/10^16}};
Do[Print[Append[q,N[prob@@q,35]]],{q,pts}];
