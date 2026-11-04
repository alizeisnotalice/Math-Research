Print["VERSION ",$Version];
ass=Element[n,Integers]&&n>=2&&c>=0&&Element[z,Reals];
Print["identity=",InputForm[FullSimplify[S+n z-2n(c+Max[z,0])==S-2n c-n Abs[z],ass]]];
Print["bound violation set=",InputForm[Reduce[n>=2&&c>=0&&lam>=c+Max[z,0]&&S+n z-2n lam>S-2n c-n Abs[z],{n,c,lam,z,S},Reals]]];
Print["unique z optimizer=",InputForm[Resolve[ForAll[{n,z},Implies[n>=2&&z!=0,-n Abs[z]<0]],Reals]]];
Print["symmetric_matrix_optimum=",InputForm[FullSimplify[Binomial[n,2](8n/3)^2(1-2(8n/3)/(8n))/64==Binomial[n,2] n^2/27,Element[n,Integers]&&n>=2]]];
