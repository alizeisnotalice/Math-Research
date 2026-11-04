(* Algebraic and semialgebraic checks; the all-set proof is in the manuscript. *)
ClearAll[x,y,l,f,g,s,n,b,d,p,q,z];
convexPositive=Resolve[ForAll[{x,y,l},0<=l<=1,
 Max[l x+(1-l)y,0]<=l Max[x,0]+(1-l)Max[y,0]],Reals];
shiftMin=Resolve[ForAll[{f,g,s},Min[f+s,g+s]==Min[f,g]+s],Reals];
normalization=Resolve[ForAll[{s,n,b,d},n>0&&0<=s<=n&&0<b<1&&d>=0,
 n+d-b(s+d)>0],Reals];
flatIdentity=Resolve[ForAll[{d,p,q},d>=0&&p>=0&&q>=0,
 Exists[z,Max[Abs[d-z],p]+Max[Abs[-d-z],q]==Max[2d,p+q]]],Reals];
flatLower=Resolve[ForAll[{d,p,q,z},d>=0&&p>=0&&q>=0,
 Max[Abs[d-z],p]+Max[Abs[-d-z],q]>=Max[2d,p+q]],Reals];
scaleResidual=Expand[(n-b)s+b s-n s];
betaResiduals=Table[FullSimplify[p Beta[p,j+1]-Factorial[j]/Product[p+k,{k,1,j}],Assumptions->p>0],{j,1,8}];
If[betaResiduals=!=ConstantArray[0,8],Exit[1]];
result=<|"kind"->"exact-semialgebraic-identities","version"->$Version,
 "positive_part_convexity"->convexPositive,"minimum_translation"->shiftMin,
 "fill_normalization_derivative_positive"->normalization,
 "two_dimensional_flat_fibre_attainment"->flatIdentity,
 "two_dimensional_flat_fibre_lower_bound"->flatLower,
 "uniform_scale_algebraic_residual"->scaleResidual,"moment_Beta_residuals"->betaResiduals|>;
Print[result];
If[!And[convexPositive,shiftMin,normalization,flatIdentity,flatLower]||scaleResidual=!=0,Exit[1]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"..","results","global_leading_wolfram.json"}],result,"RawJSON"];
