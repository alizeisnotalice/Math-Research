ClearAll["Global`*"];
Print["Version ", $Version];
check[label_, truth_] := If[TrueQ[truth], Print[label, " = True"], Print["FAILED ", label, " = ", truth]; Exit[1]];

mat3={{1+ga1 u^6,ss,aa u^3},{hh u^6,1+ga2 u^6,bb u^3},{cc u^3,dd u^3,1+ga3 u^6}};
check["symbolic shear and genuine triple term", Expand[Coefficient[Det[mat3],u,6]-(ga1+ga2+ga3-ss hh-aa cc-bb((1-ss)dd+ss(dd-cc))) ]==0];
check["quantitative stability factor", Factor[1-y(3-y)^2/4-(y-1)^2(4-y)/4]==0];
check["stability cutoff 2/7", Reduce[x>=0 && y>=0 && x y(3-x-y)>1-1/49 && Abs[y-1]>=2/7,{x,y},Reals]===False];
check["true face convex weights", Expand[(1-ss)(1+ga2 u^6)+ss(1+(ga2-hh)u^6)-(1+(ga2-ss hh)u^6)]==0];

Do[
  ss=2/5; hh=-3/7;
  gg=Table[(i^2-3i)/11,{i,n}];
  bm=Table[If[i==j || (i+j==3),0,Mod[7i-3j+i j,13]/9-2/3],{i,n},{j,n}];
  mm=IdentityMatrix[n]+u^3 bm+u^6 DiagonalMatrix[gg]; mm[[1,2]]=ss; mm[[2,1]]=hh u^6;
  val=Sum[bm[[i,j]]bm[[j,i]],{i,3,n},{j,i+1,n}]+Sum[bm[[1,j]]bm[[j,1]]+bm[[2,j]]((1-ss)bm[[j,2]]+ss(bm[[j,2]]-bm[[j,1]])),{j,3,n}];
  check[StringJoin["dense rational n=",ToString[n]], Coefficient[Det[mm],u,6] == Total[gg]-ss hh-val],
  {n,3,6}];
Print["OK"];
