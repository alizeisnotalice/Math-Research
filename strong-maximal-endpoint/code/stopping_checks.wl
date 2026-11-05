root=DirectoryName[DirectoryName[$InputFileName]];
Clear[n,k,mu,cv,v]; ell=Log[n];
eps=32/Sqrt[mu n]; delta2=cv n^(-1/4)ell^(1/4); delta=Sqrt[delta2]; x=2Sqrt[ell];
errors={x^3 eps,x^2 delta2,(1+x)(eps Log[1/eps]+delta)};
limits=Limit[#,n->Infinity,Assumptions->{mu>0,cv>0}]&/@errors;
ratio=FullSimplify[2^(k+1)/(32^(k-2)/2),Assumptions->Element[k,Integers]&&k>=3];
threshold=2Sqrt[mu]+8Sqrt[v];
chernoff=FullSimplify[(threshold-2Sqrt[mu])^2/(8v),Assumptions->{mu>0,v>0}];
result=<|"wolfram_version"->$Version,
 "FGL_error_term_limits"->(ToString[#,InputForm]&/@limits),
 "Bernstein_factor_ratio_k_ge_3"->ToString[ratio,InputForm],
 "stop_success_negative_log_n_coefficient"->ToString[chernoff,InputForm]|>;
Export[FileNameJoin[{root,"rounds","r013","symbolic_checks.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
