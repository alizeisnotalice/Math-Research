Print["VERSION ",$Version];k=16/3-64u/3;ll=128/27-1280u/27+1024u^2/9;w=u/(1+2u)+ll u^6+u^8;rho=1-2w;e=k u^3;a=4u^2;det=1-e^2;mass=a rho^2;m0=(2u/(1+2u))^2;
h=(e^2/a)/(Sqrt[rho^2+e^2/a]+rho);A=2(e h-a/e (rho h^2+h^3/3));Dlt=4w^2-mass;G=(mass-det m0)/((det-mass)(1-m0)a(2u)^6);
Print["budget/u9=",InputForm[FullSimplify[Normal[Series[(Dlt-2A)/u^9,{u,0,1}]],u>0]]];
Print["G=",InputForm[FullSimplify[Normal[Series[G,{u,0,1}]],u>0]]];
