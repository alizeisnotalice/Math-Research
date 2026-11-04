Print["VERSION ",$Version," PID ",$ProcessID];
f=(1-(1-eps)^q)/q;four=eps+(1-q)eps^2/2+(1-q)(2-q)eps^3/6+(1-q)(2-q)(3-q)eps^4/24;
Print["sigma fourth series ",InputForm[FullSimplify[Normal[Series[f,{eps,0,4}]]==four,0<q<=1]]];
Print["m coefficient ",InputForm[FullSimplify[((1-q)(2-q)(3-q)/24/.q->1/m)==(m-1)(2m-1)(3m-1)/(24m^3),m>=1]]];
Print["q zero extension ",InputForm[FullSimplify[Limit[f,q->0,Direction->"FromAbove"]==-Log[1-eps],0<eps<1]]];
Print["fifth derivative ",InputForm[FullSimplify[D[f,{eps,5}]==(1-q)(2-q)(3-q)(4-q)(1-eps)^(q-5),0<q<=1&&0<=eps<1]]];
Print["coefficient outside [0,24] ",InputForm[Reduce[0<=q<=1&&((1-q)(2-q)(3-q)(4-q)<0||(1-q)(2-q)(3-q)(4-q)>24),q,Reals]]];
Print["expansion polynomial nonpositive violation ",InputForm[Reduce[0<=q<=1 && ((1-q)<0||(1-q)(2-q)<0||(1-q)(2-q)(3-q)<0),q,Reals]]];
