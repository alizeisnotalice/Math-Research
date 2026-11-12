area=Integrate[(z-s a)/t-(z-s x)/t,{x,a,a+h}];
tests={FullSimplify[area==s h^2/(2t),t>0],FullSimplify[(area/.{s->1/Sqrt[2],t->1/Sqrt[2],h->1/2})==1/8]};Print[tests];If[!TrueQ[And@@tests],Exit[1]];
