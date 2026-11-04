(* Exact all-parameter scalar certificates only; face geometry is a separate proof. *)
Clear[A,b,t,z,q,zstar,bound,expr,certificate,flags,failures];
Print["VERSION: ",$Version];
zstar=2A t^5/(3b);bound=4A^3 t^18/(27b^2);
expr=A t^8 z^2-b t^3 z^3;certificate=b t^3(z-zstar)^2(z+zstar/2);
flags={FullSimplify[bound-expr==certificate],
 FullSimplify[certificate>=0,A>0&&b>0&&t>0&&z>=0],
 FullSimplify[(expr/.z->zstar)==bound,A>0&&b>0&&t>0],
 FullSimplify[q/(1-q/(2q-1))==2q+q/(q-1),q>=2],
 FullSimplify[2q+q/(q-1)-6==(q-2)(2q-3)/(q-1),q>=2],
 Reduce[q>=2&&2q+q/(q-1)<6,q,Reals]===False};
Print["positive substitutions v=t^3, epsilon=z^3: ",ToString[flags,InputForm]];
Print["maximum: 4 A^3 v^6/(27 b^2), attained at epsilon=8 A^3 v^5/(27 b^3)"];
failures=Count[flags,Except[True]];Print["FAILURES: ",failures];Exit[If[failures==0,0,1]];
