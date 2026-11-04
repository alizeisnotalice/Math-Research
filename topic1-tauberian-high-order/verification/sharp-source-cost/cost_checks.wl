Print["VERSION ",$Version," PID ",$ProcessID];
(* By symmetry of x,y and of +/- it is enough to consider 0<=c and x<=y. *)
fp=(x+c)(y+c)(x+y+2c);fm=(x-c)(y-c)(x+y-2c);
Print["inner branch excess=",InputForm[Factor[(fp+fm)/2-x y(x+y)]]];
Print["inner inequality violation=",InputForm[Reduce[0<=c<=x<=y && (fp+fm)/2 < x y(x+y),{x,y,c},Reals]]];
outer=Expand[((2x+t)(2x+s+t)(4x+s+2t))/2-x(x+s)(2x+s)];Print["outer branch positive polynomial=",InputForm[outer]];
Print["outer inequality violation=",InputForm[Reduce[x>=0&&s>=0&&t>=0&&outer<0,{x,s,t},Reals]]];
Print["left clip violation=",InputForm[Reduce[x>=0&&y>=0&&z< -Max[x,y]&&(x+z>0||y+z>0),{x,y,z},Reals]]];
Print["right clip violation=",InputForm[Reduce[x>=0&&y>=0&&z> Max[x,y]&&(x-z>0||y-z>0),{x,y,z},Reals]]];
(* Two targets, independently ordered P and Q. *)
i1=Integrate[c1 p1(q1-y)+c2 p2(q2-y),{y,0,q1}]+Integrate[c2 p2(q2-y),{y,q1,q2}];
i2=Integrate[c1 q1(p1-y)+c2 q2(p2-y),{y,0,p1}]+Integrate[c2 q2(p2-y),{y,p1,p2}];
i2reverse=Integrate[c1 q1(p1-y)+c2 q2(p2-y),{y,0,p2}]+Integrate[c1 q1(p1-y),{y,p2,p1}];
target=(c1 p1 q1(p1+q1)+c2 p2 q2(p2+q2))/2;
Print["target integrals, both P orderings=",InputForm[FullSimplify[{i1+i2==target,i1+i2reverse==target}]]];
Print["one target mass=",InputForm[FullSimplify[Integrate[p(q-y),{y,0,q}]+Integrate[q(p-y),{y,0,p}]==p q(p+q)/2]]];
(* All powers base-positive, so no branch-changing PowerExpand. *)
ass=Element[n,Integers]&&n>=2&&v>0&&r>0;
norm=2(v/v^n)(v r)^(n-2)(v r)^(n-2)/(n v^n r^(n-1));
Print["two endpoints normalization=",InputForm[FullSimplify[norm==2 r^(n-3)/(n v^3),ass]]];
Print["after half-integral normalization=",InputForm[FullSimplify[norm/2==r^(n-3)/(n v^3),ass]]];
H=Min[q,p q/v^2];cost=(v/v^n)(v r)^(n-2)(p r)(H r);rhs=r^n p q Min[1/v,p/v^3];
Print["large-angle normalization=",InputForm[FullSimplify[cost==rhs,ass&&0<q<=p]]];
Print["H inequalities violations=",InputForm[Reduce[v>0&&p>0&&q>0&&(Min[q,p q/v^2]>q||v^2 Min[q,p q/v^2]>p q),{v,p,q},Reals]]];
Print["large branch violation=",InputForm[Reduce[0<v<1&&c0>0&&C0>0&&0<q<=p<=C0 v&&p>=v^2&&X>=c0 q&&p q>(Max[C0,1]/c0)v X,{v,c0,C0,p,q,X},Reals]]];
Print["small branch violation=",InputForm[Reduce[0<v<1&&c0>0&&C0>0&&0<q<=p<=C0 v&&p<=v^2&&X>=c0 p q/v^2&&p q>(Max[C0,1]/c0)v X,{v,c0,C0,p,q,X},Reals]]];
