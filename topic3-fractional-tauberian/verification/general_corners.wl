(* Algebraic certificates; the measure and uniform-remainder arguments are in the paper. *)
Clear[x,u,z,v,a,b,c,d];
s=x-u;q=1-x-u;
checks={
 Expand[2x-(u^2+u+x(1-x))-(2(x^2-u^2)+s q)]===0,
 Expand[2(1-x)-(u^2+u+x(1-x))-(2((1-x)^2-u^2)+s q)]===0,
 Resolve[ForAll[{x,u},0<=x<=1/2&&-x<=u<=x,0<=s<=1&&0<=q<=1&&x^2-u^2>=0],Reals],
 Resolve[ForAll[{x,u},1/2<=x<=1&&-(1-x)<=u<=1-x,0<=s<=1&&0<=q<=1&&(1-x)^2-u^2>=0],Reals],
 Resolve[ForAll[{x,u,v},0<=x<=1/2&&-x<=u<=x&&-x<=v<=x,(x-u)(2-(1-x-u)-(1-x-v))>=0],Reals],
 Resolve[ForAll[{x,u,v},1/2<=x<=1&&-(1-x)<=u<=1-x&&-(1-x)<=v<=1-x,(1-x-u)(2-(x-u)-(x-v))>=0],Reals],
 (* Four-region inclusion-exclusion: a,b,c,d = lower-left, upper-left, lower-right, upper-right masses. *)
 Expand[(a+b+c+d)-(b+d)-(c+d)-d-(a-2d)]===0,
 4 Integrate[2x,{x,0,1/2}]+4 Integrate[2(1-x),{x,1/2,1}]===2
};
If[!And@@checks,Print[checks];Exit[1]];
out=<|"kernel"->$Version,"checks"->checks,"scope"->"exact identities and universal polynomial inequalities; geometric transfer is proved separately"|>;
Export[FileNameJoin[{DirectoryName[$InputFileName],"..","results","general_corners_wolfram.json"}],out,"RawJSON"];
Print[out];
