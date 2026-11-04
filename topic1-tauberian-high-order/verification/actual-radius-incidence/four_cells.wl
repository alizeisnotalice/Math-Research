Print["VERSION ",$Version," PID ",$ProcessID];
ass=ax<=bx && cx<=dx && ay<=cy && by<=dy && dx<ax && dy<ay;
Print["opposite-diagonal overlap violation=",InputForm[Reduce[ass && (cx>=bx || by>=cy),{ax,bx,cx,dx,ay,by,cy,dy},Reals]]];
Print["strict separation direct=",InputForm[Resolve[ForAll[{ax,bx,cx,dx,ay,by,cy,dy},Implies[ass,cx<bx && by<cy]],Reals]]];
box=0<eps<1/10 && eps<=xx<=1+eps && -eps<=yy<=1-eps && Max[Abs[xx],Abs[yy]]>=1 && Max[Abs[xx-1],Abs[yy-1]]>=1;
Print["slot violation=",InputForm[Reduce[box && !(1<=xx<=1+eps && -eps<=yy<=0),{eps,xx,yy},Reals]]];
Print["slot implication=",InputForm[Resolve[ForAll[{eps,xx,yy},Implies[box,1<=xx<=1+eps && -eps<=yy<=0]],Reals]]];
