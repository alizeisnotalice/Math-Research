ClearAll[ax,bx,cx,dx,l];
Print["Mathematica: ",$Version];
Print["identity residual: ",Expand[(ax-dx)+(cx-bx)-(ax-bx)-(cx-dx)]];
Print["domain: l>0, ax-dx>=l, cx-bx>=l; alleged ax-bx<l and cx-dx<l"];
Print["violation set: ",Reduce[l>0&&ax-dx>=l&&cx-bx>=l&&ax-bx<l&&cx-dx<l,{ax,bx,cx,dx,l},Reals]];
