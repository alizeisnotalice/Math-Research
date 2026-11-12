k=1/10;b=1/50;a=1/25;s=Sqrt[1+2/k];v=s k+Log[(s+1)/(s-1)];mu=1+ProductLog[1/E];checks={b+v<mu,a<k-b,a>b};Print[{checks,N[mu-b-v,18]}];If[!And@@checks,Exit[1]];
