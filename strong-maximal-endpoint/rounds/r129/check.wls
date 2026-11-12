den=1+d+y;r=2(2+y)/(q den)-2;s=2(1+1/den)/q-2;
ass=0<d<=1/10&&0<y<d&&0<q<=1-4d;
tests={FullSimplify[r>=2/q-1&&1/den>q,ass],FullSimplify[s-r==2d/(q den),ass],FullSimplify[Integrate[4d/(q den),{y,0,d},Assumptions->d>0&&q>0]==4d/q Log[(1+2d)/(1+d)],d>0&&q>0]};
Print[tests];If[!And@@tests,Exit[1]];
