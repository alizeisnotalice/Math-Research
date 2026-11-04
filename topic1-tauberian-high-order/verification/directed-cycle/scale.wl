ClearAll["Global`*"];
Print["VERSION ",$Version," PID ",$ProcessID];
aa=(2qq-1)/(qq-1);bb=qq aa;
checks={
 FullSimplify[bb==2qq+qq/(qq-1),qq>1],
 FullSimplify[aa(qq-1)==2qq-1,qq>1],
 FullSimplify[aa qq==aa(2qq-1)-2qq+1,qq>1],
 Reduce[qq>=3&&bb<=6,qq,Reals]===False,
 FullSimplify[(aa/.qq->2)==3&&(bb/.qq->2)==6],
 FullSimplify[(ee^qq/uu^qq)/(ee^qq/uu)==uu^(1-qq),ee>0&&uu>0&&qq>=2],
 FullSimplify[(ee^(2qq-1)/uu^(2qq))/(ee^qq/uu)==ee^(qq-1)/uu^(2qq-1),ee>0&&uu>0&&qq>=2]
};
Print["CHECKS ",InputForm[checks]];Print["FAILURES ",Count[checks,Except[True]]];
If[!And@@checks,Exit[1]];
