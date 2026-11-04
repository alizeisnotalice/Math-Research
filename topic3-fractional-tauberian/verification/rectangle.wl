ClearAll[al,w,s,t,k];
s=(w+1/w)/2; xx=(1/al-1/w)/2; yy=(1/al-w)/2;
checks=<|"intersection1"->FullSimplify[xx+(2al/w-1)yy-(s-al)],"intersection2"->FullSimplify[(2al w-1)xx+yy-(s-al)],"area"->FullSimplify[2(s-al)(xx+yy)-2(s-al)(1/al-s)],"square_gap"->FullSimplify[(1/al-al)^2/2-2(1-al)^2/al-(1-al)^4/(2al^2)],"derivative"->ToString[FullSimplify[D[2(z-al)(1/al-z),z]],InputForm],"series"->ToString[Series[(1/(1-t)-(1-t))^2/2,{t,0,10}],InputForm],"version"->$Version|>;
Export[FileNameJoin[{DirectoryName[DirectoryName[$InputFileName]],"results","rectangle_wolfram.json"}],checks,"RawJSON"];Print[checks];
