Print["Wolfram ", $Version, " PID ", $ProcessID];
base=DirectoryName[$InputFileName]; j=Import[FileNameJoin[{base,"candidate002.check2.json"}],"RawJSON"];
new=FromDigits[j["ratio_num"]]/FromDigits[j["ratio_den"]];
old=482167717377842991882283179477574801392641890014291050323/283174368139841919615520431313000857449005056000000000000;
m=65944170400277227073796193886781311/12500000000000000000000000000000000; r=6916138095063/1099511627776;
Print["exact new > old: ",new>old]; Print["root guard: ",r^2>30m/4];
Print["N=10^8 > old: ",new (1-8/10^8)^2>old];
Print["N=10^12 > floor: ",new (1-8/10^12)^2>17027531463/10^10];
Print["exact gain: ",new-old];
