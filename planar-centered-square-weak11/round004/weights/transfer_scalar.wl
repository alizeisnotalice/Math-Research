Print["PID ",$ProcessID," version ",$Version];
z=Import[FileNameJoin[{DirectoryName[$InputFileName],"transfer_34_27.check2.json"}],"RawJSON"];r=FromDigits[z["ratio_num"]]/FromDigits[z["ratio_den"]];old=4629593010326846247058837068519417805352716109376907265077/2718886251002898428308996140604808231510448537600000000000;
m=65941670400277227073796193886781311/12500000000000000000000000000000000;rho=1729001748865/274877906944;
Print["strict improvement ",r>old];Print["sqrt bound ",rho^2>30m/4];Print["finite lower > 1.7027535096 ",r(1-8/10^12)^2>17027535096/10^10];Print["gain exact ",r-old];
