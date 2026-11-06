$left42LibraryMode=True;Get[FileNameJoin[{DirectoryName[$InputFileName],"left_cost_probe42.wl"}]];$left42LibraryMode=False;
wp42=80;mu42=N[1+ProductLog[1/E],110];pa42=1/(1+Exp[mu42]);ll=N[(1/((101/100)mu42(mu42-1)/8 512))^(1/3),100];x=case42[512,ll,101/100,64,1/16];
outP42=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/left_precision42.json"];Export[outP42,<|"wolfram"->$Version,"status"->"80-digit working precision repetition; floating not certified","row"->x|>,"RawJSON"];Print[KeyTake[x,{"residual_over_f","residual_scaled_lambda2","min_action","seconds"}]];
