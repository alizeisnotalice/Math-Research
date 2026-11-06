$left42LibraryMode=True;Get[FileNameJoin[{DirectoryName[$InputFileName],"left_cost_probe42.wl"}]];$left42LibraryMode=False;
argsR42=$ScriptCommandLine;outR42=If[Length[argsR42]>=2,argsR42[[2]],"/tmp/left_raw42.json"];verifyPathR42=If[Length[argsR42]>=3,argsR42[[3]],FileNameJoin[{DirectoryName[DirectoryName[$InputFileName]],"rounds","r042","left_cost_verify42.json"}]];
rawRows42={};Do[ll=N[(ff/(aa mu42(mu42-1)/8 tt))^(1/3),60];x=case42[tt,ll,aa,64,1/16,False];AppendTo[rawRows42,x];Print[KeyTake[x,{"t","alpha","phase_fraction_alpha_J_over_D","residual_scaled_lambda2","min_action","seconds"}]],{tt,{512,8192}},{aa,{101/100,11/10}},{ff,{1/4,1}}];
verify42=Import[verifyPathR42,"RawJSON"];
Do[x=case42[row["t"],SetPrecision[row["lambda"],50],Rationalize[row["alpha"],10^-10],64,1/16,False];AppendTo[rawRows42,x];Print[KeyTake[x,{"t","alpha","phase_fraction_alpha_J_over_D","residual_scaled_lambda2","min_action","seconds"}]],{row,Select[verify42["near_contact"],#["t"]==8192&]}];
Export[outR42,<|"wolfram"->$Version,"status"->"unclipped sufficient-condition diagnostics; floating not certified","rows"->rawRows42|>,"RawJSON"];
