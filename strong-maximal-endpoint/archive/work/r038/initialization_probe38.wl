(* Actual finite Bellman initialization diagnostic; not a certified bound. *)
repo38="/Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal/strong-maximal-endpoint";
$nonlinear36LibraryMode=True;Get[FileNameJoin[{repo38,"code","nonlinear_probe36.wl"}]];$nonlinear36LibraryMode=False;
lift38[l_,v_,u_]:=Module[{aa=l-v,q,ee},If[aa<=0,Return[l]];If[u<=0,Return[Min[l,v]]];If[v<=0,Return[0.]];q=Log[v/aa]+u l/(v aa);If[q>=0,l/(1+Exp[-q]),ee=Exp[q];l ee/(1+ee)]];
run38[hh_]:=Module[{nn=100000000,mt=128,grid,pr,raw,st,ww,ex,den,ids,rows={},tm,core,idx,kind,bet,ii},
 grid=N[Range[0,Ceiling[(mt mu36+8Sqrt[mt mu36])/hh]]hh];pr=Developer`ToPackedArray[Exp[-grid]];
 For[tm=1,tm<=mt,tm++,raw=compiledStage[pr,nn,tm,hh,mu36];pr=Developer`ToPackedArray[Take[raw,Length[grid]]];
 If[MemberQ[{16,64,128},tm],st=fast36[nn,tm,#]&/@grid;
 Do[ww=If[kind==="logistic",MapThread[lift38,{st[[All,4]],st[[All,1]],1.1st[[All,2]]+bet weight36[tm]st[[All,6]]}],MapThread[Min,{st[[All,4]],st[[All,1]]Exp[(1.1st[[All,2]]+bet weight36[tm]st[[All,6]])/st[[All,1]]]}]];
 ex=Map[Max[0.,#]&,pr-ww];den=st[[All,4]]-ww;ids=Flatten@Position[(#>10^-10&)/@den,True];
 Do[core=Flatten@Position[(#<=aa Sqrt[tm Log[tm]]&)/@Abs[grid-tm mu36],True];idx=First@Ordering[ex[[core]],-1];
 AppendTo[rows,<|"N"->nn,"t"->tm,"h"->hh,"kind"->kind,"beta"->bet,"A"->aa,"maximum_positive_error"->Max[ex[[core]]],"error_area_in_tube"->hh Total[ex[[core]]],"maximum_fraction_of_coarse_gap"->If[Length[Intersection[ids,core]]>0,Max[ex[[Intersection[ids,core]]]/den[[Intersection[ids,core]]]],0.],"at_r"->grid[[core[[idx]]]],"status"->"PWL adaptive Bellman and machine mixtures; diagnostic only"|>],{aa,{1.,2.,3.}}],{kind,{"exponential","logistic"}},{bet,{0.,.01}}]]];rows];
args38=$ScriptCommandLine;out38=If[Length[args38]>=2,args38[[2]],"/tmp/initialization_probe38.json"];
rows38=Flatten[run38/@{.05,.025}];Export[out38,<|"wolfram"->$Version,"rows"->rows38|>,"RawJSON"];Print[<|"rows"->Length[rows38],"output"->out38|>];
