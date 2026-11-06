base="/Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal/strong-maximal-endpoint";
s=Import[FileNameJoin[{base,"rounds","r030","receiver_seed.json"}],"RawJSON"];r=Import[FileNameJoin[{base,"rounds","r031","checks.json"}],"RawJSON"]["receiver_repair"];
w=r["uniform_scale_with_extra_margin"] s["weights"];av=r["final_N256_inner64"]["averages_by_group"];co=Join[{1},ConstantArray[2,Length[w]-1]];q=1+10^-5;
slack=co w(av-q);rows=MapThread[<|"group"->#1,"v"->#2,"weight"->#3,"average"->#4,"complementarity_slack"->#5|>&,{Range[Length[w]],s["v_nonnegative"],w,av,slack}];
js=r["final_N256_inner64"]["scaled_cost"];gap=Total[slack];result=<|"wolfram"->$Version,"status"->"floating continuous KKT diagnostic; no certificate","primal_scaled"->js,"dual_scaled"->js-gap,"relative_gap"->gap/js,"largest_active_slacks"->TakeLargestBy[Select[rows,#["weight"]>10^-10&],#["complementarity_slack"]&,8],"active_count"->Length[Select[w,#>10^-10&]]|>;
Export["/Users/zhengzhihao/Desktop/T/work/r032_seed_kkt.json",result,"RawJSON"];Print[result];
