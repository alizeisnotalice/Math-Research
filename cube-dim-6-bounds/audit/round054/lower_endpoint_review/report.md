# F054：低端管全时间预算独立审计

**PASS，限原body32的全部31个指定低端管。** 对rho=1/64、q=6/5、mean_B∈[0,1/4]，独立证明与准确复核确认

    ||G||_(L^(6/5)(union_A E_A))
      < 10293038485441/500000000000000000
      = 0.000020586076970882。

这是全正时间Green的局部区域预算。未支付高端、中段长时间尾、管外或global g，没有全域I或新C6界。

[proof_review.md](proof_review.md) 给出完整独立链条：端点incident planes恰2^(m−1)，共2^m个ordered terms、全属min(k,m)类；对正坐标其距离≥dist(x,U_A)。其余signed planes按partition是否相同完整分两类，给gap≥3/(4sqrt2)−sqrt(k)/64。norm/所有body line下界同时有效，五种background系数非负且全部992项已支付。

真实face Jacobian是sqrt(km)，正规径向积分的maxcone fraction是1/k。u=0额外clipping只能减少积分；本审计没有把upper体积称为实际体积。k=1及m=1均完整覆盖。Minkowski后用权重(1,5,10,10,5)作q次方正union，允许重叠。

独立实际批使用**纯Fraction及整数根**，与作者Arb实现不同。读取C037三个完整dyadic M区间，不重算其定义积分。共8个预登记配置：结构1、独立Machin pi夹逼1、固定k五类、union1。最初登记的7项准备计划保留，计算前把pi夹逼单独计为第8项，未改变任何参数或上限。

唯一运行 PID **55450**，tool chunk **0be28e**，exit **0**，耗时 **0.014782125秒**。36次整数根外包均有准确幂比较，五类结果各自低于作者规范有理上界，最终独立upper约2.058607697088027e−5，准确严格小于上述claim及(251/1000)*(4/48750)。小数只用于显示；证明数据见 [results.json](results.json)，最小独立代码见 [verify_budget.py](verify_budget.py)。没有失败、重试、作者writer、Arb、旧M积分、全pairs枚举、kernel/native、网格或MATLAB调用。

18项A054与直接来源哈希、另3项F035依赖解除哈希，执行前后共21项全部不变。A054 manifest为a9cc16dff037e55071510769e1bcb1d76d9010becfe01e00ca438a8653ffbd4a，theorem为8fabd20f4579aacb2b6e3876753f86c569f5ae3508910f664500bc90106dd739，results为bb22ef64c222ad57f5ffcdfb89f8541b963517c4f433ab13a75134f1a0f2e0b3。C037原结果的旧pending状态已由冻结audit_dependency和F037证明解除，本批只读核链，没有重跑旧积分或证书。

作者默认Python初次导入缺SymPy发生于0科学配置，随后原代码在已有venv一次成功；其记录完整。本审计直接用已有venv且仅Python标准库，无自身环境失败。

原独立科学批于2026-10-25T22:33:43.063493+00:00结束，即登记后169.528245秒，处于900秒上限内。本次接续仅核保存结果、收据哈希链及21个来源pin后补齐封存，新增科学配置为0；行政封存发生在中断恢复后，总登记至封存wall时间已超过900秒，不能表述为完整审计在900秒内封存。精确时间及文件哈希见manifest.json。主代理保留最终证明和预算组合责任。
