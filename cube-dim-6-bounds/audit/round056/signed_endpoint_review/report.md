# F056：正 B 混号端点积分独立审计

**PASS，限原body32的同31个低端管、全部正时间。** 对rho=1/64、q=6/5、u=mean_B x∈[0,1/4]，独立有理证书确认

    ||G||_(L^(6/5)(union_A E_A))
      < 13404268290817/1000000000000000000
      = 0.000013404268290817。

准确独立upper约1.3404268290815612e−5，严格小于旧054上预算公式完整区间的lower endpoint。新/旧公式比值upper约0.6511327199337893，即同一局部预算公式约下降34.89%；旧公式lower不是实际G的下界，预算下降也不是C6改进。

[proof_review.md](proof_review.md) 给出逐步独立链条。对正B坐标，两组符号和非负，Cauchy–Schwarz给mixed plane距离≥sqrt(h)u；unsigned projection dominance同时给≥r，故可用max(r,sqrt(h)u)。先在真实正域替换非负integrand，再放大normal ball，避免把positivity不等式用于负坐标。真实Jacobian sqrt(km)、maxcone fraction1/k以及u=0 clipping方向均继承F054且使用正确。

mean积分在u=r/sqrt(h)处分段，随后准确radial积分得到

    R_h=J_h/J_0=(84/19)alpha_h−(7/2)alpha_h^(6/5)，
    alpha_h=1/(16sqrt(h))。

独立Fraction核系数，解析证明0<R_h<1。全负等价类按少数负号计数，m5至m1分别(5,10)、(4,3)、(3)、(1)、空；每个plane的两个有序项全部保留。Minkowski后正union仍用(1,5,10,10,5)，允许任意重叠。

| k | m | incident factor（仅展示） | 新norm upper（仅展示） |
|---|---|---|---|
|1|5|0.228242862258|6.071665726630e−6|
|2|4|0.288994630678|1.153274995069e−6|
|3|3|0.404917391977|4.948944343079e−7|
|4|2|0.603278261318|3.126179030644e−7|
|5|1|1，完整沿用054|完整沿用054|

五类background完整区间、五类原始singular/background/total来源区间，与作者054逐字相等。k=5三完整区间和已保存的加权q幂对象也逐字相等；没有重求k5、旧M定义积分或旧054完整预算。源输入不确定性保持两端，未压为midpoint。

独立实现先写[theory_before_results.md](theory_before_results.md)，再打开056算术结果，只使用纯Fraction、整数幂和整数根。登记10配置，实际一次10配置PASS；19个root证书逐项满足准确幂夹逼。所有对应作者区间均落在独立外包内，各规范有理上端比较及new upper<old formula lower均是准确比较，见[results.json](results.json)和[verify_signed.py](verify_signed.py)。

唯一科学运行PID **60298**，tool chunk **53815e**，exit **0**，耗时 **0.009299333秒**，于登记后 **139.142091秒**结束，处于900秒上限内。没有科学失败或重试。登记前首次构造来源相对路径时误用clone根，0配置即FileNotFoundError；修正project根后正常登记，未改源文件或公式。未调用作者代码、Arb、旧积分、完整pair枚举、kernel/native、grid、MATLAB或subagent，未stage/commit/push。

21来源pin执行前后及封存时全部一致。作者056 manifest为9305a1c06d472f5a712f4c31f580be5b3e653dd7edb01b751cef0e3413f48e71；其gate_release已绑定最终F054 PASS和主代理release，原超时等待批0科学配置的历史记录保留。本审计无需重跑任何旧证书。

高端、中段长时间尾、管外和global g/I仍未支付。无全域证书、无新C6界。主代理保留最终证明及后续预算组合责任。
