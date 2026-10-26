# A056：正象限混合符号端部积分预算

对原 body32 模型、rho=1/64、q=6/5、mean_B∈[0,1/4]，在054的同31个低端管及全部正时间上，得到严格上界

    ||G||_(L^(6/5)(union_A E_A))
      < 13404268290817/1000000000000000000
      = 0.000013404268290817.

该有理数经纯Fraction幂比较严格小于054原上预算公式的区间下端。新旧公式比值约0.651132719934；这是同域预算改善，没有推出新C6界。054已支付的所有background项保留，k=5的完整singular/background/total区间逐字复用，没有重算。

## 证明及范围

[theorem.md](theorem.md) 的数学内容与旧登记完全一致。关键是在实际非负B坐标域，具有h个少数符号的incident plane满足

    dist(x,P_sigma) >= max(r,sqrt(h)u),

其中r=dist(x,U_A)，u=mean_B(x)。第一项来自非负坐标的unsigned projection dominance；第二项来自非负符号组和的Cauchy–Schwarz界。先在真实域替换integrand，再删去box/positivity clipping作非负积分放大；不在负B坐标上使用正象限引理。

真实face Jacobian sqrt(km)、四维法向radial measure与maxcone fraction 1/k沿用已接受054证明。精确积分给

    alpha_h = 1/(16 sqrt(h)),
    R_h = (84/19)alpha_h - (7/2)alpha_h^(6/5), h=1,2.

全混合符号projective multiplicities为m=5,4,3,2,1时的(5,10)、(4,3)、(3,0)、(1,0)、(0,0)，加一个unsigned plane后恰2^(m-1)，每plane对应两个ordered terms。因此无删项、无全方向枚举、无新disjointness假设。Minkowski以后仍按相同权重(1,5,10,10,5)作正union。

|k|新incident norm系数（仅展示）|单管新总上预算（仅展示）|
|---|---:|---:|
|1|0.228242862258|6.071665726630e-6|
|2|0.288994630678|1.153274995069e-6|
|3|0.404917391977|4.948944343079e-7|
|4|0.603278261318|3.126179030644e-7|
|5|1，完整复用054|2.846907459672e-7|

真实Green可更小；这里严格比较的是两个有效的上预算公式。没有支付高端、middle长时间尾、指定管外、global g或全face I。条件目标4/48750仅作预算参照，新公式约占0.163364519794；它不是已闭合的全局条件。

## 续批门控与实际执行

旧registration.json和waiting_record.json原件保留。旧批绝对期限2026-10-25T22:50:48Z到达时仍处F054门控等待，0新scalar、0生产配置、0纯检查；没有静默延展旧登记。

continuation_registration.json于2026-10-26T08:52:52.329588Z登记新的固定9生产+1纯检查。准备开始于08:51:17Z；F054最终manifest和report及21个source pins实际逐项核对后，08:54:25.615120Z绑定主流程最终release。准备与门控重叠的时间单列，没有计入科学执行。

F054最终manifest SHA为805b4363b7f0e856a2c719d6def72f3a2e225154c7d5f45d4c331d976ddb123a，最终report SHA为1373a7ff9333b370b7d6e6b0ae81fbf6958b3808e25666093d3aae44d79e63b4。其source pins绑定本clone原054 manifest a9cc16dff037e55071510769e1bcb1d76d9010becfe01e00ca438a8653ffbd4a及完整源区间。gate_release.json同时绑定本续批登记和三个执行/证明文件。

科学时钟从2026-10-26T08:54:32.036581Z开始，硬期限为09:14:32.036581Z。唯一producer PID58976、tool chunk 5bac73、exit0，9/9配置用时0.003592084秒；唯一纯检查tool chunk 39178b、exit0，用时0.001707334秒，于08:54:38.362923Z完成。科学启动至纯检查完成共6.326341867秒，包含两次命令间隔。封存时producer已退出，无活句柄、失败或重试。

仅修改了evaluate.py/check_saved.py的续批登记、最终release、期限与来源记录接线；数学配置、参数及theorem均未改变。Python-flint0.9.0在768位产生向外dyadic端点；checker独立使用Fraction的五次、六次整数幂不等式认证所需方向，且不调用Arb。这里的“独立”指不同算术路径，不替代主流程对数学证明的最终责任。

旧M定义积分、旧054完整预算、k5新数值求值、kernel/native、网格、全方向枚举、MATLAB均为0；未stage/commit/push。用户建议文档只读作路线建议，没有把建议当新增定理。

## 复核入口与冻结数据

results.json保存两个R_h、k=1..4新预算、k5原区间、所有旧background完整区间以及正union。checked_results.json保存严格有理结论及10配置总数；run_receipt.json与execution_note.json给出实际工作量、时间及来源。manifest.json冻结必要文件SHA；旧登记与等待原件作为中断和续批证据一并保留。

原执行命令从cube-dim-6-bounds项目根目录运行：

```sh
python -B upper/round056/signed_endpoint_integration/evaluate.py --run
python -B upper/round056/signed_endpoint_integration/check_saved.py
```

本目录已冻结且脚本阻止覆盖生产收据。重放须复制本包与完整054依赖到新的项目布局，排除results.json、run_receipt.json、checked_results.json；保持所有源文件内容及哈希不变。gate_release.json中的审计绝对路径记录了本次原始来源位置，异机重放需要将相同哈希的F054文件映射到本机位置并显式生成新的重放release，不能篡改本次生产release或声称原收据仍绑定已改路径。复核冻结结论可直接检查保存的dyadic端点和Fraction幂比较，无需重积分旧M或重跑旧054。

下一步是主流程独立复核本正象限距离与积分替换证明，并与其余未付区域预算组合；本包自身不证明全域g、全face I或C6改进。
