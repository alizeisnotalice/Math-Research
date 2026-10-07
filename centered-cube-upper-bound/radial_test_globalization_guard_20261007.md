# Radial-test globalization：三轮有理公式守卫终态

2026-10-07；gated_radial_energy_audit。使用 L03 的符号/残差/精确数据及失败保留流程。只新增本 prefix；解析终稿 radial_test_globalization_20261007.md 保持 SHA 56d36f2bb2908ab974340c706be749a80c3758812a2cf56e22cd6b43f8cfb790，未改动。

## 1. 冻结配置与真实输入范围

原 [registration](radial_test_globalization_guard_20261007_registration.json) 在执行前冻结。三轮维数分别 (2,16)、(4,64)、(8,256)；每维取 H=(1,2,4)(n²+2)、来源振幅 ε0=1/8,1/4，共36个完整来源/receiver 参数组。原连续尺度 [1,2]，τ=1/2；新候选仅改为有理 L=1+1/n。真实完整来源仍是
\[
 f_H=(1-\epsilon_0\|y\|_2^2/(nH^2))1_{[-H,H]^n},
 \quad W=(1-\epsilon_0/3)(2H)^n.
\]
receiver E0=[−H+1,H−1]^n 始终只是 full E={M>τ} 子集，|E0|=(2H−2)^n。代码保留完整 W 与 E0 的体积比；没有以 E0 当作新来源，也没有免费背景。

原稿的内盒平均二次式与唯一真赢家1对任何 L∈(1,2]均适用，因此有理候选是同一真实来源反例的参数变体。nlog(1+1/n)≥n/(n+1)≥2/3，log(1+1/n)≤1/n 都是解析界，不声称用 float 认证 log。

精确计算 uniform radial expectation
\[
 J_n=\frac{nL/(n+1)-1+L^{-n}/(n+1)}{L-1},
\]
并独立用积分反导数与展开 (1+(L−1)t)^n 后逐 monomial 正积分交叉核。解析上也可直接看前半区间：s/L≤1−1/[2(n+1)]，故 (s/L)^n≤e^{-n/[2(n+1)]}≤e^{-1/3}≤3/4，从而 Jn≥1/8。注册守卫只要求更保守 (3/4)Jn≥3/40，不从有限6个维数断言全部 n。

## 2. 首次失败与 serialization-only 修复

第一次执行实际 tool chunk 971075、exit 1：完整来源质量 W 的非认证 float(value) 诊断转换触发 OverflowError。没有 Fraction 断言失败报告，没有成功终态结果。保留 [failure](radial_test_globalization_guard_20261007_failure.json) 和原执行 [代码快照](radial_test_globalization_guard_20261007_failed_version.py)。工具日志未记录失败前完整 check 数，所以字段为 null；不制造时间戳或累计失败检查。

随后冻结 [amendment](radial_test_globalization_guard_20261007_registration_amendment.json)。数学输入、所有数学断言、证明均不变，只把诊断 float 溢出写 null，另加原 registration hash 完整性检查。修复代码执行一次并成功；成功计数不累加首失败尝试。

## 3. 成功终态

[最终脚本](radial_test_globalization_guard_20261007.py)、[results](radial_test_globalization_guard_20261007_results.json)、[receipt](radial_test_globalization_guard_20261007_receipt.json)。

**1020 个精确检查 PASS：三轮各338，加6个登记/代码/证明完整性检查；36参数组；0旧 oracle/probe 调用；无随机种子。** 修复执行实测耗时0.009982374991523102秒。三轮从头运行当前冻结算术一次，不把失败 run 已经过的参数当独立追加收据。

逐参数检查了来源完整质量、inner receiver体积、盒内 f≥3/4、inner response>τ与完整幅度上限、平均严格随R下降、R1 versus L/2严格赢家、原后验 gap 上界、ρ²=(H/(H−1))^n≤4、实际振幅径向差≥3/40及整体 all-v 的解析 scalar 上包≤36/n、必要转换常数下包≥n/480。ρ≤2用平方的有理 inequality核，不对 sqrt 浮点比较。

以下小数仅作可读诊断；results 同时保存 Fraction 或 numerator/denominator hash与向外 dyadic128上下端点。

|n|exact公式 Jn 的小数诊断|保守局部差 (3/4)Jn|
|---:|---:|---:|
|2|0.2962962963|0.2222222222|
|4|0.3276800000|0.2457600000|
|8|0.3464394161|0.2598295621|
|16|0.3567861947|0.2675896461|
|64|0.3650313185|0.2737734889|
|256|0.3671625598|0.2753719198|

全部组的最小局部下包为2/9；n倍 global mean上包的小数诊断最大约34.6685956790<36。它们是显式完整 L1 反例的积分公式算术证书，不是对任意全局 v 的离散采样。全 v 不等式由终稿自含 L2/Plancherel 证明承担；代码未枚举 source tests、receiver网格或全 M 边界 oracle。

## 4. Hash 与结论边界

- 最终 script SHA：348efc5207753f1cb2e8350d338024042c55c8a57a710fad25f8ff721cef6b60。
- results SHA：e4af84772a57f1200a229f5b86a84e466f0b0e3f2edff73cd5fbced5f5758bc8。
- receipt SHA：367e576bf30f1ace1a45d02786bdaebf58fb363f4023bd93647c6064ac5bc749。
- 原 registration SHA：b06107cb1220bc8882e471ec0521165f1c767912fdb9ca8804f977dfe699fc1c。
- amendment SHA：28a026b193d9054de3d29147dbf4ec6ec76164b1b075793c144a07eb5c6d9f37。
- 原失败代码 SHA：73f3616a7ebb15a03a772615d46b1e04ffe859e280e1aab1d453905f3acd202d。
- failure receipt SHA：259bc16bf7b17af180e2705f239b8929278d901e20faf91ae244a6a7c7479faa。

这些是原连续中心 cube、完整 L1、唯一 inner truewinner 的解析反例构件验证。它们不认证输入 near-max，不认证原 FIRST/CPGP/LCA/history，也不是新 weak-type 反例；仅支持终稿已证明的 normalized local-test/global-test 量词损失。未归一化物理 W1、专门近极值几何转换和一般 sqrt(n)主目标都不因此被否定。

