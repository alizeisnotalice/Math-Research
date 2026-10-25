# B028：body32 全时间短势及任意正时间混合的障碍

**本轮证书通过，等待 F028 独立全群 H 交叉。** 在固定点
\(x_0=(1,1,1,9/10,27/50,0)\)，已证明对全部 \(T>0\)
\[
K_T(x_0)/T<\frac1{4100}(\ell/\pi)^6<\frac1{4100}<\frac1{4096},
\qquad \ell=3141592653/10^9.
\]
由离开真实 S2 的局部统一支配，对任意正时间概率混合 \(Q_\mu=\int K_T/T\,d\mu\)，这个失败都是 essential：
\(\operatorname{ess\,inf}_{[-1,1]^6}Q_\mu<1/4100\)。因此原 body32 的单时间平均和所有正时间概率混合，通过相应全 cube 支配产生的成本均严格大于 **64.0625**。这是方法家族障碍，不是 C6 的新下界或新上界；当前 C6 上界64不变。

## 严格证据

B024 已准确构造 \(p_T=N(T^2)/(64\pi^6D(T^2))\)，deg N=414、deg D=417、gcd=1、D在正轴严格正。临界多项式
\(P=N-64\ell^6(1/4100)D\) 有且仅有两个简单正根。完整性来自准确 Descartes 变号数2与两个不交的准确符号变化框，非近似根计数。

本轮统一使用
\(F(T)=T/4100-(\pi/\ell)^6K_T\)，导数是 \(-P/(64\ell^6D)\)。其导数符号依次为正、负、正，\(F(0+)=0\)，\(F(\infty)=+\infty\)。只需用锁定框上端的两个 K 值支付临界点下界；两次均严格通过。

|根|新 K 求值的准确时间 b|K_b 显示近似|严格付款的安全有理下界|
|---|---|---:|---|
|1|15536257/1073741824|1.436430694326986e−6|130791170401/62500000000000000|
|2|185807383/4294967296|1.0087094990298506e−5|232263757011/500000000000000000|

付款为 \(c a-(u/\ell)^6K_{b,upper}\)，c=1/4100、u=3141592654/10^9，a是 B024 原锁定下端。表中付款约2.0926587264e−6、4.6452751402e−7。所有比较使用完整 dyadic 端点和准确 Fraction；显示小数不参与认证。

`theorem.md` 给出整个正轴的单调论证，以及任意概率 μ 的完整连续性接口。特别 x0 的三个绝对坐标1、.9、.54保证 dist∞(x0,S2)≥1/20；三大参数事件和准确 basis Jacobian 给附近对所有 t 的有限 p_t 上包，故 q_T=K_T/T 对所有 T 共同有界。这个共同上包允许对无紧时间支撑、无时间矩假设的任意概率 μ 做 dominated convergence。没有从单点值直接跳到 essential 最小值。

## 算法、源锁定与实际运行

`registration.json` 在求值前锁定了两个时间、640位、20分钟、首失败停止及72份源/依赖哈希。登记 SHA：
`cbd69adabfa01e975bfb5a1e7a614b06d3fc1a438fb7a48067b3c39a51893d56`。

重要只读来源：

|来源|SHA256|
|---|---|
|B023 417项实主部表|8446b3fe84136c005d66501a64a7978e1c1adc452bc1993acd25545ddc8ff299|
|B024 最终 manifest|4ff73eaa3c4779494766794a7dad4fd014fbec225808c86eedb39ba1e22f490f|
|B024 完整 roots.json|08afc32a4a99982258f70be56fd3be40e933d3c09424602ff242811683d2c69b|
|C026 PF递推源码|5f0fe640bfe73fca40fcad55fcc55bdb6b066a5724664f1fa23abb62302013ad|
|C026 最终 manifest|10cbe75c05ec889fc3e8ac249172fe85c0240762ca2ea0e36cfc42f7ebe44b79|
|F023 独立完整符号审计报告|2851b03ea55e4ce3e95e7bbbae0b97daa9724a01f0de051aadc586255427544e|

本目录 `integrate.py` 独立实现匹配新时间的 PF 积分，没有 monkeypatch 或执行 C026 固定 T=1/25 的程序。417项、359个准确正有理平方根、最高阶3；全部正负系数进入 Arb640。运行使用 Python3.15.0、python-flint0.9.0、FLINT3.6.0，Python binary、Flint 扩展和 GMP/MPFR/FLINT 动态库均在登记中锁定。

唯一 PF 批 PID **40277**，2026-10-25 19:56:52 UTC 开始，0.005653708秒结束，exit0，tool chunk `fffb62`；无子进程、无遗留句柄。恰好2个新正时间 PF K 值，0 native、0 H、0 Gfull、0 p 新值、0额外时间、0优化、0网格、0新根隔离。没有重跑旧 Γ 的完整群计算，亦没有新增 Arb G 求值。

首次纯后核 `verify_saved.py` 因旧 B024 指数存为字符串而报 TypeError，未到核验结论、未修改结果。为保留预登记源码及历史，原文件原样保留；唯一修正是在新 `check.py` 中对两个 dyadic 指数加 `int()`。`checker_repair_receipt.json` 记录原/新源码哈希和失败 tool chunk `1b19dd`、exit1。没有重跑任何 PF 值或改变常数、时间、根框、端点。

修正后的纯后核 PID **40337**，0.473628417秒，exit0（chunk `542adf`），准确恢复全部417主部、N/D/P绑定、完整正根数与重数、π外包、旧G有限界、两份dyadic付款。其运行没有 Arb 或核求值，也没有根搜索。`checked_results.json` 是本次实际收据。保存核验器本身不能独立证明 Arb 的超越函数库实现；其区间语义及PF身份由已锁源码、库和理论承担，另一个全群H实现的交叉由F028负责。

## 复核与交付

只读复核命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round028/all_time_short_obstruction/check.py
```

它不修改冻结源和保存结果，不新求 K/H/G/p。不要把保留的旧 `verify_saved.py` 当作当前入口。`integrate.py --run` 在发现已存在执行收据时拒绝第二次运行。

F028 的PF比较输入均在 `results.json`：
`evaluations[0].K_at_T_upper` 和 `evaluations[1].K_at_T_upper`。
每个端点对象含 `lower_mantissa`、`lower_binary_exponent`、`upper_mantissa`、`upper_binary_exponent` 四字段；两个准确时间是同记录的 `T_upper`，不得用显示小数作输入。

本轮已停止新计算，最终白名单和源绑定见 `manifest.json`。所有写入仅在本目录，未修改023/024/026或其他代理文件，无MATLAB、扩招、提交或push。主代理保留最终数学审阅和集成责任；独立 F028 结论由其自己的报告记录，不写回本轮冻结原件。
