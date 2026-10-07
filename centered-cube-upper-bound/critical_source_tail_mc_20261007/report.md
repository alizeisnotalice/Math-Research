# 临界参数下完整来源占用：三轮新的原空间 Monte Carlo

2026-10-07。三轮全部完成固定预算，但**所有统计下界为零，不能判断二阶/一阶比值的维数阶数或证明 source-tail 合同**。n=64 的 outer 二阶 cross-score ESS 仅 3.6846；其原始点估计不满足真实矩一致性，显示有限样本噪声。以下保留 raw values，不投影修正、不扩大样本、不调阈值。

本次与旧 24 点有理守卫不同：实际新采样完整来源 `mu` 和条件 receiver，估计原占用 `S(y)`、`I1/W`、`I2/W`。旧 24 点、旧 oracle main 与旧 terminal run 均未重跑。唯一主运行 session `42954`，正常 exit 0，总耗时 69.222522 秒；每轮固定 256 sources、每 source 两个独立 replicas、每 replica 8 receivers，共 768 sources / 12,288 receivers。

## 1. 冻结完整来源、阈值与目标

继续使用已核 [critical arrival 输入](../critical_source_tail_arrival_20261007.md)：四个原 tensor grids 正混合，`k=(1,2,4,8)`，`w=(1/2,1/6,1/6,1/6)`，`A=16n`，一维点集 `{-Ak,...,Ak}/k`，`q_k=32nk+1`。完整质量 `W=1`。

来源一次采样：均匀整数 `0,...,5` 的前三值选 component 0，其余各选一个 component；该 component 的每坐标均匀采整数 `[-Ak,Ak]` 再除以 `k`。因此采的是完整正 `mu/W`；重合物理点上的 fractional component 贡献自然相加。既未删除重合点，也未按选中 component 另求赢家或重新归一化。

阈值为旧登记的 exact 值，不由本批输出调整：

\[
\tau_n=\tfrac12 q_1^{-n}(3/2)^{-n}2^{t_n},
\quad(t_4,t_{16},t_{64})=(4,14,47).
\]

`t_n` 来自参考内部 count 的确切 `Bin(n,1/2)` 上尾首先降至 `n^-2` 的位置；这只是固定参数映射。**完整混合赢家事件和给定 source/cone 的参考 count 事件均不是该 binomial 概率。** 本批逐坐标使用实际有限格子、source 位置、core/shell 宽度所给的 high intervals 和概率。

定义完整连续窗口原对象

\[
M(x)=\sup_{1\le L\le2}h_L*\mu(x),\quad
R(x)=\min\operatorname*{argmax}_{1\le L\le2}h_L*\mu(x),\quad
E=\{M>\tau\},\quad g=\frac\tau M\mathbf1_E.
\]

这里 `h_L=L^-n 1_{[-L/2,L/2]^n}`。完整 source occupation 及目标为

\[
S(y)=\int g(x)h_{R(x)}(x-y)\,dx,
\quad I_1=\int S\,d\mu=\tau|E|,
\quad I_2=\int S^2\,d\mu.
\]

原 receiver/source incidence 全部保留；本批不是 FIRST/soft winner 数值，也不增加任意 receiver 门。

## 2. Lebesgue receiver prior 与实际 sampler

令 `Cn=1+n log2`，对每个 source `y` 使用

\[
q_y(x)=\frac{\sup_{1\le L\le2}h_L(x-y)}{C_n}.
\]

其 core 质量为 `1/Cn`：均匀 `y+[-1/2,1/2]^n`。shell 质量为 `n log2/Cn`：`log r` 均匀于 `[0,log2]`，有符号 face 均匀选自 `2n` 个面，free offsets 均匀于 `[-r/2,r/2]`，固定 face offset 为 `sign*r/2`。规范 cone measure 下 `dx=n r^(n-1)dr dnu`，所以 shell 中 `q_y=1/(Cn r^n)`，不存在额外 `2^-n` 因子。

给定 `y`、core/shell、`r`、face/sign，对每个自由坐标以**有限 k=1 grid 在参考边长 3/2 的捕获数恰为 2**定义 high set。原均匀区间的实际 high 比例为 `p_i`，`Z_i=1-p_i+p_i T`。预固定 `T=(8,7,47/17)` 随三轮变化，整 product 倾斜的 likelihood 为

\[
L=\frac{T^{\sum_{i\in free}H_i}}{\prod_{i\in free}Z_i}.
\]

整个 conditional proposal 是 `.25 base + .75 tilted product`，只抽一次全局 branch，故

\[
w=\frac{d prior}{d proposal}
=\frac1{.25+.75L}\le4.
\]

core 的所有坐标自由；shell face 的参考 `H` 仅记录诊断，**不进入 RN**。它在条件分布中固定，不可把 face 项作为一个新的自由 Bernoulli。base branch 保留完整原支持，包括有限 source 边缘 high 比例为 0 或 1 的情况。tilted branch 先按 `p_i T/Z_i` 选 status，再在该 status 的真实 interval union 按长度均匀抽样。

在新 sampler 执行前，用 12 个新 rational configs（每 n 四个，含 source 支持边缘、core 与 shell）精确核查 partition、RN 质量 1、首 count 矩与 RN 二阶矩不超过 4；未沿用旧 9 configs 充作实现验收。每个实际 draw 又检查 partition 总长度、真实点 high status、free 坐标支持。保存逐 draw 的 branch、r/hex、face/sign、free `p_i`、high statuses、log normalizer、logL、RN 与 capture。

有限机器均匀数、interval 选择、log/RN 和 score 运算没有向外舍入认证。上述 proposal/RN 恒等式在理想连续 law 上精确；机器实现的概率误差未被 interval 包围。H02 只用于该倾斜换测度恒等式，没有使用 saddlepoint/Mills 前因子或把参考 binomial 当真实 cone 分布。

## 3. 所有 arrival 的连续真赢家与舍入范围

新实现把 sampled offsets 的浮点值冻结为 exact dyadic numbers；真实输入 `x=y+offset` 用统一二进制分母构造，避免大坐标的 `float(x-y)` 消去。对每个 coordinate/component 枚举 `[1,2]` 中所有原子到达边长，以整数 key 去重，同距离的全部 `(component,coordinate)` counts 同时更新。包含端点 1 和 2。

完整捕获质量始终为 `sum_k w_k product_i counts_ki/q_k^n`。mass/R^n 的候选比较、相等 ties、`M>tau` 和选中完整 source 点的捕获资格均做整数 cross multiplication，未用 radius grid 或浮点 near-screen 判定义。smallest-R tie 规范不变。

shell face offset 的 binary 除 2 精确，free offsets 用 rational support 检查，故实际 exact offset max 恰为 `r/2`；`q_y` 的 shell 半径与真实 receiver 构造一致。source capture 判定使用 exact `2max|offset|<=R`，含等号。支持边缘 source 点数为 9/15/10，coincident-component membership 逐 source 保留。

只对 **9 个新实际样本**（每轮三点）调用旧模块的纯 `arrival_oracle` Fraction 函数交叉核验，全部 winner R 和 M 一致；旧模块 main/旧 probes 不执行。另一个只读后处理独立重构全部 12,288 保存 receiver 的 exact x、closed capture counts、完整质量、E/capture、RN/score，全部一致；该后处理没有重新取 supremum 或采样。

`1e-10` near threshold / near max 与 `1e-12` partition boundary 的观测计数均为 0，exact 多候选 winner ties 为 0。这些计数不认证实际连续理想 draw 与冻结 dyadic 之间的误差、也不认证概率 law。完整 dyadic winner 已不依赖这些 heuristics。

## 4. 两个独立 replicas 与分布无关范围

逐 receiver 的 importance score 为

\[
Z=\mathbf1_E\mathbf1_{y\in Q_R}\frac\tau M
\frac{h_R(x-y)}{q_y(x)}w,\qquad0\le Z\le4C_n.
\]

core 的 `h_R/q=Cn/R^n`；shell 的值为 `Cn(r/R)^n`。所以在**理想连续、精确运算模型**中 `E[Z|y]=S(y)`，而 `S(y)<=Cn`。这一步把 receiver 的 Lebesgue 积分明确转为条件 cone/core proposal，不把来源概率直接当 receiver 体积。

对同一 source y 的两 replicas 采用分开的独立 RNG streams，分别独立包含 core/shell、r、face/sign、product-mixture branch 和所有 free draws。每 replica 8 draws 的均值记为 `S1,S2`，估计

\[
\widehat I_1=\frac1{256}\sum_s\frac{S1_s+S2_s}{2},
\qquad\widehat I_2=\frac1{256}\sum_sS1_sS2_s.
\]

理想模型条件于 y 时，cross product 的期望为 `S(y)^2`，避免 empirical-square 的 conditional 方差上偏。两估计固定 N 时无偏；本次 N 全完成，没有数据/成本停时后仍称均值无偏的情况。

预注册 Hoeffding union：outer 两均值 × 三轮 × 全 prefixes 1,...,256，共 1536 项，以总 delta=.025；观察范围分别 `4Cn`、`16Cn^2`。每个 S profile 用两 replicas 的 16 draws，768 profiles 再分 delta=.025；与 outer union 总覆盖至少 .95，**仅指理想精确概率模型**。以已知 `I1<=Cn,I2<=Cn^2,S<=Cn` 截取统计目标区间。没有 normal/bootstrap 区间充作高概率证书。

若外部中断，注册容许 completed prefix 的 union 区间，但主代码没有自动恢复截断 gzip/输出停时收据；本次正常固定 N 终态，不调用这条备用路径。

## 5. 原始结果、有效样本与矩一致性

| n | E / 4096 draws | E 且 source 捕获 | I1/W raw | I2/W raw | I2hat−I1hat² | outer I2 cross ESS |
|---|---:|---:|---:|---:|---:|---:|
| 4 | 3059 | 2238 | .687285852 | .554435448 | +.082073606 | 73.5337 |
| 16 | 2357 | 1510 | .156571867 | .024423802 | −.0000909474 | 10.3808 |
| 64 | 2457 | 1476 | .538028312 | .183132216 | −.106342248 | 3.6846 |

完整来源 W=1 的真实矩必有 `I2>=I1^2`。独立双 replica 的有限样本 cross estimator 不必满足这一逐样本约束；n=16、64 的负残差说明噪声，**不是数学反例，也不必然是实现错误**。原值及差值另存 `moment_diagnostics.json`，没有投影到可行矩区间后伪造通过。

| n | I1 simultaneous ideal CI | I2 simultaneous ideal CI | receiver RN ESS | receiver score ESS | outer I1 score ESS |
|---|---:|---:|---:|---:|---:|
| 4 | [0,2.97030248] | [0,14.23242567] | 1871.37 | 630.06 | 173.57 |
| 16 | [0,7.47316090] | [0,146.17668134] | 1176.35 | 44.26 | 33.74 |
| 64 | [0,27.98890694] | [0,2057.65838412] | 1099.74 | 27.08 | 24.17 |

ESS 是保存 importance contributions 的 `(sum z)^2/sum z^2` 诊断，不是独立样本数或置信度。RN sample means 为 1.005936 / 1.003232 / .996234；它们只是 sampler 数值 sanity check。所有 I1、I2 下界为零。由这些统计矩形区间形成的 `I2/I1` ratio 上界均无界；已有解析粗费 `I2<=Cn I1` 仍给线性 Cn 上界，但没有新平方根控制。

n=64 两个 replica 的 I1 raw means 为 .596565 与 .479492，source 的 S mean-estimate 最大 10.0051；单 receiver score 的权重退化仍明显。阈值和高 count tilting 已进入临界探路，但上述稀有贡献的精度不足。三个有限 n 的点比值 .806703/.155991/.340377 不作为阶数估计。

尾诊断固定 K=sqrt(n)、2sqrt(n)。`mean[(S_hat-K)+]` 是 nonlinear plug-in；理想 law 下由凸性

\[
E[(\widehat S-K)_+\mid y]\ge(S(y)-K)_+.
\]

故其期望上偏，和无偏 I2 cross estimator 不同。n=4 的 K=2 插件 .00517210；n=64 的 K=8 插件 .02081464；其余登记插件值为观测 0。没有额外 tail CI，所有观测 0 均不能断言解析 tail=0；单个 positive plugin 也不能断言该 source 的真实 S 超阈。

## 6. 可复算文件、hash 与结论范围

- [registration](registration.json)：预设来源/阈值、seeds、proposal、预算、停止与置信范围。SHA256 `cc6aabb10f51ef2e961405056641c1667884da214169e1f2ea47b962ed194446`。
- [run.py](run.py)、[results](results.json)：原源 sampler、整数 arrival oracle、全部三轮 raw results；results SHA256 `4462e37d95026400dda13ba6a79e3dccd65855e49887668bfeed0f630a018c26`。
- `n4/n16/n64_source_profiles.jsonl.gz`：每 source 的完整两个 receiver replicas、分母、counts、winner、E、source capture、p/RN、scores 与 S profile 理想区间。
- [saved-profile review](saved_profile_review.json)、[moment diagnostics](moment_diagnostics.json)：全保存点重构与 raw 矩一致性。后处理没有 oracle 最大值重跑。
- [independent review](../critical_source_tail_mc_20261007_independent_review.md)：另一个 agent 只读审阅概率/实现接口，未执行本脚本。
- [receipt](receipt.json)：主 session 42954 exit0、后处理 session4016 exit0及全部终态 hashes。

三份 source profiles SHA256：n4 `6aea08ddd02ed3f02f39c762c19f7476f85676c0f6b5c312de11b94a57712726`；n16 `2bc0923a5a0bcc2bdd200ea2f0639000c435477e6b47226fe040291d3581ff68`；n64 `430403c109bcb1e4083d323b7398306412c75f2c5e36d15db696e1ab81c63c84`。

本批完成了此前未做的临界完整来源占用 MC，并保存可重构的原空间 profiles。它没有证伪或证明一般 source-tail/二阶平方根预算；原 sharp lower-model 全资格、任意非负 L1 来源、actual geom gates 的推论仍未获得。当前范围是已冻结四层正 tensor grids 的真实连续 `[1,2]` atomic winner 与理想-law 数值压力，不将特殊族当一般结论。
