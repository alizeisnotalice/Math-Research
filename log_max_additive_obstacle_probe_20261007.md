# log-max 加性障碍：显式 bad 组件、稀释代数与 bounded-density 守卫

2026-10-07。先读 [解析稿](log_max_additive_obstacle_20261007.md)，再登记三个 bad 维数及三个只读旧-cell 组件。解析稿在执行前增加了显式 positive-width 证书，另存 registration amendment；原登记没有覆盖，维数、eta、tau、finite J 与 nu 预算不变。主执行正常 exit0，无 ongoing session，耗时 .149524 秒。旧 MC、旧 probes 和 oracle 没有执行。

**数值没有构造或评估未知 compact near-opt 来源，也没有报告真实 near-max 数值。** 本次只验显式组件的真实 score/严格常数 gap、条件稀释系数，以及旧真实 finite-J cells 上的冻结加性接口。它不证明一般维数阶数，也不否定一般 weak/source-square 或 actual history 预算。

## 1. 完整 bad 来源与连续真赢家

三轮 n=4,16,64，固定 `a=1,b=2,tau=1,eta=1/10`：

\[
\mu_{bad}=\mathbf1_{[-2,2]^n}\,dy+\eta\delta_0,
\qquad W_{bad}=4^n+\eta.
\]

完整背景全部计入来源质量。连续 `[1,2]` 下，原严格 `E=Q_2`，真实唯一赢家 `R=max(1,2||x||inf)`。背景在 E 上对所有允许尺度的响应准确为 1；捕获 spike 后信号随 R 严格下降，Q2 外不超阈。故

\[
Z=1+n\log2,\quad
S_{bad}(0)=\frac1{1+\eta}+n\log2
+\log(1+\eta/2^n)-\log(1+\eta).
\]

一般未知常数 KPhi 仅使用已证上界 `KPhi<=Z/e`，没有把它设为等号。本次核 `S0>=Z-2eta` 及 `S0-Z/e>0`。

对数常数用 Fraction 区间：ln2 的 arctanh 正级数 120 项及几何尾界；`ln(1+s)`（s<=.1）的偶数 alternating 120 项及下一项尾界；e 的正级数至 80! 及几何尾界。每次运算都向外舍入到 binary256bits；区间端点是 exact rationals。特别核 `e_lower>2`。下表小数仅显示近似，证书是 results 中的 rational lower/upper endpoints。

| n | 连续 S0 | S0−Z/e 的严格下界（约） | actual finite-J S0 | finite S0−Z/e 下界（约） |
|---|---:|---:|---:|---:|
| 4 | 3.592600001 | 2.204742170 | 3.479901125 | 2.092043294 |
| 16 | 11.904137144 | 7.456344144 | 11.429715977 | 6.981922976 |
| 64 | 45.175200285 | 28.487666608 | 43.247003276 | 26.559469599 |

finite J 是另行预固定的真实共同族 `R_l=1+l/(8n)`，l=0,...,8n；分别 33/129/513 个尺度，完整来源不变。其真赢家是最小能捕获 spike 的共同尺度，故精确 shell-sum

\[
S_J(0)=\frac1{1+\eta}
+\sum_{l=1}^{8n}\frac{R_l^n-R_{l-1}^n}{R_l^n+\eta}.
\]

672 个 shell terms 分别精确形成 Fraction，再向外包围并相加；真实有限包络 `Z_J=1+sum delta(R^n)/R_l^n` 也核 `Z_J<Z`。因此 finite score 对 continuous Z/e 的严格 gap 是独立 finite-family 证书，**不是把一个 radius grid 当连续真赢家或认证连续网格误差**。

实际连续 atomic bad 泛函另用闭式

\[
\Phi_{bad}=(2^n+\eta)\log(1+\eta/2^n)
+\eta n\log2-\eta\log(1+\eta).
\]

其值约 .368039705 / 1.199504547 / 4.526610938；除以完整 Wbad 的比值约 .001437094 / 2.79281e-10 / 1.33025e-38。这些是明确的 bad-only 数值，不能冒充全局 near-opt 常数或某个 near-max 输入。

## 2. 显式正宽度 L1 下界仅作解析 factor 证书

执行前 amendment 采用解析稿 §5.1 的 packet **半边长** `ell=1/(32n)`，质量 eta 均匀在 `[-ell,ell]^n`，完整背景不变。原稿证明多个 partial 坐标也满足 `m'/m>=1/(4ell)>n/R`，完整捕获端是真赢家；没有依赖最大坐标唯一或未证点态 mollification。

对每个 packet 来源点 y，已证

\[
S_\ell(y)\ge(1-2\ell)^n S_{bad}(0).
\]

本次数值只精确核 rational factor `f=(1-2ell)^n>=15/16`、`1/(4ell)-n>0` 和 `f*S0-Z/e>1/4`。三个 certified gap 下界约 1.985412663 / 6.733740353 / 25.749343652。**这是引用真实 packet 赢家解析证明后验常数，不是本 probe 数值求出了 positive-width 的全空间 S 或 max。** atomic 闭式与 L1 lower bound 分开保存。

## 3. 未知 near-opt 的复制稀释只核符号条件

同阈值 compact 近优 component 的存在由 sup 与 compact 截断逼近给出，见解析稿 §6；没有一个已知坐标/形状可数值代入。results 中 `WA,Kg,K,actual N` 都保留 null。条件是 `WA>0`、`Kg>=K-epsilon/2`、`K<=Z/e`，且各 translated compact sources 的 response halos 与完整 bad 背景的 halo 两两不交。

这些几何前提下，完整 source mass 一次计为 `WN=NWA+Wbad`，bad 热点的真实 S0 不变；完整 Phi 可加。于是

\[
\Phi/W_N\ge K-\epsilon/2-K\,W_{bad}/W_N.
\]

每轮预固定 epsilon=`1/8,1/32,1/128`，只核以下**条件代数系数**：用 certified Z_upper 选

\[
N W_A/W_{bad}\ge2Z_{upper}/\epsilon,
\quad\delta_{bad}\le\frac1{1+2Z_{upper}/\epsilon},
\quad\epsilon/2+(Z/e)_{upper}\delta_{bad}<\epsilon.
\]

九组都 Fraction 严格通过。对实际近优输入仍须满足解析存在/分离前提；这些参数不是执行了该输入。复制整数的充分下界写为 `ceil(2 Z_upper Wbad/(epsilon WA))`，未知 WA 故不报告 actual N。`d_N=NWA/tau+4^n+eta` 同样只是完整质量的准确符号表达。

atomic spike 的全来源质量占比上界是 `eta/Wbad * delta_bad`，趋小而始终正；其它副本和背景不被免费化。热点评价的保持来自解析的 compact support/halo 分离，不来自本次空间模拟。这个机制与 source-L1 near-flatness 或弱 Lebesgue 障碍均无冲突：异常 component 的质量占比被稀释，热点的局部 profile 不因此降低。

## 4. 旧完整 cells 的真实 bounded-density nu guard

另三轮 n=1,2,3 **只读**保存的 [log-max variation cell profiles](log_max_variation_probe_20261007.md)，不重跑 partition、旧 oracle 或 truewinner。源 nu 是完整有界密度 `nu=tau1_G dy`；G 取旧完整 receiver bounding rectangle 每轴外扩 b/2。本批 G 都是 `[-2,6]^n`，nu mass 为 2/8/32。对旧 E 的每个 cell，exact bounds 证明所有原 selected queries 在 G 内，所以 `A_R nu=tau`。

令 `p=tau/M<1`，按旧 E 完整积分：

\[
D_\nu=\tau\sum_E vol\,p=\int S\,d\nu,
\qquad Q_\nu=\tau\sum_E vol\,p^2\le D_\nu.
\]

| n | 读取全 cells / 原 E cells | V=nu mass | Dnu | Qnu |
|---|---:|---:|---:|---:|
| 1 | 10 / 1 | 2 | 1/8 | 1/16 |
| 2 | 256 / 13 | 8 | 35/256 | 275/4096 |
| 3 | 10648 / 789 | 32 | 36041/131072 | 28931241/134217728 |

这些体积、containment、质量、p 和 Q<=D 为 exact Fraction。添加完整 nu 后质量为 `1+tV`，保存 t=`1/8,1,4` 的所有值；没有把 V 从来源账中删除。

九个冻结原 E 的 loggain

\[
\tau\sum_E vol\log(1+tp)
\ge tD_\nu-t^2Q_\nu/2
\]

用 160-digit Decimal 数值通过，最小 residual 约 1.94367e-5；**loggain 不是 interval certificate**。n3,t4 的 RHS 为负 `−10478249/16777216`，这行只检查原公式，不解释成强正预算。输出只含 frozen 原 E gain：不求新完整最大值/新 E、不构造或测试 `A={S>K+u}` 本身，也不把 `Qnu<=Dnu` 组件核验当成弱障碍或 near-opt 定理的证明。此 guard 显示 capped-density 接口无需添加自由的 a^-n 因子。

## 文件与复现

- [registration](log_max_additive_obstacle_probe_20261007_registration.json)，[执行前 amendment](log_max_additive_obstacle_probe_20261007_registration_amendment.json)：原登记及解析修订记录都保留，没有依据结果改参数。
- [script](log_max_additive_obstacle_probe_20261007.py)，[results](log_max_additive_obstacle_probe_20261007_results.json)：全部 rational endpoints、finite shell terms、L1 factor、symbolic dilution 与 nu checks。
- `log_max_additive_obstacle_probe_20261007_n1/n2/n3_nu_frozen_profiles.json.gz`：旧真实 cell indices/M/R/vol、nu p 和全部 frozen gains，可重构该后处理积分。
- [receipt](log_max_additive_obstacle_probe_20261007_receipt.json)：输入/版本 hashes、终态 exit、scope 和文件 hashes。

results SHA256 `ffd882fda6b1f3cebbc5b38e3003318b2eadd94c146bdb0fbfcdf37ac52a8d18`。3bad+3nu 组件没有任何新 MC 随机样本或 seeds，不调用 whole-history/FIRST/CPGP/soft 门。理论 pointwise near-max 障碍由解析稿证明，本次数值验收其显式输入与关键收费组件，未发现或构造 general weak/actual geom 反例。
