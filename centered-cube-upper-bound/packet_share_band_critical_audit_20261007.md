# 原连续赢家的固定单格点 packet 份额带：三轮精确后处理

2026-10-07。按先保存的 [registration](packet_share_band_critical_audit_20261007_registration.json)，只读此前保存的 24 个 exact continuous `[1,2]` 完整混合来源赢家 receiver，以及此前保存的各坐标捕获 counts。本次不执行 arrival oracle、不生成新来源或 receiver、不枚举指数数量原子。运行 exit 0，耗时 0.021813 秒。三个维数各 8 个确定性记录，不是 Lebesgue 采样；其中原 `E={M>tau}` 合格点数为 4、1、1。

## 1. 固定完整来源分配

来源仍是原四个正 tensor grid 的完整混合，`k=(1,2,4,8)`，`w=(1/2,1/6,1/6,1/6)`，`A=16n`，`q_k=32nk+1`。预固定每个原 component 的一个物理格点为一个 packet：

\[
\mu_{k,z}=M_{k,z}\delta_z,\qquad M_{k,z}=w_k/q_k^n.
\]

重合格点的 component contributions 是 fractional source 分配：在该点定义 `a_{k,z}=M_{k,z}/mu({z})`，故所有贡献相加恰为完整原 `mu`，`sum a=1`。全局 packet 质量核验 `sum_k q_k^n M_{k,z}=sum_k w_k=W=1`。固定 packet 不因 receiver、赢家或 band 改变；没有为各 band 复制来源质量。

若原闭立方体赢家捕获该点，`m_{k,z}=M_{k,z}`，fullness 恰为 1（任意固定 `eta<=1` 自动合格）；否则 `m_{k,z}=0`、份额为零，不进入 band。分配为重合点 labels 不构成物理分离或几何扩散论证。

同一 component 的捕获原子数直接用既存 `N_k=product_i count_{k,i}`，其每个被捕获 packet 的份额相同：

\[
m=\sum_k N_kM_k,\qquad p_k=M_k/m,\qquad
\sum_kN_kp_k=1.
\]

所有式子以 Python `Fraction`/整数精确检查。保存的 `receiver x,R,M,tau,E` 全部与原 exact profiles 一致，且再次核对 `M=m/R^n`、`E=(M>tau)`。使用既存闭边界 counts；未改边界定义或 floating tolerance。

## 2. 严格端点与 row 预算

采用 [一般份额带稿 §1–3](general_packet_share_band_square_20261007.md) 的互不重叠规范：

\[
\theta_j=2^{-j-1},\qquad \theta_j<p\le2\theta_j,
\quad j=0,1,\ldots.
\]

特别地，`p=2^-b` 属 `j=b`，不属于 `b-1`。每个正 captured packet 恰进一个 band；同 band 的多个 component 只合并 counts。各 component 保留 `N_k,M_k,m,p_k,j`，各 band 保留 count、actual share 及 `theta_j count_j`。

记 `B(x)=sum_j theta_j count_j(x)`，则原完整质量给

\[
B(x)<\sum_kN_kp_k=1,\qquad
1=\sum_kN_kp_k\le2B(x),
\quad\boxed{1/2\le B(x)<1}.
\]

上界严格因 `theta_j<p` 且原捕获质量正；下界可等号，当全部正 captured 份额均处各自 band 的闭上端。本批没有观察到份额恰等上端；实现仍按精确端点处理，未将近端点当等号。原 E 外行只作为 implementation audit：`E`-gated 份额及预算显式为零，未把这些行加入原 receiver 合格集合。

## 3. 三轮结果

| n | 保存点 / E 点 | 全部点 B 范围 | E 点 B 范围 | 全部点占用 j 范围 |
|---|---:|---:|---:|---:|
| 4 | 8 / 4 | 0.641586304–0.807647705 | 0.696296692–0.716560364 | 2–18 |
| 16 | 8 / 1 | 0.503662142–0.774499720 | 0.503662142 | 5–66 |
| 64 | 8 / 1 | 0.500000003364–0.737128352 | 0.500000003364 | 20–256 |

24 行、96 个 positive captured-component incidences 全部通过；每行恰有四个非空 band。n=16 的唯一 E 点 bands 为 `15,32,48,64`；n=64 为 `63,128,192,256`。这些大编号来自原格点质量与完整捕获分母；没有按观测结果调整 packet、tau 或赢家。源种子沿用保存值并逐行记录，没有新随机种子、样本或统计置信区间。

每行的四个 band 不给来源 profiles 的全局四-band 限制。同一固定 packet 在不同 receivers 上的份额可进入不同 bands；各 receivers 的 labels/bands 互斥并不令来源空间上 `S_j S_k` 消失。本结果只审计完整原 source 分配下的 row 份额费，未计算 `A_j=tau theta_j sum_i |E_i,j|` 的 Lebesgue 积分，也未计算 `int sum S_j^2 dmu` 或 `int (sum S_j)^2 dmu`。一般平方函数引理的证明另见上述稿，本次有限点核验不充当其证明或跨 band 预算。

先前固定全层 packets 在这些稀疏接收点 fullness 全失，本次预固定单点 packets 使每个 capture fullness 为 1，且 complete capture fraction 为 1；与此同时 packet 质量深度/份额 bands 变大。这说明本批不同固定分配的资格差别，不证明任意来源的有效 band 数为平方根，也不证明 actual geom endpoint。

## 4. 可复算文件与范围

- [script](packet_share_band_critical_audit_20261007.py)：纯只读 exact 后处理，不 import 或执行旧 oracle。
- [results](packet_share_band_critical_audit_20261007_results.json)：全部 receiver/source/count/份额/band 重构字段，注册及输入 hashes，三轮摘要。
- [receipt](packet_share_band_critical_audit_20261007_receipt.json)：exit、耗时与终态文件 hashes。

新 results SHA256 为 `c423183271d6395a4158c280efa9fa992f1480b34ebc6f8def6ad819d55ab3c5`；输入包含旧 capture results、三份 exact profiles 与一般份额带定义稿，均在运行前按注册 hash 核验。旧文件没有修改。24 点没有 Lebesgue 覆盖概率意义；不作来源积分、阶数拟合或一般跨 band 平方预算结论。
