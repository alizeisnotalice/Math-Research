# finite-scale log-max 乘性变分：三轮完整 Lebesgue cell 守卫

2026-10-07。预登记后执行 n=1,2,3 三轮，4/6/10 个正 atom labels、2/3/4 个共同有限尺度，全部 source 质量 W=1，`sum a_j v_j=0`，`|v_j|<=1`。六个预固定 t 为 `±1/8,±1/32,±1/128`；原权重、tau 和 J 没有按结果调整。唯一执行 session37913，exit0，耗时 1.608175 秒，未调用旧 MC 或 arrival oracle。

完整 Lebesgue partition 上的 Fraction 代数全部通过；18 个 Phi 变分下界和三轮可选 q=2 检查的 160-digit Decimal 残差均非负。**log 运算没有向外舍入，数值残差不是 interval certificate 或一般引理的证明。** 这是 actual finite J 而不是连续 `[1,2]` 窗口的验算，也不作维数阶数推论。

## 输入与完整 partition

输入都取 cube corners 与一个远处半质量 atom，corner 的一个点拆成两个 coincident 正贡献。具体坐标、权重与方向见 [registration](log_max_variation_probe_20261007_registration.json)；所有 perturbations 共享完整 source/坐标及同一 J，不取 componentwise max。

coincident labels 上预设 `v_j=±1` 不能同时当成同一点的两个物理可测 v 值。实际扰动先合并该点质量，对应唯一物理方向 `vbar(y)=sum_{j:y_j=y} a_j v_j / sum_{j:y_j=y} a_j`，仍在 `[-1,1]` 内；本批拆分点的 vbar 恰为 0。合并前后同一 `mu_t`、完整响应及 Q 全不变。同点的 S_j 完全相等，所以 `sum_j a_j v_j S_j=sum_y mu({y}) vbar(y) S(y)`，两种 D 一致。标签上的相反数只是 bookkeeping，没有给同一物理来源设置两个不同的 v(y)。

| n | labels / 物理点 | J | tau | cuts 每轴 | cells | E cells | M=tau cells | 正响应 tie cells |
|---|---:|---|---:|---|---:|---:|---:|---:|
| 1 | 4 / 3 | 1,2 | 1/4 | 11 | 10 | 1 | 6 | 2 |
| 2 | 6 / 5 | 1,3/2,2 | 1/8 | 17,17 | 256 | 13 | 48 | 12 |
| 3 | 10 / 9 | 1,5/4,3/2,2 | 1/16 | 23,23,23 | 10648 | 789 | 378 | 32 |

所有 atom 的 `y_ji±R/2` 坐标面都列入 exact Fraction cuts 并去重。相邻 cuts 的开区间做 Cartesian product，cell midpoint 处各尺度的捕获 labels 固定，故完整响应在该 cell 上常数；cell 体积 exact。端面是 Lebesgue 零测集，cube capture 仍用原闭边定义。extreme cuts 外没有任何尺度响应，贡献为零；bounding rectangle 内的零响应 cells 也保存。所有 cells 的总体积与 bounding rectangle 完全一致。

原/扰动响应均使用全部 labels：`A_R(mu_t)=sum_capture a_j(1+t v_j)/R^n`。逐 cell 对共同 J 精确求 max，同值选最小 R。对 `M=0` 的外支持 cells 不计正响应 tie；其贡献准确为零。每组有 positive-volume threshold plateaus 与 ties，coincident labels 的 source profile 相等也精确核验。

## 原赢家上的一阶和平方量

取 `E={M_J mu>tau}` 严格阈值，保留其原赢家 R，令 `m=sum_capture a_j`、`m_v=sum_capture a_j v_j`。于是

\[
p=\frac{A_R(v\mu)}{M}=\frac{m_v}{m},\quad |p|\le1,
\qquad S_j=\tau\sum_{cell\subset E}
\frac{vol(cell)\,\mathbf1_{j\in Q_R}}{m}.
\]

该 S 正是题给 `tau sum_E vol h_R(y_j)/M`，不是漏掉 R^n。完整原分母 m 从全部捕获质量求取。按完整 cells 累积，检查

\[
\sum_j a_jS_j=\tau|E|,\qquad
D=\sum_j a_jv_jS_j=\tau\sum_E vol\,p,
\qquad Q=\tau\sum_E vol\,p^2.
\]

| n | 原 E 体积 | sum aS=tau E体积 | D（两种表达一致） | Q |
|---|---:|---:|---:|---:|
| 1 | 1 | 1/4 | −1/8 | 1/16 |
| 2 | 5/2 | 5/16 | −1/16 | 5/256 |
| 3 | 467/64 | 467/1024 | −35/8192 | 251/16384 |

这些都是 Fraction 等式。原 E 不含 plateau；D 是选定原赢家给出的线性收费，不能据本表称其为 ties/plateaus 下唯一的 Frechet 导数。

## 真扰动赢家与变分检查

所有扰动 atoms 正且总质量仍为 1。每个 cell、每个 t 都先精确检查

\[
A_R(\mu_t)=M(1+tp),\qquad
M_J\mu_t\ge M(1+tp).
\]

再在**完整 cells** 上求实际扰动 max 的

\[
\Phi_\tau(\mu_t)=\tau\sum_{cells}vol\,
\log_+(M_J\mu_t/\tau).
\]

这包括新增 E 及所有原 plateau，并不只对原 E 积分。每个 RHS

\[
tD-\frac{t^2Q}{2(1-|t|)^2}
\]

先用 Fraction 求取；Phi、差值和 residual 用 160-digit Decimal。每轮六个 residual 的最小值约为 n1 `.000974695212`、n2 `.001338809030`、n3 `.001205525677`，均正，完整 raw decimals 保存。

n3 的 t=1/8 有 79 个 cell 的原最小赢家改变，并有 54 个原 E cells 退出；对应扰动仍满足精确 frozen-winner domination 和积分数值下界。这避免把固定原赢家误当扰动后的实际赢家，或只选择不切换/不退出的 cells。

每轮 6 个积分变分检查、1 个可选 q=2 检查；逐 cell/t 的 exact frozen-winner domination 检查数分别为 60/1536/63888。六个 t 的 winner-changed cell 数范围为 n1 `1–1`、n2 `6–10`、n3 `16–79`。所有 counts 都是完整 partition 的 cell 数，不是采样概率。

可选 `q=2` 核验使用严格集合 `M>2tau`，数值比较 `Phi>=tau log2 |M>2tau|`，即题给 `(log q)/q ·q tau |M>q tau|`。n1 的该严格集合为空（原有一个 `M=2tau` cell）；n2/n3 的体积为 `1`、`245/64`。没有为了令可选检查非空修改阈值。

## 可复算与边界

- [script](log_max_variation_probe_20261007.py)：无随机数、旧 oracle/MC import 或运行；Fraction 几何与 algebra、160-digit Decimal log 分离。
- [results](log_max_variation_probe_20261007_results.json)：S_j、D、Q、各 t 的 Phi/RHS/residual、真实赢家变化、新增/退出 E counts 与原输入。
- `log_max_variation_probe_20261007_n1/n2/n3_cell_profiles.json.gz`：全部 cells 的 bounds/volume/captures、每尺度完整 mass/response、原实际赢家及所有实际扰动赢家；可重构完整 Lebesgue 积分。
- [receipt](log_max_variation_probe_20261007_receipt.json)：终态执行与文件 hashes。

登记 SHA256 `f010e1099d823c124c99224551eca8e4ac69eaac9ff755edb679d9a98f249edc`；results SHA256 `35f5afa00dc2e8135e72221ffc8455681228afbfee65c966e0a4b3f68ca0cf96`。本次完整 finite-cell 核验不等于只抽有限 receiver 点，但仍只是三组有限 positive atomic inputs 的数值/代数验收，不代替任意来源/连续窗口或一般 log-max 变分定理。
