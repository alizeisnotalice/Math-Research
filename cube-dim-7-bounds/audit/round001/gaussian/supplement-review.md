# Gaussian 完整 I 解析预算补充独立审计

2026-10-27。审计对象为主集成 clone `cube-dim-7-bounds/upper/round001/gaussian/` 的 `proof.md`、`check.py`、`check-results.json`、`results.json`、`sources-hash.json`，以及 `src/check_exact_budget_round001.py` 和 `audit/round001/exact_budget_results.json`。上述七文件已全文读取并锁定哈希。本补审不修改先前封存的 theory 审计工件；其全部九文件在检查前后哈希相同，原 manifest 也重新核对通过。

**结论：接受 Gaussian 解析完整 I 上界。** 用 768-bit Arb 独立重评后得到的精确有理上端严格小于主集成 256-bit 的锁定上端，因此该旧上端可继续保守使用。代入条件性 g=1/179000，精确有理计算确认成本六次幂小于 128^6。此 g 下界尚未证明；这里没有新 C7 上界结论。MC 结果仍只属于 candidate diagnostics，其完成状态或样本均值不是解析证据。

## 1. Gaussian subordination 与共同时间

设独立 S_v>0 的 Laplace 变换为 exp(−sqrtλ)，Σ=Σ_v S_v vv^T。给定 Σ 的 Gaussian 协方差是 2t²Σ，其特征函数是 exp(−t²ξ^TΣξ)。对 S 平均，因 t>0，得到 exp(−tΣ_v|v·ξ|)，正好是原 64 个 weight-one 未单位化方向的共同时间 Cauchy 模型。这里没有把 t 改为 64 个可独立选择的时间，也未单位化方向。

条件密度为 `(4πt²)^(-7/2) detΣ^(-1/2) exp(−Q/(4t²))`。t 积分以 u=Q/(4t²) 作非负 Tonelli 换元，常数为 Γ(3)/(4π^(7/2))=`1/(2π^(7/2))`，且 Q 的幂为 −3。对 x≠0 逐个 SPD shape 均有限；对期望可保留扩展正值，因而也兼容 S3 的 G=∞。这一步由本维模型重新推导，不以旧 n6 theorem 或历史数值状态为证。

## 2. 全 boundary 恒等式、Minkowski 与 896

在七维 cube 的每个开 face 上，z=rx 的 Jacobian 为 r^6；各 face 相交边界为零六维测度。先对 r 积分，有

`∫R7 exp(−Q(z))dz = Γ(7/2)/2 · ∫∂K Q(x)^(-7/2)dσ(x)`。

另一方面 Gaussian 积分等于 π^(7/2)sqrt(detΣ)，故 boundary 因子严格为 H=`2π^(7/2)/Γ(7/2)=16π³/15`。这个恒等式对每个 SPD Σ 成立，不要求一个 shape 有 cube 对称性。

q=7/6 时，单 shape 的 `Q^-3` 的 q 次幂是 Q^(-7/2)。全 boundary 的 Lq 范数取根后，乘 detΣ^(-1/2)，剩下 determinant 指数 `−1/2+1/(2q)=−1/14`。因此 Minkowski 给 `||G||Lq(∂K)≤C H^(1/q) E detΣ^(-1/14)`。为避免先假设期望可积，作者限制在递增的事件 `{N^-1 Id≤Σ≤N Id}`，随后 MCT；这些事件穷尽所有 SPD Σ，论证合法。

只有在完成期望后，原 law 才按坐标符号与置换对称。原 frame64 的方向线和等权保持不变，因此 `∫∂K G^q=14·64I=896I`。作者没有对每个随机 shape 的单一正 face 强加对称。故

`I≤C^(7/6)H M^(7/6)/896 = M^(7/6)/(105·2^(25/6)·π^(13/12))`。

本审独立 SymPy 恒等式确认此系数，Arb 则从未简化的 C、H、896 路径评价，与简化式重叠。

## 3. Gamma/Dirichlet、AM–GM 与矩上端

S 的 inverse-Levy 变换 Y=1/(4S) 的确是 Gamma(1/2, rate1)。R=ΣY_v 服从 Gamma(32,1)，P=Y/R 服从 Dirichlet(1/2,…,1/2)，且独立。A(P)=Σ vv^T/P_v，Σ=A/(4R)，于是 `detΣ^(-1/14)=2sqrtR detA^(-1/14)`；R 的半阶矩给 D=`2Γ(65/2)/Γ(32)`。点值表达式 T_x 同样缩放为 2sqrtR，原 proof 的径向消元正确。

Cauchy–Binet 给 `Σ_J c_J=det(VV^T)=64^7`。权 α_J=c_J/64^7 的 i 边际由 determinant derivative 得 `Σ_{J∋i}α_J=7/64`，不是 1/64 或 7 的其他倍数。Weighted AM–GM 因此给 `detA≥64^7 ∏P_i^(-7/64)`，再取 −1/14 次幂，得到

`Z≤64^(-1/2)∏P_i^(1/128)≤1/64`。

Dirichlet beta 积分中 64 个 exponent 的和为 1/2，分母参数由 32 改为 65/2，故

`E∏P_i^(1/128)=Γ(32)/Γ(65/2) · [Γ(65/128)/Γ(1/2)]^64`。

与 D 相乘后两个 Gamma 因子消去，而 2/sqrt64=1/4，最终严格为

`M≤MΓ=(1/4)[Γ(65/128)/Γ(1/2)]^64`。

该 MΓ 是确定性解析期望上界；并不依赖 MC 均值、observables、缓存恢复或八个点诊断。完整 I 的所有时间和所有 boundary 奇角都已包含，没有 cap 或小时间遗漏。旧 n6 theorem 只作为作者的方法来源，不承担这里任何七维证明步骤。

## 4. 实际独立 768-bit 复评与精确有理比较

本审先登记固定任务，再运行独立 `check.py`。未导入、未执行任何作者 checker；未读 observables.npz，未重跑 MC，未优化 M，未新增 Green 节点。与作者共享 python-flint/Arb 库，不声称不同算术内核；不同精度与 log-Gamma 安排提供独立复评。实际版本为 Python3.15.0、SymPy1.14.0、python-flint0.9.0。一次完整执行 exit0，无实现修复或数学断言失败，脚本计时 0.01583466600277461 秒。

MΓ 用 `exp(64(logΓ(65/128)−logΓ(1/2)))/4` 评价，I 用 `C^(7/6)H MΓ^(7/6)/896` 评价。所有输入参数是精确整数或有理数，精度 768 bits。区间上下端各经 `.lower().fmpq()` / `.upper().fmpq()` 序列化为精确有理数，完整数值留在 results.json。复评 I 表达式约为 `9.79365339029564320881646517954×10^-6`。

锁定旧上端为

`74319631963157785256240159489393050185362460367580814299903508439387034491195 / 7588550360256754183279148073529370729071901715047420004889892225542594864082845696`。

本审用 Fraction 精确证明 `0<I768_lower<I768_upper<I256_upper`；results.json 保留旧上端减新上端的正有理差。还核对主集成 exact-budget 文件序列化的成本六次幂和 gap，均与本审独立有理运算完全一致。此比较不以近似小数或“两个区间显示相同”代替严格有理大小关系。

## 5. 条件 g 接口和结论边界

令 g_target=1/179000。用锁定旧 I256 上端而不是更紧新端，精确验证

`7^7 I256_upper^6 / g_target^7 < 128^6 = 4398046511104`。

左端约 4278838307752.45756749，对应条件成本约 127.41512341408573；这些小数仅帮助阅读，接受依据是保存的精确有理差严格正。新 768-bit I 上端的条件成本还严格更小。

此处只证明“若 g≥1/179000，则结合已接受有限源比较可得成本小于128”。八个 face 点的样本 minimum、MC 的 M mean 或经验 U≈19.7 均不证明全局 essential g。主集成 exact-budget JSON 已明确 global_g_target_proved=false 和 new_C7_upper_bound=false，本审逐项核对并保持该状态。

Gaussian 原 proof 的 anchors/S3 辨识与有限诊断在本次阅读范围内没有发现与已接受理论冲突；本补审不新增点值积分认证，也不独立认证 MC 的工程执行。可通过解析完整 I 预算减少后续 I 认证负担，但全域 g 门槛仍是开放义务，不能把本补审称为新 C7 上界。

复现：从本独立 clone 根运行

```
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-7-bounds/audit/round001/gaussian/check.py
```

本补审仅写 gaussian 子目录；源、sealed theory 与全部已封存文件保持不动。无 commit/push、MATLAB 或子代理。
