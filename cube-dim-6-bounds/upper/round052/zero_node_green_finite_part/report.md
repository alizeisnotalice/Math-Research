# B052：零 Fourier 节点的标量 Green 有限部分

本轮闭合了一个新理论接口：对任意正、满秩的六维有限 Cauchy frame，真实 Green 值有限当且仅当点在真实 $S_2$ 外。在完整、已验证的有理 fan 上，固定点的六节点 Hermite DD5 可以直接提取常数项，得到

\[
G(x)=\frac{\Gamma(x)}{128\pi^5},\qquad x\notin S_2.
\]

这里 $\Gamma$ 是给所有非零节点赋 sign 值、零节点赋 0、所有高阶归一 jet 赋 0 后的完整加权线性泛函。它不是“sign 在零点可微”的说法，也不支持零节点上的 Green 梯度或 Hessian。

本轮没有求任何完整空间点的 $G,H,p$，没有运行 native 引擎、旧 baseline、完整 tuple 枚举、MATLAB 或新时间参数。公开/内部 $C_6$ 上界仍为 64；本轮不是新界，也不是全域 Green floor。

## 严格证明与适用边界

完整证明见 [theorem.md](theorem.md)，关键链条如下。

1. **一般正 frame 的解析控制。** 满秩六列给显式 coercivity 常数，得到全空间 $p_t\le C_\infty t^{-6}$，处理大时间。若 $\operatorname{dist}(x,S_2)\ge\delta$，参数纤维中“大参数”列的秩必须至少 3；选出的独立三元组才可用于 coarea。其三个 Cauchy 尾的 $t^3$ 与真实三维商密度的 $t^{-3}$ 抵消，给 $p_t\le M_\delta$。并未把任意三个大参数误当独立。文中给出 Gram/Jacobian 常数，以及只用有理列公共分母的可核粗上界。
2. **硬拒绝真实 $S_2$。** 选包含该点的独立两列，补为六列基；限制其他参数在 $[-t,t]$，正推前给 $p_t(x)\ge c_xt^{-2}$。因此整个 $S_2$ 上 $G=+\infty$，包括原点。本轮给出全部 unordered pairs 的精确 Gram 残差证书规范，包含 dependent pairs；没有借用 body32 的特殊分层替代新 frame 的真实几何。
3. **先完整正 cutoff，后取极限。** 完整 fan 的 $H_\epsilon$ 有 $(2\pi)^{-6}$ 归一。固定节点的 DD5 Laurent 负幂只有 $-1,-3,-5$。由 off-$S_2$ 的有限性，完整有限和的三个负幂系数必须准确为零；常数项再乘 $\pi/2$，得到 $1/(128\pi^5)$。没有先对单 cone 取 Green 积分。
4. **零节点不可跳过。** 含一个或多个零节点的完整 tuple 仍可有非零常数项，必须保留。所有六节点均零的正式测试给 $H_\epsilon=1/(5\epsilon^5)$、常数项 0；它是明确拒绝例，绝不能据此认证 $G(0)=0$。
5. **连续性与导数分开。** 本轮证明 $G$ 在真 $S_2$ 外连续。固定点有限部分的假定高阶 jet 为零，只是代数常数项提取；它不是移动点的微分规则。现有 C013 generic gradient 的适用域未修改。

完整 fan 必须覆盖全空间，权重和实际 frame 对称性必须核定。一般 frame 不自动拥有 signed-permutation 群；若用已有全群合并权，不能再乘 stabilizer 或负配对系数。对新 sparse-four frame，本轮仅给可用的数学 scalar 接口，具体点的 all-pairs 证书和新引擎均留待另行授权与审计。

## 唯一实际验证批

注册先于计算，入口 [registration.json](registration.json)。唯一作者批 PID **52237**，exit **0**，工具句柄 chunk **17ea3f**，未产生持续会话或子进程。数学计算用时 **0.584824916 秒**（该 exec 返回用时约 0.732 秒），上限 900 秒。恰好 12 个预登记六节点 multiset，没有补点、失败重试或第二批。

每个 multiset 使用以下相互核对的精确表达：

- 六阶 Hermite 插值矩阵的最高系数行；
- 独立局部 pole-residue 权重；
- confluent DD5 三角递推；
- SymPy 对真实 atan 归一 jets 的 Laurent 展开；
- 符号变量 $T$ 下的 resolvent 有理恒等式，而非时间采样。

另在同批内核归一恒等式 $(2\pi)^{-6}(\pi/2)=1/(128\pi^5)$。保存结果见 [results.json](results.json)，每项包含准确节点、Hermite/residue 权重、真实 atan jets、四个 Laurent 通道及密度有理式系数。以下通道依次为 $(\Gamma,L_1,L_3,L_5)$：

| case | nodes | exact channels |
|---|---|---|
| 1 | $-7,-3,-1,2,4,9$ | $(1891/3326400,0,0,0)$ |
| 2 | $-5,-2,0,1,3,8$ | $(101/46800,0,0,0)$ |
| 3 | $-4,0,0,1,1,5$ | $(353/3600,-1/20,0,0)$ |
| 4 | $-3,0,0,0,2,2$ | $(-443/10800,1/18,0,0)$ |
| 5 | $-2,0,0,0,0,3$ | $(97/6480,-7/216,1/18,0)$ |
| 6 | $0,0,0,0,0,2$ | $(1/32,-1/16,1/12,0)$ |
| 7 | $0,0,0,0,0,0$ | $(0,0,0,1/5)$；拒绝 Green 值解释 |
| 8 | $2,2,2,2,2,2$ | $(0,0,0,0)$ |
| 9 | $-3,-3,-3,-3,-3,-3$ | $(0,0,0,0)$ |
| 10 | $-2,-2,-2,1,1,4$ | $(17/972,0,0,0)$ |
| 11 | $-2/7,-2/7,0,3/5,3/5,11/6$ | $(-119500150175/15574317726,0,0,0)$ |
| 12 | $-5,-5,0,0,5,5$ | $(-3/6250,1/625,0,0)$ |

这些是 **12 个正式标量后核**，不是真实 frame 的 12 个空间核值。没有新完整 $\Gamma$。C010/C013 仅只读核 phase/normalization；F023 的完整历史 $\Gamma$ 及负幂取消仅作已冻结来源，没有重新运行其 writer 或 checker。

## 来源、复核和未做工作

17 个只读来源的绝对路径和 SHA256 在注册文件中。特别地：

- C010 source：e857e480efae68ec99fd82f97c64f4756b6bff465476e3a788d9c283f5b751b7；
- C013 source：516d69876325810868728ff450db7a8844dcd6f7258223e6495db819b9bda7e2；
- F023 report：2851b03ea55e4ce3e95e7bbbae0b97daa9724a01f0de051aadc586255427544e；
- B023 theorem：9cde00be134dc2ede43c401808068a634bbc30b6da244090931196408498889e。

实际使用 Python 3.15.0、SymPy 1.14.0。Python executable 和 SymPy module 的 SHA 在 [run_receipt.json](run_receipt.json)。python-flint 0.9.0 仅登记为已有环境，本批算术未调用 FLINT/Arb；新理论允许未来对有限常数项做 Arb 外包，当前没有实现完整引擎。

最小只读复核命令会重算相同 12 个正式 scalar 并与保存结果比较，不写任何文件：

~~~text
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python /Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round052/zero_node_green_finite_part/check_scalars.py
~~~

作者只实际运行过一次 --record，上面只读命令作为后续复核入口，未作为第二作者批执行。所有已注册数学检查通过；理论证明提交独立审计及主代理最终审核。外部 F052 审计不属于本轮已执行结果。

下一最小工作是独立审定本定理和拒绝规范，再另行授权一个有完整 frame/fan 与 off-$S_2$ 证书绑定的 scalar 引擎验收。这里不自动开始、不预付完整核、不推广到导数、不改变现有 $C_6$ 门槛。
