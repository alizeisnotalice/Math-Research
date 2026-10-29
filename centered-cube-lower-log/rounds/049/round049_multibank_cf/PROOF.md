# 有限多 bank / 多尺度 anchor-CF 比较

先运行 `probe.py`，完成同一 GF(2⁷)、三个真实 trace basis、三个尺度、同一 query 的原模型/共享相位解耦模型有限和；然后给出以下扩展。这里的“真实”指确实按域乘法和 trace 产生标签，并非为小例证明目标 Sidon 参数。一般命题不需要 Sidon 或联合满秩。

## 1. 索引与保留的全部依赖

令有限索引集 I 包含任意 (bank,scale) 对。所有 bank 用同一 GF(q)、同一本原元 g、同一个 j、同一个 θ；bank b 可用不同独立 trace basis β_(b,i)。端点按各自 basis 与位移计算。假设存在统一界

\[
\|H_\iota(\cdot,\theta)\|_W\le U,\quad |H_\iota|\le2,
\qquad \|B_\iota(\cdot,\theta)\|_A\le V,\quad |B_\iota|\le V.
\]

若全部 basis 长度为 b，可取 U=2(b+1)；不同长度取相应最大值。B 可共享，也可作已知固定物理 Fourier 平移，循环 Fourier 系数只乘单位模因子，范数界不变。端点平滑峰须留格，与47轮相同。

写 a_ι=p/(KL_ι)>0、a*=max a_ι、W=U+V、r=Wa*<1。定义

\[
Y_\iota=1+a_\iota[H_\iota(j,\theta)+B_\iota(j,\theta)],
\quad
\widetilde Y_\iota=1+a_\iota[H_\iota(j_H,\theta)+B_\iota(j_B,\theta)].
\]

比较模型仅令 j_H,j_B 条件独立。**所有 bank、所有 scale 的 H 仍用同一个 j_H，所有 B 仍用同一个 j_B；θ也仍共享。** 所以任何跨 bank 的 trace 线性关系、联合输出纤维、rank43等缺陷均被完整保留。Y、Ỹ≥1−a*(2+V)>0。

零态上的 H 使用同一 trace 公式赋值。对任意 GF函数 F 与循环函数 Ψ，Gauss–Walsh 估计仍是

\[
|E F(g^j)\Psi(j)-EF(g^j)E\Psi(j)|
\le \lambda\|F\|_W\|\Psi\|_A,
\quad\lambda=\sqrt q/(q-1).
\]

F 可以是不同 bank、不同端点 H 的任意乘积。其范数用次乘法估计，不把这些 H 当成独立随机变量。每个 basis 拉回的 Walsh 字符都是同一个域上的加性字符；系数合并只会降低绝对和范数。因此 bank 的联合 trace rank 完全不进入这条比较。48轮实际四bank的 rank43不妨碍应用。

## 2. 任意有限 score 家族和一个合法 anchor

取家族中的 anchor α。令

\[
A_\alpha=EY_\alpha=E\widetilde Y_\alpha>0.
\]

等式由纯H、纯B的边缘均值相同精确得到。真实 Haar x、实际周期源和平移后的零均值 Fourier 项给 A_α=1；固定 θ 或有限相位混合时不能强设1。复幂始终用正数实对数定义。

对于所有索引一起收费：

\[
S_1=\sum_{\iota\in I}|\eta_\iota|a_\iota,
\qquad S_2=\sum_{\iota\in I}|\eta_\iota|a_\iota^2.
\]

置

\[
D_\alpha(\eta)=\lambda UV\exp\!\left(\frac{WS_1}{1-r}\right)
\left[
\frac{(1+a_\alpha W)(S_1^2+S_2)}{(1-r)^2}
+\frac{2a_\alpha S_1}{1-r}
\right].
\]

在合法 anchor 倾斜 dP_α=(Y_α/A_α)dP 下，定义完整向量 logY 的 CF φ_α(η)，解耦模型对应 φ̃_α。对 n 个 iid **倾斜后的坐标**，

\[
\boxed{|\varphi_\alpha(\eta)^n-\widetilde\varphi_\alpha(\eta)^n|
\le \min\{2,nD_\alpha(\eta)/A_\alpha\}.}
\]

证明与47轮逐字兼容，但 I 现在为完整 bank×scale 集。将
$Y_\alpha(\prod_\iota Y_\iota^{i\eta_\iota}-1)$
按 H/B 展开；绝对系数由

\[
(1+a_\alpha w)\left[\prod_\iota(1-a_\iota w)^{-|\eta_\iota|}-1\right],
\quad w=Uu+Vv
\]

支配。纯 H、纯 B 项在两模型下完全相同；所有混合项用上述 Gauss 协方差界。混合导数给 D_α。r<1使范数级数绝对收敛并对θ一致，故可交换求和、有限域期望及θ积分。减去 anchor 使η=0误差精确为0。再除同一 A_α，并用两个真正 CF 的模≤1及代数 n次幂差恒等式，得到结论；无需复对数。

因此费用依赖所有 bank/scale 的总ℓ¹频率，而不是“每bank有一份免费频率”。若 η=t/√v、所有 L≥Lmin 且 Σ_I|t_ι|≤T，在47轮接受参数下同样得到

\[
nD_\alpha=O_{p,L_{min}}\!\left(\frac{\log B}{\sqrt K}(T^2+T\sqrt v)\right).
\]

有限索引个数本身不另收费，前提是这个**总**ℓ¹预算及统一范数前提成立；逐坐标频率有界不能替代总预算。

## 3. 固定 mixture anchors：两种不同 n 坐标构造

对若干合法归一化 anchor law P_α，固定概率权重 ν_α≥0、Σν=1。权重不得依赖被观察的 bank 赢家、样本或 η；两个模型必须用相同权重。

**一个全局 component 标签的混合：**
$Q^{(n)}=\sum_\alpha\nu_\alpha P_\alpha^{\otimes n}$。
其 CF 是 Σν_α φ_α(η)^n，不是某一个 bank 的 CF。于是

\[
|\Phi_Q-\widetilde\Phi_Q|
\le\sum_\alpha\nu_\alpha\min\{2,nD_\alpha/A_\alpha\}
\le\min\{2,n\sum_\alpha\nu_\alpha D_\alpha/A_\alpha\}.
\]

也可用 max_α D_α/A_α 代替加权费用，不出现未归一化的 anchor 数倍。

**每坐标重新混合：**
$Q^{\otimes n}$、$Q=\sum_\alpha\nu_\alpha P_\alpha$ 的 CF 是 $(\sum_\alpha\nu_\alpha\varphi_\alpha)^n$。它通常与前一构造不同；单坐标先凸组合，再用 n次幂差，同样得到最后一个加权界。

若原先给的是非归一化正 anchor $h=\sum\omega_\alpha Y_\alpha$，须先用
$\nu_\alpha=\omega_\alpha A_\alpha/(\sum\omega A)$
作单坐标倾斜权重；两模型 A_α 相同，所以此权重也相同。若 n坐标先混合未归一化 anchor products，则权重相应为 $\omega_\alpha A_\alpha^n/\sum\omega A^n$。Haar下 A_α=1，原权恰好保留。

讨论真实 input mixture tilt 时，必须先核对该 tilt 与所用合法 anchor law 的实际关系；若它是各bank anchor的固定混合，就须保留全部 component 及原权。一个 bank 的 anchor 不能代表整个 mixture tilt。特别是原 input density 与其某尺度卷积平均不是同一个函数，不能因为两者质量1而直接替换倾斜。

## 4. 不覆盖的非线性量

已有 CF 比较的随机向量是
$S_{b,\ell}=\sum_{k=1}^n\log Y_{b,\ell}(X_k)$。
全局 bank mixture 的平均为 $\sum_b\omega_b e^{S_{b,\ell}}$，其对数是 log-sum-exp，**不是 S 向量的线性 score**。不能把它塞进某一个 η坐标，或只用一个 bank 的 anchor 直接套公式。

同理，max_L、max_bank、mixture阈值或罕见事件下界没有自动包含在这个有限频率 CF 界中。若要由联合CF处理非线性事件，仍须建立反演/平滑、频率总预算、反集中和尾部误差；anchor的凸组合性质本身不给这些费用。这里也不把 rank43当成独立性，不给Gaussian迁移或 C_n 新下界。

## 5. 已运行的有限和核对

`probe.py` 用 GF(2⁷)、模0x83、g2、b3、三个basis `[1,2,4]`、`[3,8,16]`、`[5,32,64]`。三尺度10/9、7/6、11/9、p1/2、ε1/128、M4；相位1/7与4/9作共同有限混合。两模式分别是所有bank共用 Fourier，以及已知物理平移0、1/9、−1/7。

所有 H 从同一 g^j trace 产生；全部54个端点行只有秩7，绝没有把bank标签独立化。逐域元素Walsh变换实测最大norm1.8888889≤8。每相位原模型枚举127个状态，比较模型枚举127²个 (j_H,j_B)，各自保留全部bank/scale。三个anchor、四个完整η向量、n1/100/10000与两种 mixture构造共96案例全部满足界；最大单anchor n坐标CF差 .0032835773，耗时 .024秒。

这些是 NumPy 浮点有限和核对，不是向外证书；有限相位模型使用其实际 A，未把它当作已认证 Haar quadrature。完整结果与脚本SHA保存在 `RESULTS.json`。复现只需已有 NumPy 与标准库：

```sh
python3 probe.py
```

最小交付为本证明、`probe.py` 和 `RESULTS.json`。没有读取或修改旧轮数据、没有大规模Gaussian重复实验。
