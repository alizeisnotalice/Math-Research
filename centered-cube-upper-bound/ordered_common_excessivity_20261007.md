# Ordered T_r 的共同过剩性：单个正生成元不足，原谱的正 Green 证据

2026-10-07。独立研究；只编辑本稿及同前缀守卫，不改主账。结论：固定 ΣB_i 或固定 Σ(I−G_c,i) 的单一全局正势条件，均不足以推出 T_r u≤u。证据使用原 B/G_c 的真实正 Green 势与空间远场，非 projection 模型、非多 c 共同目标四矩刚性。势不在 L¹，因此不称原有限占用／饱和障碍的反例；不否定 ordered weak endpoint。

## 1. 准确 family、生成元与时钟

固定 c∈(0,1]、物理尺度1，B(v)=v/log(1+v)−1，B_i=B(−∂_i²)，G_c,i=(I+cB_i)⁻¹。原正族为
\[
T_r=\prod_i[(1-r)I+rG_{c,i}]
=\prod_i[I+c(1-r)B_i]\prod_i[I+cB_i]^{-1}.
\tag{1}
\]
对 r≤R<1，其瞬时正衰减生成元是
\[
\mathcal A_r=\sum_i\frac{cB_i}{I+c(1-r)B_i}
=\frac1{1-r}\sum_i(I-G_{c(1-r),i}),
\quad \partial_rT_r=-T_r\mathcal A_r.
\tag{2}
\]
每一项是正概率跳核的 bounded jump 生成元；总率 n/(1−r)。转移 T_sT_t⁻¹ 的真实正表达是原 K 的相应重参数，不把逆 T_t 当正算子。所有符号在共同 Fourier 域上核实，正性由原 G 的正表示给出。

hazard 时钟 z=−log(1−r) 下生成元变为 Σ_i(I−G_{ce^{-z},i})，每轴率1，跳核随 z 收缩。这是非齐次 evolution，不能仅由初始固定生成元的 excessive 类覆盖它。

本轮核对已读 J05 的生成元／共同参数合同，以及原 NB--A/M 与 frozen positive G 表示；不引用未核外部 theorem。下界目录 A/B 的任意联合来源仍保留；本稿的 Gaussian 只是正势反证中的一个明确允许来源，不限制原 weak 输入类。

## 2. 单一条件到底缺哪些正项

记 e_k(B)=Σ_{|A|=k}∏_{i∈A}B_i、G=∏_iG_c,i。准确代数为
\[
\boxed{I-T_r
=G\sum_{k=1}^n c^k[1-(1-r)^k]e_k(B).}
\tag{3}
\]
因此 ΣB_i u=e_1(B)u≥0 只给第一项正；其余 mixed actions 没有自动符号。若这些 actions 是有限 signed measures，得到安全的固定正缺陷上包
\[
\sup_{0≤r≤1}(T_ru-u)_+
≤G\mathcal D_B,
\quad\mathcal D_B=\sum_{k=2}^n c^k[e_k(B)u]_-.
\tag{4}
\]
其 L¹ 费用是 ||D_B||，未证明它由原 source W 控制，不把有限性冒充 source-once √n 预算。所有 e_k(B)u≥0 是一个充分条件，不能由第一项正推出。

对于 A_0=Σ_i a_i、a_i=I−G_c,i，
\[
I-T_r=\sum_{k=1}^n(-1)^{k+1}r^ke_k(a).
\tag{5}
\]
当 A_0u≥0，亦仅得到
\(\sup_r(T_ru-u)_+≤\sum_{k≥2}[(-1)^ke_k(a)u]_+\)，没有 W 质量界。

对足够好的固定 u，最精确的必要充分条件是
\[
u-T_ru=\int_0^rT_v\mathcal A_vu\,dv≥0\quad\text{对每个 }r.
\tag{6}
\]
若要求所有起始参数的转移均保持 u 过剩，则 \(\mathcal A_vu≥0\) 对全部 v 是必要（右导数）且充分（正 evolution/Duhamel）。仅要求从0出发 T_ru≤u，(6) 比逐 v 正更弱；不混淆这两个量词。逐坐标 B_i u≥0 足以给全 v 条件，但一般 fixed obstacle 不提供它。

## 3. 原正 Green 势及严格 Fourier 远场方法

取任意固定 n≥3，正径向 Schwartz 来源
\[
\rho(x)=(4\pi)^{-n/2}e^{-|x|^2/4},\qquad \widehat\rho(\xi)=e^{-|\xi|^2},\quad\int\rho=1,
\]
采用 \(\widehat f=\int e^{-ix\cdot\xi}f\)。定义
\[
L_B=\sum_iB_i,\quad A_0=\sum_i(I-G_{c,i}),
\quad u_B=L_B^{-1}\rho,\quad u_A=A_0^{-1}\rho.
\tag{7}
\]
两个逆均指正 Green 势：先取 \((\varepsilon+L)^{-1}\rho=\int_0^\infty e^{-\varepsilon t}e^{-tL}\rho dt≥0\)，再令 ε↓0。原 B 是 Bernstein function，故 L_B 半群正；A_0 是原 bounded 正跳过程生成元。

近零 symbols 分别是 |ξ|²/2+O(|ξ|⁴)、c|ξ|²/2+O(|ξ|⁴)；n≥3 使 |ξ|⁻² 局部可积，远频 ρhat Schwartz 衰减。因此 ρhat/L 的绝对积分有限，ε 逆以 dominated convergence 一致收敛到有界连续非负 u。确有 L_Bu_B=ρ≥0、A_0u_A=ρ≥0。它们远场分别以2Φ_n、(2/c)Φ_n为首项，Φ_n(x)=|x|^{2−n}/[(n−2)|S^{n−1}|]；故不在 L¹，而 n≥5 时还在 L²。无需把形式 Fourier 导数当矩存在性论证。
非 L¹ 也可不依赖远场直接证明：正 Green 表示与 Tonelli 给 \(\int u=\int_0^\infty\int e^{-tL}\rho\,dxdt=\int_0^\infty1\,dt=\infty\)，保质量成立且点态有限已在上面确认。

本稿的远场余项为真正 O(|x|^{-n-4})，不是数值拟合：用径向 χ=1 于零邻域切分 Fourier。所有 symbols 远离零光滑，乘 ρhat 后每阶导数可积，高频积分按任意次数分部为 rapid decay。近零除去写出的常数、多项式与 homogeneous degree2 项后，余项 R 满足所有 multi-index 的 annular 界 |∂^αR(ξ)|≤C_α|ξ|^{4−|α|}。在 |ξ|≈2⁻ʲ 的 dyadic annulus，积分为
\[
O\big(2^{-(n+4)j}(1+2^{-j}|x|)^{-N}\big),\qquad N>n+4.
\]
按 2⁻ʲ|x| 大小求和得 O(|x|^{-n-4})。余项常数允许依赖固定 n,c,r 和来源，未声称随 n 一致。这些 derivative bounds 来自 B 在0的收敛解析级数、分母正定二次首项，以及有限次乘除；不是从单个 O(|ξ|⁴) 值界推断导数。截断 homogeneous 项与其全分布 inverse 的差在 x≠0 rapid decay（高频高阶导数可积）。多项式部分乘 Schwartz ρhat 也 rapid decay。故以下 principal inverse 是精确的远场首项。

## 4. 两个真实空间方向的失败

写 b_i=B(ξ_i²)、a_i=cb_i/(1+cb_i)。B(v)=v/2−v²/12+O(v³)，直接展开得到
\[
\frac{1-\prod_i(1-ra_i)}{\sum_i b_i}
=rc-c^2r(1-r/2)\frac{\sum_i b_i^2}{\sum_i b_i}
 -(r^2c^2/2)\sum_i b_i+O(|ξ|^4).
\tag{8}
\]
其中非局部 degree2 主项为 \(-c^2r(1-r/2)\sum_iξ_i^4/(2|ξ|^2)\)。对 A_0 则
\[
\frac{1-\prod_i(1-ra_i)}{\sum_i a_i}
=r-(r^2c/2)\frac{\sum_{i<j}ξ_i^2ξ_j^2}{|ξ|^2}+O(|ξ|^4).
\tag{9}
\]
多项式／常数项不贡献远场。结合 §3 的 remainder，
\[
\begin{split}
(I-T_r)u_B(x)&=-\frac{c^2r(1-r/2)}2\sum_i\partial_i^4Φ_n(x)+O(|x|^{-n-4}),\\
(I-T_r)u_A(x)&=\frac{r^2c}{4}\sum_i\partial_i^4Φ_n(x)+O(|x|^{-n-4}).
\end{split}
\tag{10}
\]
第二式使用 Δ²Φ_n=0 于 x≠0，故 mixed fourth derivatives 的和为 \(-\tfrac12\sum_i\partial_i^4Φ_n\)。

对一般 n≥3 的 Laplace Green 常数 C_n>0，直接微分给
\[
\sum_i\partial_i^4(C_n|x|^{2-n})
=C_n(n-2)n(n+4)|x|^{-n-2}
\left[(n+2)\sum_i\eta_i^4-3\right],\quad\eta=x/|x|.
\]
对任意 n≥3，axis 的角向括号为 n−1>0，平衡对角为 −2+2/n<0。以 n=3 为具体常数例：远轴 x=Re_1 上该值为 \(21/(2\pi)R^{-5}>0\)；平衡对角 x=R(1,1,1)/√3 上为 \(-7/\pi R^{-5}<0\)。对任意固定 c>0、r∈(0,1)，余项是更低阶，因此
\[
\boxed{T_ru_B(Re_1)>u_B(Re_1),\qquad
T_ru_A(R(1,1,1)/\sqrt3)>u_A(R(1,1,1)/\sqrt3)}
\tag{11}
\]
对全部足够大 R 成立。同样的两种严格符号对任意固定 n≥3 都成立，n≥5 时上述正势还是 L²（仍非 L¹）。两者分别否定题设两种 single-generator positive condition 的普遍迁移；均使用原实际谱、同一个固定 c、完整 T_r。

## 5. 对原障碍与 weak 目标的边界

非零正全局 Green 势有无限 L¹ 质量。反而，若 conservative 生成元的 L¹ 域中 u≥0、Lu≥0 且 ∫Lu=0，则 Lu=0；对本原谱，其 Fourier 仅在零点可能为零，L¹ Fourier 连续性迫使 u=0。因此不能把“全空间 Lu≥0”当作原有限占用障碍的现成条件。原障碍是在 Ω 上有 killed operator、饱和 μ、外域正跳项；本稿没有证明原饱和 u 违反某个 actual 门，也不强行以不加证明的有限域 exhaustion 认证它。

本稿 (11) 已足以排除普遍 positive-potential transfer。要用固定障碍控制完整 ordered family，必须另给 (6) 或 all-future 正条件／mixed defect 的真实 source-once 预算。这里的 D_B/D_A 仅准确量化所欠，不声称其质量≤Wpolylog；原 weak 端点仍可能由其它机制成立。主 cube R_dagger、完整 FIRST fullfuture 与全部 actual history 合同不变。

## 6. 专属三轮守卫预登记

只准备针对新 symbol/空间主项的 deterministic rational guard：三轮 n=8/32/128，各取 c=1/2,3/4,1 及 r=0,1/4,1/2,3/4,1。从 log(1+v)/v 直接有理反演生成 B 系数，不硬编码 principal sign；核 (8)(9)、真实 polynomial identity (1)(3)(5)、轴与平衡对角 fourth-derivative sign，及实际 B 的小 Fourier 点有理 alternating-log interval。来源为上述固定正 Gaussian，未拟合 source。

这些只认证解析式的组成和原 symbol 的有限点外包，不能给出 (11) 的有限 R 门槛，不能认证一般过剩性、原 FIRST/history、有限 L¹ 饱和障碍或弱范数。连续远场由 §3–4 解析负责；数值不替代余项证明，不重跑旧 projection、共同四矩或 rational-mesh 守卫。

终态：新 `ordered_common_excessivity_exact_guard_20261007.py` 与 `_results.json` 已保存。三轮 n=8/32/128 分别1560/3720/12360项，共17641项 Fraction 精确谓词 PASS；每轮45个 formal symbol cases、45个原 B 小 Fourier 点外包（ξ=εv，ε=1/(32n)），另核 axis/diagonal 归一化空间主项的符号。r=0 的零式与 r=1 的 T endpoint 也保留，不把 r=1 的瞬时 unbounded 生成元当 bounded 参数。

代表性 c=1,r=1/2、Fourier v=(1,…,1) 的 degree2 系数：n8 的 D_B/D_A 为 −11/16、−7/16；n32 为 −35/16、−31/16；n128 为 −131/16、−127/16。这些 Fourier 系数不等于相同空间方向的 Green 首项符号，后者由 (10) 的 inverse 和角向微分得到，守卫单列核查，没有混淆两种方向。

首次运行在大整数转十进制 hash 时触及 Python 字符串位数限制，未保存成功结果、没有数学谓词失败。改为有符号及长度前缀的 binary SHA 后，仅重启该新脚本一次，全部完成；未修改阈值/tolerance，不重跑其它实验。无随机种子，无 live handle；JSON 不储存巨大有理展开，只存主项小分数与外包 hash。原 B 的有限点 tolerance 是预登记 symbol probe，不是均匀余项或空间有限 R 证书。

最终可复用结论：两类 single-generator positivity 的普遍过剩性转移均被原正 Green 势否定；必要补项由 (3)–(6) 精确列出。主 cube/general weak 目标未由此否定或闭合，原实际 R_dagger 不变。
