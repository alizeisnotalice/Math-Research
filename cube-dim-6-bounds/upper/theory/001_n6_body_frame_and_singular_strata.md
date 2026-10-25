# 第001轮：六维 body Cauchy frame 的通用上界接口与奇异层

日期：2026-10-25。本文是六维上界的首轮理论输入，不是 `C_6`
新上界证书。所有积分取 Lebesgue 测度；方体均为中心、轴平行且六个
坐标共用一个半径。

## 模型及有限准确诊断

取一个不定向 body 方向的代表集

\[
 V=\{(1,s_2,\ldots,s_6):s_j\in\{-1,1\}\},\qquad |V|=32,
\]

并令每个方向独立地配尺度为 \(t\) 的一维 Cauchy 变量。于是
\(\psi(\xi)=\sum_{v\in V}|v\cdot\xi|\)，\(L=-\psi(D)\)，且
\(P_t\) 是真实的保正、保质量、自伴卷积半群。由于

\[
 \sum_{v\in V}vv^T=32I_6,
 \qquad \sqrt {32}|\xi|\le\psi(\xi)\le32|\xi|,
\]

其 \(L^2\) 生成元域为 \(H^1(\mathbb R^6)\)，形式域为
\(H^{1/2}(\mathbb R^6)\)。

`body_frame_rank_audit.py` 实际枚举了 \(\binom{32}{6}=906192\)
个六行候选作为浮点筛查；其最高圆整数候选值是160，但该扫描**不是**
全局最大行列式证明。用于以下结论的准确证书只是编号
\([0,3,5,14,22,25]\) 的六行：SymPy Bareiss 行列式精确等于 \(160\)。
全部 32 个单方向、496 个二元组及4960 个三元组分别具有秩1、2、3。
完整运行输出在 `body_frame_rank_audit.json`。

取这个六列基 \(B\)，其余 26 项只作概率卷积，故真实密度满足

\[
 0<p_t(x)\le {1\over160\pi^6}t^{-6},\qquad
 p_t(x)=t^{-6}p_1(x/t). \tag{1}
\]

这是一个严格可用的热核输入；它没有认证 \(G=\int_0^\infty p_tdt\)
在方体上的下界或 face 积分。

对 \(\phi\in C_c^\infty\)，生成元的二差分表示中的每个方向都满足
\[
 \|\phi(\cdot+sv)+\phi(\cdot-sv)-2\phi\|_{L^r}
 \le\min\{s^2|v|^2\|D^2\phi\|_{L^r},4\|\phi\|_{L^r}\},\quad r=1,\infty.
\]
所以 \(L\phi\in L^1\cap L^\infty\cap C_0\)。这既使后文的全空间源配对
合法，也给固定大 cutoff 的 \(\|L\chi_R\|_\infty=O(R^{-1})\)。另外，
\(\psi\widehat\phi\in L^1\) 给 \((P_t\phi-\phi)/t\to L\phi\) 于 \(L^\infty\)，
是正部源测试极限所需的具体收敛，而不是未说明的点态生成元断言。

## 六维的真实奇异集合

令
\[
 S_2=\bigcup_{v,w\in V}\operatorname{span}\{v,w\}.
\]
这里包含方向线（令两个方向相同即可）。对本模型有

\[
 G(x)<\infty\ \text{且连续于 }\mathbb R^6\setminus S_2,
 \qquad G(x)=+\infty\ \text{于 }S_2. \tag{2}
\]

证明要点如下，且与四维的“线外”论证不同。设 \(K\) 是与 \(S_2\)
距离正的紧集。取足够小的 \(\delta>0\)。若至多两个 Cauchy 坐标的绝对值
超过 \(\delta\)，其余 30 个向量和与某个二平面之间的距离由
\(30\sqrt6\delta\) 控制，不能落入 \(K\)。若大坐标的方向仍只张成一个
二平面，同样矛盾。因此落入 \(K\) 时可挑出三个线性无关的大坐标。
有限诊断保证任何三个不同 body 方向都线性无关；再贪心补成六向量基。
这三个受限 Cauchy 密度各至多 \(t/(\pi\delta^2)\)，补基的其余三个
密度各至多 \((\pi t)^{-1}\)。对有限三元组及补基求和，余变量只作概率
卷积，给

\[
 \sup_{x\in K,0<t\le1}p_t(x)<\infty. \tag{3}
\]

上述事件界的密度依据：在K的一个与S2仍有正距离的开邻域U内，完整参数空间的指示函数1_{X∈U}不超过所有“三指定参数大于δ”事件指示函数之和。对每个三元组分别选补基B，对六个基参数作可逆线性换元，其Jacobian为|det B|。固定余参数后，该受限law的密度上界是三个受限Cauchy sup与三个未受限sup之积除以|det B|；对余参数Tonelli积分，其总概率为1。有限事件求和得到U内a.e.密度上界，各事件可以使用不同补基。真实p_t连续，因此上界成立于U每一点。这不把零概率的条件X=x当作事件概率，也不把不同基坐标误当同一独立参数组。

结合 (1) 的长时可积支配，DCT 给出 (2) 的第一部分及连续性。

反过来，若 \(x=av+bw\) 且 \(v,w\) 独立，选包含它们的 body 基。
令其他变量的缩放和落在该基坐标 \(O(t)\) 盒内；该事件的概率与 \(t\)
无关且严格为正。两条指定坐标若非零各贡献 \(\asymp t\)，其余四条
在零附近贡献 \(\asymp t^{-1}\)，故 \(p_t(x)\gtrsim t^{-2}\)。某系数为
零只会提高发散阶，原点由缩放和正性同样发散。这证明 (2)。

所以六维的高维断点已被明确定位：奇异层是有限个二维子空间，而不是
四维的方向线。它们虽是零六维测度，不能从逐点 Green 公式或 face 采样中
忽略；但它们不妨碍下面的全空间 \(L^1\) 正部接口。

## 可从共享分析重新建立的六维接口

下列推导只使用有限正权 spanning Cauchy frame、\(n=6>1\)、(1) 和
Lebesgue 卷积。它是对资料包 `paper_finite_source.md` 的逐维重检，不把
三、四维接受标签迁移成六维结论。

对任何非负、有界、紧支撑 \(f\) 及 \(\kappa>0\)，离散障碍逼近给
\(f=\sigma-Lu\)，其中

\[
 u\in L^1\cap L^\infty\cap H^1,\quad 0\le\sigma\le\kappa,
 \quad\sigma=\kappa\ \text{a.e. on }D=\{u>0\},\quad
 \int\sigma=\int f. \tag{4}
\]

真正需要重查的两个指数均在六维成立：
\(\|P_sf\|_\infty\le\min(\|f\|_\infty,C\|f\|_1s^{-6})\) 使离散势的
右 Riemann 和有界；频率尾仍以 \(\psi_t\ge(1-e^{-1})Q\) 控制，局部
Rellich 给极限。质量回收仍由 \(\|L\chi_R\|_\infty=O(R^{-1})\) 及
\(u\in L^1\) 完成。这里没有使用 (2) 的线外连续性。

令 \(H_c=(G-c)_+\)。由 (1)，短长时拆分给 \(G\in L^1_{\rm loc}\) 及
\(H_c\in L^1\)，更精确地 \(G\in L^{6/5,\infty}\)。半群分布计算和
正最小原则给全空间

\[
 LH_c=\mu_c-\delta_0,\qquad \mu_c\ge0,\quad
 \mu_c(\mathbb R^6)=1. \tag{5}
\]

小 cutoff 估计现在为 \(\|L\chi_\epsilon\|_1=O(\epsilon^5)\)，故
\(\mu_c(\{0\})=0\)。这一步直接处理所有二维奇异平面及其交叠，未从
它们的补集延拓任何点值结论。

若后来认证 \(G\ge g_0>0\) a.e. 于 \([-1,1]^6\)，取 \(c=g_0/6\) 并令

\[
 K={H_c\over 64(g_0-c)}.
\]

则 (5) 与 (4) 的有限源卷积比较在同一个 \(D\) 外同时控制所有 \(r>0\)，
先由正有理尺度、再由 \(L^1\) dilation 连续性扩展到实尺度，最后以单调
截断扩展至一般 \(L^1\) 输入。因此这是满足题设量词的条件性上界接口。

令 \(q=6/5\)，\(I=\int_{[0,1]^5}G(1,y)^qdy\)。六维 face/radial
Jacobian 计算给条件性成本公式

\[
 C_6\le6^{6/5}{I\over g_0^{6/5}}. \tag{6}
\]

它需要一个连续或 essential 的 cube 下界、一个严格 face 积分上界，以及
二者所对应同一真实 \(G\)。本轮没有这些量，故 (6) 不是数值上界、不是
对 \(C_6\) 的改进。

## 六维接口的逐式证明及 whole-face 短时预算

这一节补出前节所用的全部分析步骤。设 \(C=(160\pi^6)^{-1}\)，
\(F=\|f\|_\infty\)，\(M=\|f\|_1\)。对 \(0<t\le1\)，定义

\[
 V_t=t\sum_{j\ge0}P_{jt}f,\qquad
 u_{t,0}=0,\qquad u_{t,k+1}=(P_tu_{t,k}+t(f-\kappa))_+.
\]

由 (1) 和 Young 不等式，
\(\|P_sf\|_\infty\le\min(F,CMs^{-6})\)。递减右 Riemann 和给

\[
 \|V_t\|_\infty\le F+{6\over5}(CM)^{1/6}F^{5/6}=:A_f. \tag{7}
\]

固定 \(t\) 时 \(V_t-tf\in C_0\) 且 \(P_tV_t=V_t-tf\)。因此迭代递增，
受 \(V_t\) 控制，并收敛至最小固定点 \(u_t\)，且
\(0\le u_t\le(V_t-\kappa t)_+\)。后者在无穷远趋零，故每个固定 \(t\)
的 \(u_t\) 有本质紧支撑、属于 \(L^1\cap L^2\)；此处没有要求 \(V_t\)
本身可积或这些支撑一致。

令 \(L_t=(P_t-I)/t\)、\(\sigma_t=f+L_tu_t\)。固定点在正支与零支
分别给出 \(\sigma_t=\kappa\) 和 \(0\le\sigma_t\le\kappa\)。因为 \(u_t\)
可积，半群质量守恒给 \(\int\sigma_t=M\)。将等式和 \(u_t\) 配对，得到

\[
 \langle u_t,-L_tu_t\rangle+\kappa\|u_t\|_1=\int fu_t.
\]

故 \(\|u_t\|_1\le MA_f/\kappa\)、\(\|u_t\|_2^2\le MA_f^2/\kappa\)、
\(\|\sigma_t\|_2^2\le\kappa M\)，并且
\(\|L_tu_t\|_2\le\|f\|_2+\sqrt{\kappa M}\)。

写 \(\psi_t=(1-e^{-t\psi})/t\)。在 \(t\le Q^{-1}\)、\(\psi\ge Q\) 时，
\(\psi_t\ge(1-e^{-1})Q\)，所以 \(\widehat u_t\) 的该频率尾有统一
\(O(Q^{-2})\) 的 \(L^2\) 质量。在 \(\psi\le Q\) 时，\(\psi_t\) 与
\(\psi\) 一致可比。频率截断、局部 Rellich 和尾估计给某个 \(t\downarrow0\)
子列上 \(u_t\to u\) 于 \(L^2_{\rm loc}\) 和 a.e.，同时 \(u_t,\sigma_t\)
弱 \(L^2\) 收敛。Fatou 保留 \(L^1\)、\(L^\infty\) 界；
\(L_t\phi\to L\phi\) 于 \(L^2\) 传递分解，得到

\[
 f=\sigma-Lu,\quad u\in L^1\cap L^\infty\cap H^1,\quad
 0\le\sigma\le\kappa. \tag{8}
\]

在 \(u>0\) 的 a.e. 点最终 \(u_t>0\)，故在任意有界测试域上由 DCT 和弱
收敛有 \(\sigma=\kappa\)。先以正大截断给 \(\int\sigma\le M\)，再由
\(\|L\chi_R\|_\infty=O(R^{-1})\)、\(u\in L^1\) 在 (8) 中令 \(R\to\infty\)，
得到 \(\int\sigma=M\)。这证明 (4)，特别是 \(\kappa|D|\le M\)。

为构造正部源，令 \(U_T=\int_0^Tp_tdt\)、\(W_T=\int_T^\infty p_tdt\)。
有 \(\|U_T\|_1=T\)、\(\|W_T\|_\infty\le CT^{-5}/5\)。选择
\(T=(2C/(5\lambda))^{1/5}\)，给

\[
 |\{G>\lambda\}|\le2(2C/5)^{1/5}\lambda^{-6/5}. \tag{9}
\]

另一方面 \(\int_EG\le\int_0^\infty\min(1,C|E|t^{-6})dt<\infty\)。因此
\(G\in L^1_{\rm loc}\)、\(H_c\in L^1\)，包括所有 \(S_2\) 邻域及远尾。
对 \(\phi\in C_c^\infty\)，同一拆分给
\(\int G|L\phi|<\infty\)，半群导数给 \(LG=-\delta_0\)。

Tonelli 给 \(P_tG=\int_t^\infty p_sds\)，从而 a.e. 有 \(G=P_tG+U_t\)。
正部的 1-Lipschitz 性给 \(P_tH_c\ge H_c-U_t\)。将该式与非负测试函数
配对、除以 \(t\) 并令 \(t\downarrow0\)，得到 \(LH_c+\delta_0\ge0\)。
故 \(\mu_c=LH_c+\delta_0\) 是正 Radon 测度。以 \(\chi_R\uparrow1\) 测试：

\[
 \int\chi_Rd\mu_c=1+\int H_cL\chi_R\longrightarrow1,
\]

故其总质量为一。令 \(h=\min(G,c)\)，则 \(\mu_c=-Lh\)。小 cutoff 的
\(L^1\) 缩放是 \(\|L\chi_\epsilon\|_1=\epsilon^5\|L\chi\|_1\)，故
\(\mu_c(\{0\})=0\)。这完成 (5)，没有把平面奇异性排除出分布论证。

设以后已认证 \(G\ge g_0\) 于 cube，\(c=g_0/6\)，并令
\(N=64(g_0-c)\)、\(K=H_c/N\)。Fourier 变换严格给
\(-\psi\widehat K=N^{-1}(\widehat\mu_c-1)\)。所以对 (8) 的 \(u\)，

\[
 K_r*Lu={1\over Nr}(\mu_{c,r}*u-u),
 \qquad K_r(x)=r^{-6}K(x/r). \tag{10}
\]

在同一个 \(D^c\) 上，右端非负（有限正测度卷积与 \(u\ge0\) 的 Fubini
结论）。取 \(\kappa=\lambda/\int K\)，先合并正有理尺度的零集，再由
\(L^1\) dilation 连续性推广到全体 \(r>0\)，得坏集包含于 \(D\)。单调
截断 \(f_j\uparrow|f|\) 推广到一般 \(L^1\) 输入。这给出题设的全尺度、
一般输入比较，而不是只给有限尺度或光滑输入。

最后，齐次性 \(G(ry)=r^{-5}G(y)\) 和每 face 的径向 Jacobian \(r^5\)
给 \(\int H_c=(5/6)c^{-1/5}J\)。cube 对称性给
\(J=6\cdot64I\)。除以 \(N\)、代入 \(c=g_0/6\)，恰得 (6)。

### whole-face 短时间 \(L^{6/5}\) 预算

这条估计替代逐一二维奇异平面的帽覆盖。先取开 slab
\(x_1>1-\eta\)（\(0<\eta<1\)）。32 个 body 变量的第一坐标系数都为一，
所以在该 slab 的每一个参数 fiber 上，至少一个参数绝对值不小于
\((1-\eta)/32\)。已精确认证的 determinant-160 基包含代表方向；坐标符号
群在 body 方向线上传递，故每个 active 方向均可补成一个绝对行列式160的基。
把该大参数与其余五个基参数取密度上界，再对32个可能方向求和，先由
Tonelli 得到 a.e. 的 slab 密度界
\[
 p_t(x)\le {32^3\over160\pi^6(1-\eta)^2}t^{-4}.
\]
真实密度 \(p_t\) 连续，故此 a.e. 界延拓到整个 slab；令 \(\eta\downarrow0\)
才得到 slice 上的

\[
 \|p_t(1,\cdot)\|_\infty\le C_f t^{-4},\qquad
 C_f={32^3\over160\pi^6}={1024\over5\pi^6}. \tag{11}
\]

另一方面第一坐标边际的分布是尺度 \(32t\) 的 Cauchy。Tonelli 先给几乎
处处的边际等式。为得到固定 slice 的充分上界，对 \(s\to1\) 的区间平均
使用该等式，并对连续非负 \(p_t(s,y)\) 用 Fatou，得到

\[
 \|p_t(1,\cdot)\|_1\le {32t\over\pi(1+32^2t^2)}\le {32t\over\pi}. \tag{12}
\]

对 \(q=6/5\) 插值 (11)--(12)：

\[
 \|p_t(1,\cdot)\|_{6/5}
 \le C_f^{1/6}(32/\pi)^{5/6}t^{1/6}. \tag{13}
\]

body frame 对后五坐标逐个翻转对称；正正交象限的 \(L^{6/5}\) 范数为全
face 的 \(32^{-5/6}\)。令 \(T_\epsilon=\int_0^\epsilon p_tdt\)，Minkowski
和 (13) 因此给可直接用于连续证书的严格预算

\[
 \|T_\epsilon(1,\cdot)\|_{L^{6/5}([0,1]^5)}
 \le {6\over7}32^{-5/6}C_f^{1/6}(32/\pi)^{5/6}\epsilon^{7/6}
 ={6\,2^{5/3}\over7\,5^{1/6}\pi^{11/6}}\epsilon^{7/6}. \tag{14}
\]

定义截时势 H_ε(x)=∫_ε^∞p_t(x)dt；Tonelli给G=H_ε+T_ε，不作∞−∞。所以若未来的有限连续程序严格证明 \(\|H_\epsilon(1,\cdot)\|_{6/5}\le B\)，
则 \(I^{5/6}\le B+\)右端(14)。本轮没有运行该五维连续认证，故 (14)
只消除了一个分析缺口，不产生 \(I\) 或 \(C_6\) 的数值界。

## 首批验收结论与下一步

本批完成的严格结论是：body frame 的 Gram、一个 \(|\det|=160\) 基、(1)、六维
奇异层 (2)，以及在该模型下可重新闭合的条件性一般 \(L^1\)/全尺度接口。
它排除了把四维线奇异诊断照搬到六维的路线。

下一批应在这个固定 body frame 上优先认证截时 \(H_\epsilon\) 的五维连续
下界与 \(L^{6/5}\) face 积分，再用 (14) 一次性补回短时间项。该 whole-face
预算不要求逐个 \(S_2\) 帽；它把奇异平面及其交叠纳入全空间 \(L^{6/5}\)
控制。只有连续证书的覆盖、插值/区间误差和 (14) 同时闭合，才能形成
\(g_0,I\) 并进入成本比较。

## 来源与复现

输入来源的 SHA-256：

* `papers/n3/main.tex`: `cf87e51a19d212cbdde601968352f5b9a4ea1b2e2d5d10ac908b7072fd80715c`。
* `papers/n4/main.tex`: `1e173ee28eed26945f5f94523c902d4c1496f9218cc3fd488376ee01a76b984b`。
* `audit/second_review/theory/paper_finite_source.md`: `5e5cd1cac6cd3172eeff584696661066c61c3b7a666ed287bc1cbd5e93325608`。
* `audit/second_review/theory/audit.md`: `f683a1f7e7e1b01e20063a93e4488102e690c8b73ed1e9865707e83f93b849c6`。

有限诊断的实际命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/upper/theory/body_frame_rank_audit.py
```

实际环境为 Python 3.15.0、SymPy 1.14.0、python-flint 0.9.0、NumPy
2.5.4、SciPy 1.18.1；该批未调用 MATLAB。
