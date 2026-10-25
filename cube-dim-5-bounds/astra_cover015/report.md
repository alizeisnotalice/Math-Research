# astra_cover015：32-clique 独立审计、共享包障碍与全局源截减

## 结论与范围

固定 Lebesgue 测度，研究
\[
 M f(x)=\sup_{r>0}(2r)^{-5}\int_{x+[-r,r]^5}|f(y)|\,dy.
\]
所有坐标始终使用同一个半径；解析结论适用于一般 \(L^1\) 输入与全部正尺度。
本轮**没有证明 \(C_5<32\)**。成果是四项可独立核验的结论：

1. 独立确认冻结的 `upper_cover014` 中任意不等半径 32-clique 的 core 几何与局部 \(127/4\) 估计；退化及零体积边界没有发现缺口。
2. 构造一个新的 993-cube 家族：32 个满重叠 cliques 全部共用同一个中央 cube；所有见证平均可同时接近阈值，全部源恰为32重，见证并集/质量趋于32。但 014 的“索引不交完整 clique + 剩余重叠数”收费，在这一固定合法家族上的精确最优值仍是32。
3. 证明全局源截减：只移除至多 \(D_{\rm avg}+D_{\rm mult}\) 的质量，就能同时得到“源支撑在32重区”与全部 cube 的质量上限。这是一次性全局处理，适用于所有有限家族，不重复对每个包扣除同一质量。
4. 给出一个明确、尚待证明的 capped-core 稳定估计，并证明它一旦成立将推出显式全局改进；包括从有限覆盖到一般 \(L^1\)、全部尺度的完整传递。

没有将局部 whole-mass 界、见证并集、有限 LP 最小值或候选稳定估计记为新的 \(C_5\) 上界。

## 1. 来源与独立审计

实际核对了以下冻结 SHA-256：

| 文件 | SHA-256 |
|---|---|
| `upper_cover013/report.md` | `f787895653d2d382a28321ed44e0191e31a0ce47e93710dffb8c8cbc85d0a91d` |
| `upper_cover014/report.md` | `e863d180d9d125b980b9e9c269d2611d1f4b6fd45222e1171efde1cd97309f67` |
| `upper_cover014/registration.json` | `41007eb61d05e9f375347d0137aae91a56696e827103f330252b7e4bb881cd62` |
| `upper_cover014/verify.py` | `e633ba51410941ed97f0db49959da6d28a5049b344fba054823132df942032c6` |
| `upper_cover014/results.json` | `3df5428219086ff73927f61ebeae4d58c69fab7a27d5a54b8533b9fcacf355bb` |

直接阅读了 [Aldaz v2 原文](https://arxiv.org/html/1811.10478v2) 的 Definition 1.1 与 uniform/fixed-background 区分。此处只沿用其 \(N\le32\)、弱型32基准；其对所有背景测度的尖锐性不等于固定 Lebesgue 尖锐性。两条定向检索没有新增可认证的五维严格上界来源；不声称穷尽文献。

### 1.1 闭正交象限的边界

设32个闭 cubes \(Q_i=B_\infty(c_i,r_i)\) 共点，且
\(\|c_i-c_k\|_\infty>\max(r_i,r_k)\) 对所有 \(i\ne k\)。将共同点移到0。
同一闭正交象限内两个中心在每个坐标的距离都不超过两绝对值的较大者，故不可能共处该象限。
每个中心至少占据一个闭象限；若有零坐标，就占据至少两个。32个互不重复的中心恰好耗尽32个闭象限，因此没有零坐标，每个开象限恰有一个中心。

对任一坐标 \(j\)，取符号只在 \(j\) 不同的一对中心。其余坐标同号，不可能负责严格中心排除。因此
\[
 |c_{i,j}-c_{k,j}|>\max(r_i,r_k).
\]
这也保证两个一维区间互不包含；它们的交长恰为
\[
 r_i+r_k-|c_{i,j}-c_{k,j}|<\min(r_i,r_k)\le R,
 \qquad R=\max_i r_i.
\]
共同 core \(K\) 的每个坐标宽度均小于 \(R\)，故能放入半径 \(R/2\) 的 cube。
这里既不要求半径可比，也不要求共同点位于 core 内部。若 core 零体积，则 \(L^1\) 源在其上的质量为零。

### 1.2 局部界的各量词

记 \(m=\|f\|_1\)、\(t=\|f1_{K^c}\|_1/m\)、\(U=\bigcup_iQ_i\)，并假定全部见证平均严格大于 \(\lambda\)。最大 cube 给 \(m>\lambda(2R)^5\)，且
\[
 \lambda|U|<\int f\sum_i1_{Q_i}\le(32-t)m.
\]
当 \(t\ge1/4\) 时得 \(\lambda|E_\lambda\cap U|\le127m/4\)。
当 \(t<1/4\) 时，按阈值 \(\lambda/2\) 分解 core 源和外源；core 源的所有好半径都由总质量上界约束，其 halo 体积满足
\[
 \frac\lambda m|\{M(f1_K)>\lambda/2\}|
 \le\left[\tfrac12(\lambda(2R)^5/m)^{1/5}+(2(1-t))^{1/5}\right]^5
 <(17/10)^5.
\]
外源使用已证32基准，贡献小于 \(64t m<16m\)。准确比较
\((17/10)^5+16<127/4\) 完成证明。
这条证明使用全输入质量，不能逐包直接相加；这一限制没有因不等半径推广而消失。

## 2. 一个共同中央 cube 的32个满重叠 cliques

本节为新的解析族，计算只在一个预登记有理参数上执行。
取 \(2/3<a<1\)，中央 cube 为
\[
 Q_0=[-1,1]^5.
\]
另取所有半径 \(1/2\) 的 cubes，其中心属于
\[
 \{\pm a/2,\pm3a/2\}^5\setminus\{\pm a/2\}^5.
\]
共有 \(4^5-2^5=992\) 个小 cubes，加中央 cube 共993个。

### 2.1 全族确为 Besicovitch，且每个32-clique都含中央 cube

不同小中心至少有一坐标相距 \(a>1/2\)。任一小中心至少有一坐标绝对值 \(3a/2>1\)，所以中央中心与它也严格互相排除。

在任一点，每个坐标最多被四个一维小区间中的两个包含，因为 \(a>1/2\)。若全部五个坐标各有两个候选，这些候选必须是相邻中心对，每一对都含至少一个内层中心 \(\pm a/2\)。因此至少一个全内层组合已被删除，小 cubes 的重叠至多31。
故所有32-cliques都必须含 \(Q_0\)。要达到31个小 cubes，五个坐标的候选对各自必须恰有一个内层和一个外层中心；于是恰有32个满 cliques，分别对应 \(\sigma\in\{-1,1\}^5\)。

令
\[
 p_\sigma=a\sigma,\qquad h=(1-a)/2,
 \qquad K_\sigma=p_\sigma+[-h,h]^5.
\]
这32个盒子正是32-clique cores：中央 cube 与相应31个小 cubes 的共同交。它们两两不交且位于中央 cube 内。

在每个 core 上均匀放质量 \(1+\eta\)，\(\eta>0\)，得有界、紧支撑、绝对连续输入。每个小 cube 只含其对应 core，中央 cube 含所有32个 cores。因此
\[
 m=32(1+\eta),\qquad
 \int_{Q_i} f=(1+\eta)|Q_i|\quad\text{对全部993个 cubes},
\]
以及
\[
 D_{\rm mult}=0,
 \qquad \frac{D_{\rm avg}}m=\frac{32\eta}{1+\eta}.
\]
取 \(\eta=0\) 时恰好是全部质量上限饱和的 critical 输入；正 \(\eta\) 保证严格好见证。

### 2.2 见证并集趋于后段32，仍不等于接收域

若补回删除的32个内层小 cubes，四个一维区间连续覆盖一个区间，完整网格并集就是
\([-(1+3a)/2,(1+3a)/2]^5\)。被删除内层 cubes 的并集含于中央 cube，故补回它们不改变原993个 cubes 的并集。因此
\[
 U=[-(1+3a)/2,(1+3a)/2]^5,
 \quad \sum_i|Q_i|=1024,
 \quad \frac{|U|}m=\frac{(1+3a)^5}{32(1+\eta)}\longrightarrow32.
\]
这使后段体积重叠、平均余量与源重复计数都接近等号。

实际接收域有可量化缺口。对计算中的
\(a=65535/65536,\eta=2^{-24}\)，考虑32个外角盒
\[
 T_\sigma=\{x:3/2\le\sigma_jx_j\le(1+3a)/2\text{ for every }j\}.
\]
若某个半径在阈值1上好，由总质量知 \(r<(1+\eta)^{1/5}<101/100\)。对 \(x\in T_\sigma\)，这个半径不能碰到任何 \(K_\tau\), \(\tau\ne\sigma\)。若只碰到匹配 core，其质量又强迫
\(r<(1+\eta)^{1/5}/2<(1+\eta)/2\)，但距离至少为 \((2-a)/2>(1+\eta)/2\)，矛盾。
故全部 \(T_\sigma\) 都在真实 \(E_1\) 外，并且
\[
 \sum_\sigma|T_\sigma|=(3a-2)^5,
 \qquad \frac{(3a-2)^5}{m}>\frac1{33}.
\]
对任何 \(A\subset E_1\cap U\)，其 \(D_{\rm rec}/m>1/33\)。因此本构造**没有**使三个真实缺陷同时趋零；它清楚保留了 receiver 约束。

### 2.3 014 的完整包收费在这一固定家族上恰好不能改进32

014 的收费形式为
\[
 \int f\left[\frac{127}{4}\sum_P1_{D_P}+N_0\right],
 \qquad D_P=\bigcup_{i\in P}2Q_i,
\]
其中被选完整32-cliques的 cube 索引互不重复，\(N_0\) 来自剩余 cubes。
因为所有满 cliques 共用 \(Q_0\)，最多只能选一个。

若不选，收费恰为 \(32m\)。若选任意一个，\(2Q_0\subset D_P\) 包含所有源，所以包收费为 \(127m/4\)。被选 core 的剩余重叠为0；其他31个 cores 的剩余重叠均为31。因此总收费/质量恰为
\[
 \frac{127}{4}+\frac{31\cdot31}{32}=\frac{1977}{32}>32.
\]
33种合法选项（不选或选32个中的一个）的精确最小值为32。
即使允许这些选项的随机混合/凸组合，也不能降低这一数值。

这是对**这一固定合法家族和这一特定收费形式**的准确排除。
没有证明所有其他覆盖选择也失败；没有排除分割源、部分包、不同接收分配或利用 \(D_{\rm rec}\) 的改进。

## 3. 小亏损仍不能直接控制膨胀源收费

另一新族区分了“先截减小质量”和“直接给整源计费”。取
\[
 r_j=q^j,\ p_j=(br_j,0,0,0,0),\ c_{j,\sigma}=p_j+ar_j\sigma,
 \quad j=0,\ldots,J,
\]
其中 \(q=1/8,b=11/4,a\) 充分接近1。
各包有32个等半径 cubes；包并集在第一坐标分离。膨胀包
\(D_j=p_j+[-(2+a)r_j,(2+a)r_j]^5\) 嵌套，且恰好含所有 \(k\ge j\) 的 cores。

记 \(v_j=(2r_j)^5,S=\sum_jv_j\)。在 core \(j\) 放质量
\[
 \mu_j=(1+\eta)v_j+\varepsilon S1_{j=J}.
\]
则全族合法、全部源32重，且
\[
 \frac{D_{\rm avg}}m=\frac{32(\eta+\varepsilon)}{1+\eta+\varepsilon},
 \qquad \frac{|U|}m=\frac{(1+a)^5}{1+\eta+\varepsilon},
\]
而
\[
 \frac{\sum_j\int_{D_j}f}{m}
 =\frac{\sum_j(j+1)\mu_j}{m}
 \ge\frac{\varepsilon(J+1)}{1+\eta+\varepsilon}.
\]
因此可同时使前两项分别趋于0与32，并使 raw 膨胀收费任意大。
这证明仅补上“小 \(D_{\rm avg}\)”也不能直接给整份输入的膨胀收费一个统一界。

计算只取 \(J=64,a=8191/8192,\eta=2^{-24},\varepsilon=2^{-12}\)：
\(D_{\rm avg}/m<1/100\)，\(|U|/m>31.98\)，收费大于 \(128/127\)。未采样真实接收域。

若把每个 core 质量截到 \(v_j\)，则 retained mass 为 \(S\)，删除质量恰为 \(D_{\rm avg}/32\)，且
\[
 \frac{\sum_j(j+1)v_j}{S}<\frac1{1-q^5}=rac{32768}{32767}.
\]
该有限保留质量 LP 的原始可行解 \(z_j=v_j\) 与“每包选择一个 cube 容量约束、对偶权重1”的对偶解目标相等，所以最优性完全精确。
这一特例显示先删除小污染的价值；它不解决第2节的共享 critical cliques。

## 4. 所有有限家族均适用的全局容量截减

**引理（已证）。** 给定有限闭 Besicovitch cube 族，令
\(N=\sum_i1_{Q_i}\le32\)、\(H=\{N=32\}\)，并设每个 \(Q_i\) 对 \(f\ge0\) 的平均大于 \(\lambda\)。记
\[
 D_{\rm avg}=\sum_i\bigl(\int_{Q_i}f-\lambda|Q_i|\bigr),
 \quad D_{\rm mult}=32m-\int Nf,\quad B=D_{\rm avg}+D_{\rm mult}.
\]
则存在 \(0\le g\le f\)，使
\[
 \operatorname{supp}g\subset H\quad\text{a.e.},\qquad
 \int_{Q_i}g\le\lambda|Q_i|\ (\forall i),\qquad
 q:=\|f-g\|_1\le B. \tag{4.1}
\]
若 \(m_g=\|g\|_1\)，其总容量空余
\[
 S_g:=\lambda\sum_i|Q_i|-32m_g
\]
还满足精确恒等式与界
\[
 S_g=32q-B,\qquad 0\le S_g\le31B. \tag{4.2}
\]

**证明。** 先把 \(f\) 限制到 \(H\)，移除质量至多 \(D_{\rm mult}\)，因为 \(32-N\ge1\) 于 \(H^c\)。依次处理各 cube：若当前保留源在 \(Q_i\) 内超出 \(\lambda|Q_i|\)，就在该 cube 内按常数比例降低保留源，恰好达到容量。该步删除量不超过原先的 \(\int_{Q_i}f-\lambda|Q_i|\)。后续仅删除源，故不会破坏此前的容量上限。相加得 (4.1)。

由于保留源支撑在 \(H\)，\(\sum_i\int_{Q_i}g=32m_g\)，所以 \(S_g\ge0\)。另一方面
\(\lambda\sum_i|Q_i|=32m-B\)，代入即得 (4.2)。证毕。

这一算法对非原子 \(L^1\) 输入直接成立；有限原子源也可按质量比例截减。不存在对坐标独立半径、有限源近似或固定尺度的额外依赖。
它不保持原 cubes 对 \(g\) 的严格好见证；任何后续使用局部定理的步骤必须重新核验这一点。

## 5. 一个精确的剩余目标与其完整传递

以下 **capped-core 稳定估计尚未证明**：是否存在固定
\(0<\alpha<1,\varepsilon>0,K\ge0\)，使对任意有限闭 Besicovitch cube 族及任意非负 \(g\in L^1\)，只要
\[
 \operatorname{supp}g\subset H,\qquad \int_{Q_i}g\le\lambda|Q_i|,
\]
就有
\[
 \lambda|U\cap\{Mg>\alpha\lambda\}|
 \le(32-\varepsilon)\|g\|_1
 +K\bigl(\lambda\sum_i|Q_i|-32\|g\|_1\bigr)? \tag{5.1}
\]
这里 \(U=\bigcup_iQ_i\)。必须同时允许共享 cubes、任意半径比、任意 core 源分布；第2节正是容量完全饱和时不可遗漏的合法配置。

**若 (5.1) 成立，则严格推出**
\[
 \boxed{C_5\le32-\frac{\varepsilon}{1+31K+32/(1-\alpha)}}. \tag{5.2}
\]
这个传递式是本轮证明，不是 (5.1) 的证明。

### 5.1 有限真实覆盖上的代数

设 \(A\subset E_\lambda\cap U\) 被原好 cubes 覆盖，令
\(D=32m-\lambda|A|=D_{\rm rec}+B\)。由第4节得 \(g\) 与 \(h=f-g\)，且
\(\|h\|_1\le B\le D\)、\(S_g\le31B\)。次线性与已证32基准给
\[
 \begin{aligned}
 \lambda|A|
 &\le\lambda|U\cap\{Mg>\alpha\lambda\}|
   +\frac{32}{1-\alpha}\|h\|_1\\
 &\le(32-\varepsilon)m+
       \left(31K+\frac{32}{1-\alpha}\right)D.
 \end{aligned}
\]
代回 \(D=32m-\lambda|A|\) 并移项，即得 (5.2) 的有限覆盖版本。
这里小源质量只在全空间弱型32中收费一次。

### 5.2 有限严格 Besicovitch 覆盖足以处理全部 \(L^1\) 与尺度

为避免把“有限候选中心”当成一般输入的终点，给出所需覆盖事实。
固定紧集 \(A\subset E_\lambda\)，取 \(0<\lambda_1<\lambda_0<\lambda\)。
联合连续性 \((x,r)\mapsto(2r)^{-5}\int_{Q(x,r)}f\) 在 \(r>0\) 成立；可由 \(L^1\) 逼近连续紧支撑函数证明。
由紧性可选统一 \(\rho>0\)，使每个 \(x\in A\) 至少有一个半径 \(r\ge\rho\) 的平均大于 \(\lambda_0\)。所有平均至少 \(\lambda_0\) 的半径又由总质量统一从上控制。

令 \(r(x)\) 为这些半径的最大值。它存在，且在 \(A\) 上上半连续。
每一步在剩余紧集上选择最大 \(r(x)\)，并删除所选**开** cube 中的接收点。
剩余集仍紧；被选半径不增；不同选中中心的距离至少是较早半径，因而至少 \(\rho\)。在有界紧集内只能选有限多个，故得到覆盖 \(A\) 的有限开 cubes。

最后把全部所选半径稍微缩小同一比例。有限开覆盖与 \(A\) 的紧性保证缩小后的闭 cubes 仍覆盖 \(A\)；平均连续性和 \(\lambda_1<\lambda_0\) 保证它们对阈值 \(\lambda_1\) 仍严格好。原先的中心距离至少为较大旧半径，因此严格大于较大新半径，得到所需闭 Besicovitch 家族。

对该族应用第5.1节，令 \(\lambda_1\uparrow\lambda\)，再由 Lebesgue 内正则性令紧集穷尽 \(E_\lambda\)，即得一般 \(L^1\) 的 (5.2)。没有删去小尺度或大尺度。
这只完成候选估计的传递；真正未闭合的仍是 (5.1)。

## 6. 有限核验、复现与交付状态

每一批在运行前登记，首批至多三个新几何/LP实例，第二批只有一个新几何及其有限选择问题。所有几何、容量、体积和不等式都用整数或 Python `Fraction`；未使用浮点优化、半径网格或接收点采样。

- 首批：新不等半径 clique 的496对中心与80条象限边；65包的2080对包分离和4225个膨胀源入射；精确保留质量 LP 原始/对偶证书。实际 PID **34173**，0.0275秒，退出码0。
- 共享包批：**492528**对中心、**59049**个精确开矩形 cell、全部32个满cliques与源入射、33个完整包选择。增加外角接收缺口有理核验后最终 PID **34935**，0.3491秒，退出码0。全部120秒 checkpoint均可恢复；该时限没有当成总算力上限。

复现：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/astra_cover015/verify.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/astra_cover015/verify_shared.py
```

两个脚本均只需 Python 标准库；`--output-dir` 可写入独立目录。结果与运行收据分开，使结果 SHA 可复现。

| 证书 | SHA-256 |
|---|---|
| `results.json` | `d42af94518c50b5fb3729efe2b7be96b0d298714264e9e80f2f6bd27b5c8f710` |
| `verify.py` | `d9882239ca13eaf98912b8c416ea6ccc1379d39c4829a0545352fa44adf28958` |
| `registration.json` | `fb5e0fcc896cc597d2307f37cab4363cacd2629a7d9055f83f4224b1e10cba82` |
| `shared_results.json` | `d62af1ff3f90ffde44849f9d550e4ea667ccfd77a6a73612a10bf1c637269599` |
| `verify_shared.py` | `a8b16953e0333725545d7587b81ec24e4c8e29e73bdec081fd73c116839d273a` |
| `registration_shared.json` | `e7da25ad17dcceb90806c5d4a2678c9609252b6212ddce172e1098ed9684cb37` |

必要交付为本报告、两个登记、两个验证器、两个结果 JSON 与真实运行/checkpoint收据。独立 clone 的实际 `AGENTS.md` 已读；遵从本次明确分工，只修改 `cube-dim-5-bounds/astra_cover015/`，未提交或推送，未调用共享 MATLAB，未修改主 clone 或三四维材料。

当前真正剩余问题：证明一个能处理共享 critical cliques 的全局接收估计，例如 (5.1)，或另构造能利用真实 \(D_{\rm rec}\) 的选择/分配。局部 \(31.75\)、小平均余量、源32集中及整个包的索引匹配本身都不足以完成此事。
