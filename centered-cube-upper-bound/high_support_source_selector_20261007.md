# 高 support 原来源选择器：Gamma 大偏差支支付与 regular 空间核心

2026-10-07。此稿研究固定 \(c,\ell\) 的原有序核，不改变主账。新的严格结果是一个保留原来源与任意 receiver 选择器的正核支付支；没有证明剩余 regular 时钟区的弱端点。三轮新证书只核该支付支的原 Exp/Gamma 潜变量概率与 Bernoulli mask 质量，没有重新执行旧数值。

## 1. 原对象、已付支与范围

固定
\[
 B(v)=v/\log(1+v)-1,\qquad
 G_i=(I+cB(-\ell^2\partial_{ii}))^{-1},\qquad
 S=\sum_{i=1}^n(I-G_i).
\]
令
\[
 K_r=\bigotimes_{i=1}^n[(1-r)I+rG_i],
 \quad H_r=e^{-rS},\quad D_r=K_r-H_r,\quad 0\le r\le1.
\tag{1}
\]
来源始终是障碍分解的同一 \(\nu_b\ge0\)，质量 \(W_b\)，支撑于
\(\Omega=\{u>0\}\)。原障碍合同为
\[
 Su=\nu_b-\mu_b,\qquad 0\le\mu_b\le\kappa,\qquad
 \int\mu_b=\int\nu_b=W_b.
\tag{2}
\]
源分布可任意相关，不假设其坐标独立。由正 Markov 核，
\[
 |D_r\mu_b|\le\kappa,\qquad
 \{x\in\Omega^c:\sup_r D_rSu>2\kappa\}
 \subset\{x\in\Omega^c:\sup_r D_r\nu_b>\kappa\}.
\tag{3}
\]
不使用已被 spike 排除的势能上界。

按实际非零跳跃坐标 support 分解 \(D_r=\sum_A D_{r,A}\)。
原 \(G_i\) 连续核使该空间面分解合法；以下潜变量证明也可以直接按 Poisson 标签分解成立。这里的标签是算子内部跳跃坐标，不是来源原子的坐标结构或 FIRST 标签。

取 \(K=\lfloor n^{1/3}\rfloor\)、\(b=K/(4n)\)。已有
late_multiface_source_budget_20261007.md 的两支为：全部 \(r\) 的
\(|A|\le K\) 正包络费用至多 \(W_b\)，以及 \(r\le b\)、
\(|A|>K\) 的正包络费用至多
\(\tfrac12(2/3)^K W_b\)。此稿不重复验证或将这些费用重复登记。
真正剩余对象为
\[
 B_{\rm high}(x)=\sup_{b\le r\le1}
       [D_r^{[>K]}\nu_b(x)]_+,\qquad x\in\Omega^c.
\tag{4}
\]

## 2. 同一个原 Gamma 时钟下保留差核

写
\[
 \Pi_T^{(i)}=\exp[-cT B(-\ell^2\partial_{ii})].
\]
原 \(B\) 是 Bernstein 函数，故 \(\Pi_T^{(i)}\) 是正、对称、守恒 Markov 卷积核，且
\[
 G_i=\int_0^\infty e^{-T}\Pi_T^{(i)}\,dT,\qquad
 G_i^j=\int_0^\infty
       \frac{T^{j-1}e^{-T}}{(j-1)!}\Pi_T^{(i)}\,dT.
\tag{5}
\]
这是真实原 resolvent 的 Exp/Gamma 从属表示，而非 reset 替代。

对 \(A\subset[n]\)、\(k=|A|\)，令
\[
 \Pi_{\mathbf T}^A=\prod_{i\in A}\Pi_{T_i}^{(i)},\qquad
 p_{n,k}(r)=r^k(1-r)^{n-k},\qquad
 F(z)=\sum_{j=0}^\infty\frac{z^j}{(j+1)!j!}.
\]
在共同参考测度 \(e^{-\sum_{i\in A}T_i}d\mathbf T\) 下，
\[
 D_{r,A}
 =\int_{\mathbb R_+^k}e^{-\sum T_i}
 \left\{p_{n,k}(r)-e^{-nr}r^k\prod_{i\in A}F(rT_i)\right\}
 \Pi_{\mathbf T}^A\,d\mathbf T.
\tag{6}
\]
这是先完整合并原 Poisson 重复跳跃后得到的公式，保留跨次数消去。
其中 \(H\) 项来自每个 active 坐标至少一跳：
\[
 e^{-r}\sum_{j\ge1}\frac{r^j}{j!}G_i^j
 =\int_0^\infty e^{-r-T}rF(rT)\Pi_T^{(i)}\,dT.
\]
有限 support 总和与所有积分对正来源可分别 Tonelli；带符号差通过两项各自的有限质量定义。

以下稀有时钟概率在 \(K\) 的独立 Exp 参考测度下计算。此内部表示不要求输入来源独立；\(H\) 的时钟分布不是 Exp 乘积，不能将其概率直接写成相同 Gamma 尾概率。

## 3. 可付的真实 source-weighted 稀有时钟支

定义同一个共同标记
\[
 \mathcal B_A=\{\sum_{i\in A}T_i\le k/4\}
               \cup\{\sum_{i\in A}T_i\ge4k\}.
\tag{7}
\]
分别在式 (6) 的 \(K,H\) 两项中使用此标记，得到
\(D_{r,A}^{\rm bad}\) 和其补集上的 \(D_{r,A}^{\rm reg}\)。
没有用另一个时钟区替换负项。对非负来源，
\[
 D_{r,A}^{\rm bad}\nu_b\le
 K_{r,A}^{\rm bad}\nu_b
 =p_{n,k}(r)\,G_A^{\rm bad}\nu_b,
\quad
 G_A^{\rm bad}
 =\int_{\mathcal B_A}e^{-\sum T_i}
          \Pi_{\mathbf T}^A\,d\mathbf T.
\tag{8}
\]
\(G_A^{\rm bad}\) 是固定的联合时钟 subprobability 卷积核。
其时钟标记使 active 坐标相关；不将它归一化成新的来源再收 \(W_b\)。

令
\[
 q_k=\Pr\{\Gamma(k,1)\le k/4\}
           +\Pr\{\Gamma(k,1)\ge4k\}.
\]
Chernoff 给出
\[
 \Pr\{\Gamma(k,1)\le k/4\}
 \le(e^{1/4}/2)^k\le(2/3)^k,
\]
\[
 \Pr\{\Gamma(k,1)\ge4k\}
 \le(2/e^2)^k\le(1/3)^k.
\tag{9}
\]
前者取 \(e^{-T}\) 的指数矩，后者取 \(e^{T/2}\)；
\(e^{1/4}\le4/3\)、\(e^2>6\) 已足够，不寻求最优常数。

定义一个与 receiver 参数无关的固定正包络
\[
 P_{\Gamma,K}
 =\sum_{|A|>K} \beta_{n,|A|}\,G_A^{\rm bad},\qquad
 \beta_{n,k}=(k/n)^k(1-k/n)^{n-k}
            =\max_{0\le r\le1}p_{n,k}(r).
\tag{10}
\]
在 \(k=n\) 处按 \(0^0=1\) 解释。原 Bernoulli 概率质量
\(\binom nk\beta_{n,k}\le1\)，所以
\[
 \|P_{\Gamma,K}\nu_b\|_1
 =W_b\sum_{k=K+1}^n\binom nk\beta_{n,k}q_k
 \le C_\Gamma(K)W_b,
\]
\[
 C_\Gamma(K)
 =3(2/3)^{K+1}+\tfrac32(1/3)^{K+1}
 \le17/18\quad(K\ge2).
\tag{11}
\]
这里逐 mask 取最大值之所以没有 \(\sqrt n\) 或 \(2^n\) 损失，
是因为该 mask 的真正时钟 subkernel 质量 \(q_k\) 随 \(k\) 指数下降。
没有将此论证套到 regular 核。

因而对任意可测选择器 \(r_*(x)\in[0,1]\)、任意
\(0\le\theta(x)\le1\) 及任意可测 receiver 集合 \(E\)，有
\[
 \int_E\theta(x)
   K_{r_*(x)}^{[>K],{\rm bad}}\nu_b(x)\,dx
 \le\int P_{\Gamma,K}\nu_b
 \le C_\Gamma(K)W_b.
\tag{12}
\]
这是真实同源、同 receiver 参数的 source-weighted 预算。
不要求逐源点整个 dual column 有界，不用条件弱范数平均，
也未按输出选择的 mask 或时钟重新领取源质量。

## 4. 回代剩余目标时保留符号

对所有 \(r\)，
\[
 D_r^{[>K]}\nu_b
 \le P_{\Gamma,K}\nu_b+D_r^{[>K],{\rm reg}}\nu_b.
\tag{13}
\]
因此
\[
 \{x\in\Omega^c:B_{\rm high}(x)>\kappa\}
 \subset
 \{P_{\Gamma,K}\nu_b>\kappa/2\}
 \cup
 \{x\in\Omega^c:
       \sup_{b\le r\le1}D_r^{[>K],{\rm reg}}\nu_b>\kappa/2\}.
\tag{14}
\]
第一集合的 \(\kappa\)-体积至多 \(2C_\Gamma(K)W_b\)，
第二集合仍未付。此处损失 2 来自阈值分配，不能把它遗漏。

也可用近赢家 \(r_*(x)\) 线性化式 (4)，再取
\(\theta=\kappa/B_{\rm high}\) 于其超水平集。
式 (12) 直接支付 bad 部分；regular 部分须保留如下有符号来源侧积分：
\[
 \int\nu_b(dy)\sum_{|A|>K}
 \int_{\{k/4<\sum T_i<4k\}}e^{-\sum T_i}
 \int_E\theta(x)
 \left[p_{n,k}(r_*(x))
       -e^{-nr_*(x)}r_*(x)^k\prod_iF(r_*(x)T_i)\right]
 \Pi_{\mathbf T}^A(dx-y)\,d\mathbf T.
\tag{15}
\]
边界等号在时钟绝对连续测度下无质量。一般 \(L^1\) 来源可用有理参数
sup 和近赢家，随后单调极限；不需要未经证明的实际历史 winner 定义。
式 (15) 未按 \(y\)、\(A\)、\(\mathbf T\) 或 receiver 分拆收费。

## 5. regular 原空间核心的结论与限制

式 (15) 中 \(r_*(x)\) 留在空间积分内部，\(\Pi_{\mathbf T}^A\nu_b\)
由同一来源产生。它不是任意独立 Bernstein 通道，也不是可以直接
对其应用 posterior 驻点的正多项式；原 \(H\) 的整个负项必须保留。
总时钟 regular 条件只限制 \(\sum T_i\)，并不控制每个坐标的时钟。
其后原空间卷积和 \(\Omega^c\)、\(\nu_b\) 支撑条件才可能提供进一步费用。

我未证明下列任一新强合同：regular 核的整体弱型、
其来源加权 signed 选择器积分的 polylog 预算、
或由 regular 时钟约束得到冻结 \(H\) 的共同正支配。
也未将潜变量密度比的坏构型称为原空间弱型反例；
Markov pushforward 可以显著降低潜变量区分。

已查重的限制继续有效：任意 Bernstein 通道的 posterior-only
弱端点、总 test-energy 和全空间时间配对合同已被各自守卫否定；
原 \(N=(K-H)(I+S)\) 的 face-TV 下界也不能用于否定本式的域外弱端点。
Rota 投影后的条件最大值不能无损转移弱范数。
本轮没有用这些已知障碍冒充新进展。

capped_reference_absorption_20261007.md 降低参考合同所需强度：
若未来从式 (15) 产生同一原来源的一个**整体**非负 reference，
其独立弱预算为 \(A W_b\)，并有正确实际行 cap 和可付误差，
则可通过截帽层蛋糕与主账自吸收收费，无需 reference 的强 \(L^1\)。
此引理仍不证明 regular reference 的弱预算，不能把每个时钟或每个
参数的条件弱预算平均后代入。固定 \(c,\ell\) 的结果也没有自动接入
移动尺度、原 FIRST/CP/GP 等 actual geom 门。

因此新完成的是式 (12) 的非零有偿子支和精确式 (15)；
一般 regular 空间核心、\(A_{\rm ord}\) 的 polylog 端点及 actual geom
转移仍未闭合。不能因为 \(C_\Gamma(K)\) 小就删除式 (15)。

## 6. 新三轮证书与收据

预登记 high_support_source_selector_registration_20261007.json 后执行
high_support_source_selector_guard_20261007.py，同一 session 77660 完成，
exit 0。n=8/32/128，K=2/3/5，exp 正 Taylor 项数128/512/1536；
分别56/240/992项，共1288项
PASS_EXACT_RATIONAL_INTERVAL。没有旧数值重跑、浮点拟合或实际来源采样。

证书用
\[
 \Pr\{\Gamma(k,1)\le x\}
 =1-e^{-x}\sum_{j=0}^{k-1}x^j/j!
\]
及其上尾公式，对每个 \(K<k\le n\) 以 Fraction 外包
\(e^{-k/4}\)、\(e^{-4k}\)，核验两尾、每 mask 的实际 Bernoulli
峰质量、有限包络质量和式 (11)。正 Taylor 下界与首遗漏项的几何
上界先包住 \(e^x\)，再取精确有理倒数，CDF 消去也保留外包方向。

| n | K | checks | \(C_\Gamma(K)\) |
|---:|---:|---:|---:|
| 8 | 2 | 56 | \(17/18\) |
| 32 | 3 | 240 | \(11/18\) |
| 128 | 5 | 992 | \(43/162\) |

原 \(\Pi_T^A\) 的 Markov 守恒和同源选择器 Tonelli 支配由解析证明负责，
未在证书中用 reset 矩阵替代空间核。结果没有认证式 (15) 的弱界或
actual FIRST 的非空事件。

注册 SHA256：
4e4cca666cae278b7078492d8d8160c1e7c19992516dd901cdbb97da36e08b99。

脚本 SHA256：
01387b5f3d472b065cc729f035569ae97489de4fbc619217b99bdfa4cf2b206b。

结果文件 high_support_source_selector_results_20261007.json 保存每个
Gamma 尾的有理区间、各层峰质量与三轮包络质量外包；脚本和结果
至此冻结。

