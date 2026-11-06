# 给 Pro 的研究提示词｜非中心强极大算子端点常数的次指数渐近

研究基线为 2026 年 11 月 5 日。配套《强极大算子端点常数研究笔记》。将下面全文作为研究任务。内部研究材料，公开检索仅用公开术语、定理和书目信息。

## 1. 对象与任务

在 $\mathbb R^n$ 上用 Lebesgue 测度。对非负可测函数定义非中心强极大算子
$$
M_S f(x)=\sup_{x\in R}\frac1{|R|}\int_R f,\qquad R=\prod_{j=1}^n I_j,
$$
其中 $I_j$ 是有界正长度区间，各边长独立。设 $\Phi_n(t)=t(1+\log^+t)^{n-1}$，对数为自然对数，并设 $\Phi_n(0)=0$。所求最优常数
$$
C_n=\sup_{\substack{f\ge0\\0<\int\Phi_n(f)<\infty}}
 \frac{|\{M_S f>1\}|}{\int\Phi_n(f)}.
$$
原阈值 $\alpha>0$ 已通过 $f\mapsto f/\alpha$ 缩放为 1。示性函数保留阈值记法时，分母为 $|E|\alpha^{-1}(1+\log(1/\alpha))^{n-1}$，其中 $0<\alpha<1$。文献中的 $t[1+(\log^+t)^{n-1}]$、Luxemburg 范数、Tauberian 常数及中心矩形的常数各自独立登记。公式、代码和报告统一使用上面的 $\Phi_n$。

**两条主线同时推进。** 第一条是系统遍历、归并测试族，找出能突破双曲示性函数次项的明确机制，并将最佳候选升级为严格下界。第二条是继续优化适用于任意输入的上界，尽量逼近实际下界。用 $\log C_n-n\log B$ 衡量进展。第一层目标是改变 $\sqrt{n\log n}$ 的系数或量级。更高目标是确定 $\log C_n$ 的第二项、判断 $C_n\asymp B^n$，最终争取渐近等价。保持这些目标的逻辑区别。有限维实验用于校准和发现结构。

## 2. 已证基线及复核入口

记 $\mu=\log B$，其中 $B$ 是 $\log B=1+B^{-1}$ 的唯一正解，$B=3.5911214766\ldots$，$\sigma=\sqrt\mu=1.1306920636\ldots$。内部推导得到
$$
\lim_{n\to\infty}C_n^{1/n}=B,\qquad
\limsup_{n\to\infty}C_n/B^n\le\sqrt B.
$$
令 $E_L=\{x\in[-1,1]^n:\sum_j-\log|x_j|\ge L\}$，$H_n$ 是仅用 $\mathbf1_{E_L}$、任意 $L\ge0$ 和 $\alpha\in(0,1)$ 时的最优归一化比值。内部推导得到
$$
\log H_n=n\mu-\sigma\sqrt{n\log n}
+\sigma\{\log(\sigma\sqrt{2\pi})-1\}\sqrt{n/\log n}
+o(\sqrt{n/\log n}),\qquad H_n\le C_n.
$$
第三项系数约为 $0.0472272$。文献首创性单独核查。复核基线后，把主要算力用于突破误差。

### 下界的可复现公式

在 $x_j\ne0$ 处取从 $-r_jx_j$ 到 $x_j$ 的区间，$0\le r_j\le1$。均匀取点给出独立的增量 $Z_{r_j}=-\log|U_j|$，$U_j\sim\mathrm{Unif}[-r_j,1]$。可写成 $Z_r=Y+(-\log r)J$，其中 $Y\sim\mathrm{Exp}(1)$，$J\sim\mathrm{Bernoulli}(r/(1+r))$，二者独立，$r=0$ 时按极限取 $Z_0=Y$。因此
$$
m(r)=1-\frac{r\log r}{1+r},\qquad
v(r)=1+\frac{r(\log r)^2}{(1+r)^2}.
$$
$m$ 在 $r=B^{-1}$ 唯一最大，此时 $m=v=\mu$。用 $\Pr(\sum Z_{r_j}\ge q)=\alpha$ 定义 $q$。在方体内可按随机序约化到这些区间。$L\to\infty$ 时比值趋于
$$
\frac{\alpha e^q}{(1+\log(1/\alpha))^{n-1}},
\qquad
|E_L|=2^ne^{-L}\sum_{k=0}^{n-1}\frac{L^k}{k!}.
$$
固定 $n,\alpha$ 后先令 $L\to\infty$。有限 $L$ 的上界可用包含超水平集的盒 $[-R,R]^n$，$R=2/\alpha-1$。其体积比至多 $\max(e^{q_{\max}},R^n)$。$\alpha\ge1/2$ 时 $R\le3<B$，$\alpha<1/2$ 用锐 $L^p$ 与大振幅截断统一排除。令 $D=n\mu-\sum m(r_i)$、$V=\sum v(r_i)$，则 $\sqrt V\ge\sigma\sqrt n-O(\sqrt D)$。因此把均值和方差合并优化时，收益至多 $O(z^2)$，其中 $z=(\sum m(r_i)-q)/\sqrt V$。这一句限于均值与方差两项，概率尾的余项要另算。统一中偏差下，$1-\alpha=(1+o(1))\overline\Phi_{\mathrm N}(z)$，$\overline\Phi_{\mathrm N}$ 为标准正态上尾。优化 $\sigma\sqrt n\,z+n\overline\Phi_{\mathrm N}(z)$ 给出上述展开，最优量级 $1-\alpha\sim \sigma/\sqrt{n\log n}$。数值计算优先用 Gamma 与二项卷积，别用朴素蒙特卡洛估小尾。 等边长参数 $r$ 的精确实现为
$$
\Pr\!\left(\sum Z_r\ge q\right)=
\sum_{k=0}^n {n\choose k}\pi_r^k(1-\pi_r)^{n-k}
\overline F_{\Gamma(n,1)}(q+k\log r),\qquad
\pi_r=\frac{r}{1+r}.
$$
当 Gamma 尾的自变量非正时取 1。用对数域求和、单调二分求分位数，并记录位数与截尾误差。分块参数先用独立块卷积，非等权离散化须对分位数作双向误差夹逼。

### 一般上界与新定位引理

Grafakos–Montgomery-Smith 的锐强型公式为 $\|M_S\|_{p\to p}=c_p^n$，$(p-1)c_p^p-pc_p^{p-1}=1$。设 $b_p=p\log c_p$，则
$$
p(1-e^{-b_p/p})=1+e^{-b_p},\qquad b_p=\mu+\frac{\mu}{2p}+O(p^{-2}).
$$
逐坐标截断估计可独立使用。对 $g\ge0$、$s>0$、$0<r<1$，
$$
|\{M_Sg>s\}|\le
\frac{(2/(1-r))^n}{s(n-1)!}
\int_{g>sr^n}g\,\bigl[\log(g/(sr^n))\bigr]^{n-1}.
$$
证明入口是 $M_S\le M_1\circ\cdots\circ M_n$，以弱 $(1,1)$ 常数 $2$ 控制 $M_i(g-rs)_+$，再交换分布函数积分。固定 $a>0$，对 $f\le e^a$ 用 $p_n=1+(n-1)\log(1+a)/a$，对 $e^a<f\le e^n$ 用大而固定的 $p$，其余用上式。分配阈值 $1-2n^{-2},n^{-2},n^{-2}$，最后令 $a\downarrow0$，得到 $\sqrt B$ 前因子上界。

可先使用这个一般输入的定量定位。固定 $A>0$，记 $\rho_n=\sqrt{\log n/n}$，$g=f\mathbf1_{\{f>e^{A\rho_n}\}}$，则
$$
|\{M_S g>2n^{-1/2}\}|
\le B^n\exp\!\left[-(A-\sigma)\sqrt{n\log n}+O_A(\log n)\right]\int\Phi_n(g).
$$
验证方法是把 $g$ 在 $e^{\sqrt n}$ 再切一次。中段取 $p=\sigma\sqrt{n/\log n}$，阈值 $n^{-1/2}$，对 $u=\log t$ 用 $(p-1)u-(n-1)\log(1+u)$ 的凸性取两端。大段以 $r=1/2$ 用截断核，其系数的 $n$ 次方根为 $O(n^{-1/2})$。所以对目标 $\log C_n\le n\mu-c\sqrt{n\log n}+o(\sqrt{n\log n})$，取 $A>\sigma+c$ 后，核心转为有界振幅 $f\le e^{A\rho_n}$、阈值 $1-2n^{-1/2}$ 的统一估计。此定位没有给出核心估计。幅值接近阈值与矩形平均接近饱和需要分别处理。

## 3. 测试族遍历

每一族先写出输入、参数、实际测度、矩形见证、阈值和分母，再给出无量纲评分。先复现对照组，然后按下列顺序进入主攻。

| 层级 | 测试族 | 参数与关键问题 |
|---|---|---|
| 对照 | 单盒、乘积输入、双曲单层 | 对积型输入复核 $M_S(\bigotimes f_i)=\bigotimes M_1f_i$。用单盒校准体积和小阈值。以 Gamma 公式复现 $H_n$。 |
| A1 | 双曲多层与软振幅 | $f=h(\sum-\log|x_i|-L)\mathbf1_{[-1,1]^n}$，从单调三至五段阶梯开始。优化层宽和高度，特别比较 $t-1$ 在 $1/n$、$1/\sqrt{n\log n}$、$\rho_n$ 三个尺度。同步推导连续极限。 |
| A2 | 加权和分块 | $\sum w_i Z_i$、块大小与块间聚合，测试 max、k-of-m、分块阈值及交并。要求在保持平均增益 $n\mu-o(n)$ 的前提下，验证方差、下尾或对超阈值幅值的代价是否改变。 |
| A3 | 相关层与多中心 | 多个平移或各向异性伸缩的双曲层做并、交、叠权，测试具有显式重叠控制的格点和块图。对每个模式同时记录输入重叠与接收集合重叠。 |
| A4 | 周期、格点及有限单元输入 | 先用可积的有限盒截断，再比较周期块的边界与内部。可扫二值、三值、跨坐标相关编码和盒内薄层。高幅尖峰以实际密度和支持体积入模。 |
| B | 自由网格与反例搜索 | 从 2–4 维搜索结构并跨维放大，和匹配的独立、打乱、等权、等测度对照比较。小维改进本身不代表高维次项改进。 |

**A1 的解析入口。** 令 $u=\sum-\log|x_i|-L$，固定维数再令 $L\to\infty$。对于尾部可积且极限交换已证明的轮廓，输入模数的归一化极限为 $\int_{\mathbb R} e^{-u}\Phi_n(h(u))\,du$，矩形平均对应 $\mathbb E\,h(u+\sum Z_{r_i})$。分母积分若不收敛，该轮廓须另作截断。先以固定 $\mathbf r$ 解有限阶梯的凸问题，再检查自由矩形和非单调轮廓的适用边界。

统一记录
$$
L(f,\alpha)=\frac{|\{M_S f>\alpha\}|}{\int\Phi_n(f/\alpha)},\quad
G_n=n\mu-\log L(f,\alpha),\quad
G_n/\sqrt{n\log n}.
$$
若分子采用见证集合，则以认证下界替代真实体积并明确标注。每条拟合都与同维数的 $H_n$ 数值基线比较，记录 $n$、阈值、尺寸、随机种子、网格分辨率和精度。分别检查固定维数的截断极限与跨维外推。优先筛选能证明 $G_n\le(\sigma-\varepsilon)\sqrt{n\log n}$ 或 $G_n=o(\sqrt{n\log n})$ 的机制。

## 4. 可认证搜索

用有限网格单元 $Q$ 的非负常值 $a_Q$ 表示输入。对给定的矩形列表 $\mathcal R$ 与严格余量 $\eta>0$，求解
$$
\min_{a_Q\ge0}\ \sum_Q |Q|\Phi_n(a_Q),\qquad
\sum_Q |Q\cap R|a_Q\ge(1+\eta)|R| \quad(R\in\mathcal R).
$$
这是凸目标与线性约束。认证接收测度为 $|\bigcup_{R\in\mathcal R}R|$。离线用有理数复核交体积和所有平均，取分母的定向上界，用并集的定向下界计算证书。网格未覆盖的矩形只会遗漏可用见证。有限盒、步函数、矩形严格余量和紧集穷尽可将证书上确界提升到 $C_n$，这个归约没有维数一致的收敛率。可采用“优化振幅、重选见证”的交替迭代与多起点。数值局部最优和证书下界分开登记。以原子或 Dirac 只作几何代理，任何发表级候选均需回到可积密度并结算实际 Orlicz 模数。

## 5. 上界主攻

1. **先做阈值附近的核心估计。** 使用定位引理剔除大振幅。分别证明示性函数、高度略超阈值的近示性函数、分层一般输入的结论，报告每次提升的损失。示性函数的 Tauberian 估计若无法定量转为一般输入，就保留为独立定理。
2. **研究锐 $L^p$ 证明的稳定性。** 回到一维 rising-sun 分布不等式与 Hölder 等号，追踪切片上的缺口。尝试在 $p$ 随 $n$ 增长的范围内，把近饱和导向概率反集中、积分重叠代价或结构近似。核对坐标最大算子复合允许逐步适配，实际矩形必须同时确定所有边。任何缺口都要说明针对哪一个算子。
3. **利用直接矩形几何。** 尝试选择覆盖、联合矩阵或矩形线性化，将输入和接收测度同时保留。先在有限矩形字典寻找可审计的对偶证书，再分析常数是否维数一致。单纯重选 $p$ 或把坐标损失相乘，除非导出新的统一误差界，均计为基线优化。
4. **用下界检验上界假说。** 若猜测普遍有负的 $\sqrt{n\log n}$ 次项，必须与多层和相关测试对照。若发现 $C_n\gtrsim B^n$，转向收紧前因子、判断归一化极限及极值结构。允许上下界路线各自得出部分结论。

## 6. 文献核查与输出

- Jessen–Marcinkiewicz–Zygmund，*Note on the differentiability of multiple integrals*，1935，[DOI](https://doi.org/10.4064/fm-25-1-217-234)。
- Fava，*Weak type inequalities for product operators*，1972，[DOI](https://doi.org/10.4064/sm-42-3-271-288)。de Guzmán，*An inequality for the Hardy-Littlewood maximal operator with respect to a product of differentiation bases*，1974，[DOI](https://doi.org/10.4064/sm-49-2-185-194)。Córdoba–Fefferman，*A geometric proof of the strong maximal theorem*，1975，[期刊页](https://annals.math.princeton.edu/1975/102-1/p07)。
- Grafakos–Montgomery-Smith，*Best Constants for Uncentred Maximal Functions*，1997，[论文](https://stephenmontgomerysmith.github.io/preprints/uncentered.pdf)。[DOI](https://doi.org/10.1112/S0024609396002081)。
- Qin–Wu–Xue，*The limiting weak-type behaviors of the strong maximal operators*，2021，[期刊页](https://doi.org/10.1007/s43037-021-00154-6)。[预印本](https://arxiv.org/abs/2104.03497)。区分固定 $n$ 时 $\alpha\downarrow0$ 的 $2^n/(n-1)!$ 与本任务中 $n\to\infty$、阈值趋近 1 的机制。
- Hagelstein–Parissis，[Solyanik estimates in harmonic analysis](https://arxiv.org/abs/1310.3771)。Hagelstein，[Current perspectives on the Halo Conjecture](https://doi.org/10.1515/ans-2023-0157)。查证全文、定理号及归一化，记录与本题相同的量和可迁移的量。公开优先权仅据已核实的文献表述。

每轮提交六项内容。① 日期与结论状态。② 严格改进或障碍的公式。③ 测试族表，列出证书与数值各自的最优值及参数。④ 可复核推导、可复现实验设定或反例。⑤ 目前上下界及其差距。⑥ 下一轮最有区分力的两项动作。最终交付一般输入上界、最优下界或已证障碍，附完整证明和出处。把仅在测试族内成立的障碍、数值现象和一般定理分栏。每轮以突破基线或排除明确机制为收尾标准。
