# R050：八阶实际稳定性之后的全高支撑 gross 小量

2026-11-05。条件性新命题，供独立审核。固定 \(n\ge3\)。仓库、论文及其他代理文件只读。

本稿明确以 R049 冻结稿

/Users/zhengzhihao/Desktop/T/topic1-tauberian-work/scratch/r049-second-variation-upper/competitive_eighth_stability.md

SHA256 **5c4d7b371a0c33c39449afd686de599c49570aaa639d4f7db7d3c1f166a040f4**

的稳定性结论为新增条件。它尚未完成独审前，本文只作条件推论，不升级真实八阶精确系数。R043、R044 的几何则复用其已审正式版本。

## 1. 同层对象和条件的准确口径

沿 R049 (1) 的八阶竞争实际输入序列，一次选取同一新 \(N=8\) 族。保持原有限严格见证、实际源 \(F_i\subset G_i\)、质量根 \(R_i=(|G_i|/a)^{1/n}\)、原／有效／中点 CDF 盒和全部外部 pair。记
\[
 a=v^n,\quad S=\sum_i r_i^n,\quad
 \mathscr D_i=(x_i+W_i^0)\setminus F_i,\quad
 D_i'=(x_i+W_i')\setminus F_i.
\]
正常集合 \(\mathcal N\)、正常实际层 \(Y_{i,N}\)、全部有界实际层 \(Y_{i,\rm all}\)、最后大管 \(K=K_{\rm all}\) 及正常 clean 槽分配均为 R045/R049 的**同一组实际对象**。特别
\[
 Y_{i,N}\subset Y_{i,\rm all}\subset D_i'\cap\mathcal C_i(K),\qquad
 R_N=\frac va\sum_i|\mathscr D_i\setminus Y_{i,N}|,
\]
\[
 Z(K)=\frac va\sum_i|D_i'\setminus\mathcal C_i(K)|
       \le R_{\rm all}\le R_N.
\tag{1}
\]
\(\mathcal C_i(K)\) 是整个两近面管，不能换成某个单独角臂。

在本次新八层统一 CDF 及固定工具之前，取
\(\tau_*=1/128<1/100\) 和 R043 的固定窗口 \(R_*=100(n+1)\)。
R049 的新增稳定性条件给
\[
 \boxed{
 R_N+R_{\rm all}+\mathfrak E+\mathcal F+\mathcal F_{\tau_*}
           +v^6D_N\le C_*v^8S,}
\tag{H49}
\]
\[
 \mathfrak E=\sum_iR_i^{n-1}\Delta_i^{\rm end}.
\]
这里 \(R_N\) 的估计由稳定量同时保留 \(B_{\rm bad}\) 及
\(R_N\le R_{\rm all}+B_{\rm bad}\) 得到，不从 \(R_{\rm all}\) 单独补造。常数 \(C_*\) 使用 R049 证明中**先固定的一次**高支撑份额 \(\eta_0>0\)；下文另引任意 \(\eta>0\) 时不改变这个构造或 \(C_*\)。

正常实际层的原几何同时给逐源上界
\[
 |Y_{i,N}|\le C a v^5R_i^n.
\tag{2}
\]
这是正常窗内有限角数、薄轴 \(O(v^4R_i)\)、伴随轴 \(O(v^3R_i)\) 和其他原宽 \(O(vR_i)\) 的实际体积，不要求真实缺集逐源与标准臂接近。

令 \(q(e)\) 为原实际相交 pair 的临界轴数，\(I_e=|X_i\cap X_j|\)。深处 pair 的所有者总取 \(R_i\ge R_j\) 一端，同根按预定编号。R043 的真实缺槽坏根集记为 \(\mathcal B\)；由其定义及 (H49)，
\[
 S(\mathcal B)\le C[D_N+v^{9/2}S]=O(v^2S).
\tag{3}
\]
不再把任何没有定量权控制的 LT 异常根加入这个集合。

**条件命题。** 在上述同层实际输入和 (H49) 下，
\[
 \boxed{\displaystyle \sum_{e:\,q(e)\ge3}I_e=o(v^8S).}
\tag{4}
\]
和含全部深处及外部边界高支撑 pair，每个无向 pair 只计一次。结论是 gross 交量，不是再次减去完整 \(\sum J_i\) 的净价。二维没有这些 pair，结论以空和为零成立，毋须套用高维引理。

## 2. 缺集分解给出的局部 epsilon 估计

先明确实际源质量比较：
\[
 \mu_i=|F_i|\ge(1-Av)|G_i|,\quad
 \varepsilon_i\mu_i=|W_i|-\mu_i\le|\mathscr D_i|.
\]
故对小 \(v\)，
\[
 \varepsilon_iR_i^n
 \le C|\mathscr D_i|/a
 \le Cv^5R_i^n+\frac C a|\mathscr D_i\setminus Y_{i,N}|.
\tag{5}
\]
不能把粗先验 \(|\mathscr D_i|\le C\varepsilon_i aR_i^n\) 反过来使用；(5) 来自实际 \(\mu_i,g_i\) 和原盒包含。

对任意指标子集 \(A\)，令 \(\operatorname{Adj}(A)\) 为与某个 \(i\in A\) 发生实际 \(X\) 相交的端点集合，含 \(A\) 本身。原相交度和相交根比给
\[
 \sum_{j\in\operatorname{Adj}(A)}R_j^n\le C S(A).
\]
由 (5)、正常层的 (2) 及 (H49)，有统一估计
\[
 \sum_{j\in\operatorname{Adj}(A)}\varepsilon_jR_j^n
 \le Cv^5 S(A)+\frac C v R_N.
\tag{6}
\]
只在第二项用**全局实际 residual** 的上界，不要求其逐根分布均匀。

对 \(A=\mathcal B\)，(3) 给
\[
 \sum_{\operatorname{Adj}(\mathcal B)}
          \varepsilon_jR_j^n=O(v^7S).
\tag{7}
\]
对外部边界集合 \(\partial I\)，原 \(S(\partial I)=O(v^5S)\) 给
\[
 \sum_{\operatorname{Adj}(\partial I)}
          \varepsilon_jR_j^n
 \le O(v^{10}S)+O(v^7S)=O(v^7S).
\tag{8}
\]
所有这些常数与下文尾参数 \(\eta\) 无关。

## 3. 好所有者：实际稳健面路径给九阶以上

对 \(i\notin\mathcal B\)，R043 使用同一固定完整窗口中的长度量
\[
 \mathcal E_i=
 \sum_{\substack{c,d\in\mathcal I_i(R_*)\\
                  \text{相关固定稳健面对}}}\Gamma_{cd}
       +\sum_{c\in\mathcal I_i(R_*)}\Delta_c^{\rm end}.
\]
该几何清单原有切向裕量可取 \(1/16\)，所以这些 pair 也属于本稿先固定的 \(\tau_*=1/128\) 候选。原半盒法向缝、半径比及深处条件不变。固定窗口的数目与重复次数有界，故 (H49) 把 R043 原来的六阶均值加强为
\[
 \sum_{i\notin\mathcal B}\mathcal E_i r_i^{n-1}
 \le C(\mathcal F_{\tau_*}+\mathfrak E)=O(v^8S).
\tag{9}
\]
没有将一个原始面邻居的缝直接当成窗口内所有稳健缝。

R043 的完整面对路径、实际 helper 身份和原中心严格分离，仍逐边给
\[
 I_e\le C_{n,A,T_0}
          [v^{3q}r_i^n+v^{q-1}\mathcal E_i r_i^{n-1}],
 \qquad q=q(e)\ge3,\quad\operatorname{owner}(e)=i\notin\mathcal B.
\tag{10}
\]
该式使用所有临界轴的实际端点传递，不是把一条薄横截面占用率乘到高维盒上。对有限 \(3\le q\le n\)，有限度与 (9) 给
\[
 \sum_{\substack{q(e)\ge3,\ \mathrm{deep}\\
                  \operatorname{owner}(e)\notin\mathcal B}}I_e
 \le C[v^9S+v^2O(v^8S)]=O(v^9S).
\tag{11}
\]
这里 \(q=3\) 的首项正是九阶，\(q\ge4\) 更小。

## 4. 坏所有者：先分掉实际质量和面对例外

取两端都深、所有者 \(i\in\mathcal B\) 的 \(q\)-临界边，沿 R044 写
\[
 0<p_1\le\cdots\le p_q\le C_pv,\qquad
 p_l=L_l/r_i,\quad A=\prod_{l=1}^q p_l,\quad
 P=\prod_{l=2}^q p_l,\quad H=\min(p_1,A/v^q).
\]
有 \(I_e\le C A r_i^n\)。这些是实际最大根所有者上的实际交长，未替换为平均半径。

原 R044 目标缺量例外
\(\varepsilon_j\ge c_0A/v^q\) 给
\[
 I_e\le Cv^q\varepsilon_jr_i^n.
\]
对坏所有者的目标邻居用 (7)，有限反向度及尺度可比给
\[
 \sum_{\text{这些目标例外}}I_e=O(v^{q+7}S)=O(v^{10}S)
 \quad(q\ge3).
\tag{12}
\]
这一步特意没有沿用原全族 \(v^q\sum\varepsilon=O(v^{q+5})\)：正常层的 \(av^5\) 只乘坏所有者邻目标的 \(O(v^2S)\) 权，剩下缺量由 \(R_N\) 付体积控制。

其余边若 \(\alpha_i+\Delta_i^{\rm end}/r_i\ge c_1H\)，原关系
\(A/H=\max(P,v^q)\le Cv^{q-1}\) 给
\[
 \sum_{\text{面对／端移例外}} I_e
 \le Cv^{q-1}(\mathcal F+\mathfrak E)
 =O(v^{q+7}S)=O(v^{10}S).
\tag{13}
\]
这里仍只使用实际面对点态半径容量；没有引入新的 \(J_i\)。

## 5. 坏所有者的实际矩形：截宽与共同重数

对上述两类例外之外的边，R044 在原有效盒内给实际缺源矩形：
\[
 d_1\in[c_3H,2c_3H]r_i,\qquad
 d_l\in[c_2p_l,2c_2p_l]r_i\quad(2\le l\le q),
\]
其余轴取到两面均 \(>Kv^3r_i\) 的核心。参数 \(c_2,c_3,c_0,c_1\) 先按原固定几何次序选定。

为明确实际尺度接口：原始面对端距
\((R_i+R_j)/2-L_l\le R_i-L_l\)，原指定见证半径大于 \(R_i\)，端移受非例外条件控制，故其各临界交长至少是
\(\min(p_lr_i/2,w^0_{jl})\ge\gamma p_lr_i\)，\(\gamma>0\) 仅依原宽上下界；非临界轴捕获目标完整原宽。因此同一个见证看见目标实际 \(F_j\) 的固定比例 \(a r_i^n A/v^q\)。它与完整 \(G_i\) 排他，而真正浅面的实际容量上界只容许 \(Ca r_i^nH\) 的足够小固定倍数，矛盾。这证明整块确在 \(D'_i\)，不只是一条角点不等式。

固定 \(K\) 后取 R044 的 \(T_K>2K/c_2+1\)。若 \(p_2>T_Kv^3\)，矩形至多第一轴近某面，其余轴同时离两面远，故避开整个 \(\mathcal C_i(K)\)。其真实成本
\[
 X_e=\frac va|\mathcal R_e|
 \ge c r_i^n\frac{HP}{v^{q-1}}.
\tag{14}
\]
同一个所有者的全部相交边（跨所有 \(q\)）至多固定 \(D\) 条，故
\(\sum_{\mathrm{owner}=i}X_e\le D Z_i(K)\)。
横向核心对小 \(v\) 保留固定比例；其比例常数不依固定 \(K\)，只有起效阈值依赖 \(K\)。\(q=n\) 的横积为空积仍有效。

因此，对任意**另行给定的固定** \(\eta>0\)，R044 的逐根两分支净价给
\[
 \sum_{\substack{\mathrm{nonexception},\,q\ge3,\ p_2>T_Kv^3\\
                   \mathrm{owner}\in\mathcal B}}I_e
 \le \eta Z(K)+
       C_\eta\sum_{q=3}^n v^{\beta_q}S(\mathcal B),
 \qquad \beta_q=\frac{q(2q-1)}{q-1}.
\tag{15}
\]
这是从同一个实际矩形体积得出的估计。总度数已统一分份，不能每个 \(q\) 分别用一整份 \(\eta Z\)。

因 \(\beta_3=15/2\)、\(\beta_q\ge28/3\) 对 \(q\ge4\)，(3) 给长窗总量
\[
 \le \eta Z(K)+O_\eta(v^{19/2}S).
\tag{16}
\]
对短窗 \(p_2\le T_Kv^3\)，不用任何矩形预算，排序给
\[
 I_e\le C_Kv^{q+4}r_i^n.
\]
现在只对坏所有者求和，用 (3) 而不是旧全根估计，得
\[
 \sum_{\substack{\mathrm{nonexception},\,q\ge3,\ p_2\le T_Kv^3\\
                   \mathrm{owner}\in\mathcal B}}I_e
 =O_K(v^9S).
\tag{17}
\]
特别 \(q=4\) 的短窗是 \(v^8S(\mathcal B)=O(v^{10}S)\)，不残留一个全根 \(O(v^8S)\)。

## 6. 外部边界的实际 epsilon 项

两端不都深的 \(q\ge3\) pair，使用 R014/R044 原单源界
\[
 I_e\le Cv^q
       [\varepsilon_i+\varepsilon_j+v^{q/(q-1)}]r_i^n.
\tag{18}
\]
对 epsilon 部分必须使用 (8)。由有限度和尺度可比，
\[
 \sum_{\text{边界 }q\text{-边}}I_e
 \le O(v^{q+7}S)
       +C v^{q+q/(q-1)}S(\partial I).
\]
第二项为
\(O(v^{q+q/(q-1)+5}S)\)，其最小指数在 \(q=3\) 为 \(19/2>8\)。于是全部高支撑边界交量为
\[
 O(v^{19/2}S)+O(v^{10}S)=o(v^8S).
\tag{19}
\]
若只用原全族 \(\sum\varepsilon=O(v^5S)\) 会在三临界处停在八阶；本稿没有这样丢掉目标分解。边界端若没有正常层，相关缺集自动全部属于 \(R_N\)，公式 (5)–(8) 仍成立。

## 7. 先固定稳定构造，再令尾份额趋零

合并 (11)–(19)，对每个固定 \(\eta>0\)：
\[
 \sum_{q(e)\ge3}I_e
 \le \eta Z(K)+O_K(v^9S)+O_\eta(v^{19/2}S)+O(v^{10}S).
\tag{20}
\]
对同一个 R049 实际构造，(H49) 给 \(Z(K)\le C_*v^8S\)，其中 \(C_*\) 只使用早先固定的 \(\eta_0\)，与这里可任意选的 \(\eta\) 无关。故
\[
 \limsup_{v\to0}\frac{\sum_{q(e)\ge3}I_e}{v^8S}
 \le C_*\eta.
\]
先取 \(v\to0\)，随后 \(\eta\downarrow0\)，即得 (4)。不需任何 \(\eta(v)\)、增长深度、增长管宽或量化 LT 速率。

式 (15) 并非从另一份完整 \(J_i\) 购买交量：\(Z\) 是同一实际缺集的既定几何量。这里先用 R049 已完成的正账本控制该几何量，再把 R044 的几何估计作为一个任意固定 \(\eta\) 的事实，推出 gross 比值趋零。不在最终源价格中重复扣除 \(\eta Z\) 或 \(J\)。

## 8. 范围与计算

本条件命题只提纯高支撑 gross；不估计二轴的实际八阶源形状损失，也不把 R041 的具体下界系数认作最优。三近面的源层裁剪／重叠即使发生在二轴禁层之间，仍可能是八阶主量，不能因 (4) 删除；那是另一项实际源几何问题。

本稿复用已实际成功的 R043 branches.wl/log、R044 young_branches.wl/log/json 与 R046 固定层账本，不重跑已有有限代数。唯一新增的幂次插入全部逐项写在 (7)–(20)，可直接核查；计算程序不能替代目标实际分解、矩形截宽、原指定见证或统一重数。

条件依赖：

* R049：上述确切 SHA 的 (4)、(15)、固定稳健 \(\mathcal F_{\tau_*}\) 估计及正常层逐源体积；
* R043 正式源 SHA 75d0f89e79e8de4d00a6d3926532357f08a60a1ed5cbbc20b6abfd83e82ab4ac 的同正常槽、坏根界、完整面对路径 (20)；
* R044 正式源 SHA ef3a6df6ee4581e990fd55fbe76c2b03f43b2e5bdbd0cb27f7eba18231cf0cc4 的实际矩形、例外、逐根净价及全边界单源界。

状态：本稿给出 \((\mathrm H49)\Rightarrow(4)\) 的完整实际证明，待独立审核；在 R049 正式审完之前不取消条件。
