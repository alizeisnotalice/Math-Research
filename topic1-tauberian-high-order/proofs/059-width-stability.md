# R059：八阶竞争输入的实际端盒宽度均方根稳定性

2026-11-05。新命题及完整证明，供非作者审查。只写本独占 scratch；仓库、已发表论文和旧冻结稿只读。本稿没有提交、推送或修改主定理，**不据此宣布九阶全局渐近**。

## 1. 结论与准确输入

固定整数 (n\ge2)、固定 (B<\infty)，记

\[
a=v^n=2^nt,\quad A_v=(1+v)^n,\quad
P_7(t)=L_n(t)+K_nav^6-(n+26/9)K_nav^7.
\]

考虑任意实际有限正测度输入及外部严格 halo 紧集的序列

\[
v\downarrow0,\quad V=|E|>0,\quad M=|K|,\qquad
1+M/V\ge P_7(t)-Bav^8.                                      \tag{1}
\]

允许输入数、尺度跨度、原指定半径及指标数随 (v) 改变。使用现有正文的**同一个新八层族**：先固定原 (P,S) 图的 (I)、八层源/目标及全部深度常数，随后统一强 CDF；不以新结论重选源、半径、正常集合或 CDF。全篇 (S=\sum_i r_i^n) 仍是原有限族，未改成有效根的和。下列未注明层号的指标和均取本次第八层外部源；标记内部 (I) 不定义源盒。外部为空时所有结论按空和成立。

沿用

\[
g_i=aR_i^n,\quad \mu_i=|F_i|,\quad
b_i=|W_i|^{1/n},\quad
\varepsilon_i=|W_i|/\mu_i-1,\quad
\mathfrak q_i=\sum_{l=1}^n\log^2(w_{il}/b_i),
\]

其中 (w_{il}) 是中点盒 (W_i) 的边长。三个盒满足

\[
W_i'\subset W_i\subset W_i^0,
\quad c vR_i\le w_{il}^{\bullet}\le C vR_i,
\quad c g_i\le\mu_i\le g_i,
\quad R_i\asymp r_i,
\quad 0\le\varepsilon_i\le C.                              \tag{2}
\]

常数由已固定的维数、八层深度和 CDF 决定。实际还供应 (\mu_i/g_i\ge1-A_8v)，本证明只需其中统一正下界。三个盒的对应边长之差不超过总端移 (\Delta_i\)。定义与现正文完全相同的**最大相对误差**

\[
\xi_i=\max_{1\le l\le n,\,\bullet\in\{0,\prime,\mathrm{mid}\}}
        \left|\frac{w_{il}^{\bullet}}{vR_i}-1\right|.
                                                                  \tag{3}
\]

**命题。** 对上述任意序列，同一第八层有

\[
\boxed{
v\sum_i\mathfrak q_iR_i^n+
\frac{\mathcal D_8}{a}+v\frac{\mathcal D_7}{a}
\le C_{n,B}v^8S,}                                             \tag{4}
\]

并且

\[
\boxed{\quad v\sum_i\xi_i^2R_i^n\le C_{n,B}v^8S,\qquad
\sum_i\xi_iR_i^n=O_{n,B}(v^{7/2}S).\quad}                    \tag{5}
\]

所有构造常数先固定，随后取充分小 (v)。前缀 little-o 的模量可以依序列；常数无须依赖该模量，最终阈值允许依序列。没有逐源 (\xi_i=o(1)) 或统一逐源 (O(v^{7/2})) 假设，也不丢弃异常源。若已知任意固定八阶中心的竞争条件，它蕴含某个固定 (B) 的 (1)，故同样适用。

## 2. 实际可用接口及预算范围

本次实际读取的现仓库正文及 SHA256 为：

| 文件 | SHA256 |
|---|---|
| `paper/sections/global.tex` | `7b6154f74b56e4cb8ea133c99d4c240ceec526004930c510bc29ebcb25b5f10a` |
| `paper/sections/sharp.tex` | `8eae7fa41e8bc1f9883d51330b01183001295ca0b623da12205c807c6bbe7c39` |
| `paper/sections/eighth.tex` | `175b314a9f076c39583a8926921cdd8f99921b3579f8c2c9860df4221006f151` |
| `paper/main.tex` | `c04e39835f9c5b2bb3ae5081aea92a64807a8223a419bd14cb28bc5651275b4d` |

另实际读取 R049 `competitive_eighth_stability.md`，源 SHA `5c4d7b371a0c33c39449afd686de599c49570aaa639d4f7db7d3c1f166a040f4`，及 R045 正常层的实际定义。这里使用的供货只如下：

1. **同层真实质量。** 各层
   \[
   \mathcal D_h=\sum_i g_{i,h}-M'\ge f_h\ge0,\quad
   \mathfrak d_h^0\le f_{h-1}+aP_\partial,\quad f_0=A_0vM',
   \]
   且
   \[
   \sum_i(g_{i,h}-\mu_{i,h})
      =\mathcal D_h+\mathfrak d_h^0+\Delta_h^{\rm mass}.       \tag{6}
   \]
   第一层用 (\mathfrak d_1^0\le f_0)。

2. **真实基线的净重复量。** 在 Bonferroni、实际环境和删除处理后，每层仍有 (\mathcal D_h/(2a)) 的净正系数，只有未删除源缺量 (\mathfrak d_h^0) 乘固定 (b v/a)。对应标签为 `eighth:ledger`；系数 (b) 在辅助测试窗之前固定。

3. **同族前缀价格。** 以
   \[
   Q(v)=K_nv^6[(1+v)^{-4}+(10/9)v(1+v)^{-2}],\quad
   \mathcal T=V+M,
   \]
   第七、八层可取 (E_h=Q(v)\mathcal T+O(v^8S))。前六层 (E_h=O(v^6S)) 已足够本稿；实际第六层供货更强。只用 `eighth:two-prefix-prices` 和粗先验，不需要精确八阶价或第六层的 little-o 系数。

4. **原八阶稳定性的数值数据。** 记
   \[
   \mathscr D_i=(x_i+W_i^0)\setminus F_i,\quad
   Z_i=|\mathscr D_i\setminus Y_{i,N}|,
   \quad \mathfrak E=\sum_iR_i^{n-1}\Delta_i.
   \]
   原正常保证层满足 (Y_{i,N}\subset\mathscr D_i)，且
   \[
   R_N=(v/a)\sum_iZ_i=O(v^8S),\qquad
   \mathfrak E+\mathcal F+P_I+P_\partial=O(v^8S).            \tag{7}
   \]
   它们是 R049／`eighth:exact-stability` 已供的同层数据，未在新证明中将它们重新当收入相加。

5. **正常层逐源真实体积。**
   \[
   |Y_{i,N}|\le C a v^5R_i^n.                               \tag{8}
   \]
   这是原正常有界窗的几何上界，不是逐源 raw 小量：每个正常目标的两个临界交长 (O(v^3R_i))，实际正高度至多 (Cv^4R_i)，排序底变量宽 (Cv^3R_i)，其余 (n-2) 个源变量宽 (CvR_i)。有限相交度给总量 (C v^{n+5}R_i^n)。只保留真实目标质量、减非负侵蚀和剪核心都会减小该层。非正常源或边界源的 (Y_{i,N}) 可为空，(8) 仍成立。这里用的是已锚定的旧正常保证层，不以新双反馈模型代替它。

6. **原小源删除的数值上界。** 八阶粗证明给
   \[
   (v/a)\Delta_{8,\mathrm{in}}^{\rm mass}
   \le Ca v\sum_i\varepsilon_i r_i^n+Cv^n\mathcal F,
   \quad (v/a)\Delta_{8,\mathrm{out}}^{\rm mass}\le Cv\mathcal F,
   \]
   边界项至多 (CavS(\partial I))。保留旧
   (\sum\varepsilon_i r_i^n=O(v^5S))、(S(\partial I)=O(v^5S))，
   再用 (7)，即得
   \[
   \frac{\Delta_8^{\rm mass}}a
       \le C(v^{n+5}+v^8)S=O(v^7S).                         \tag{9}
   \]
   最坏维数为 (n=2)。无需调用二维 sharp inside 的准确常数；原全半径 dyadic 上界已经足够。环境不是删除质量，在 (9) 中没有混入环境费用。

上述接口确实适用于同一个实际 N8 及其第七前缀。不得以已发表状态免除这个对应检查。
事实上固定 (B) 后，最终 (Bv\le1)，(1) 蕴含现正文 `eighth:critical` 的固定 (B_0) 分支，故三个前缀和粗先验的假设确实满足。

## 3. 独立于 (J_i) 的各向异性实际盈余

设

\[
f_i(m)=(R_i+m^{1/n})^n,\quad
J_i=f_i(|W_i|)-f_i(\mu_i),\quad
\mathcal A_i=|X_i|-f_i(|W_i|).
\]

由 (X_i=W_i+[-R_i/2,R_i/2]^n) 的实际边长，准确有

\[
|X_i|=\prod_l(R_i+w_{il})
      =f_i(\mu_i)+J_i+\mathcal A_i.                         \tag{10}
\]

正常共同价格、反馈和补集消费的是 (J_i)，不是 (\mathcal A_i)。现证明对完整实际边长的统一定量下界

\[
\mathcal A_i\ge c v\mathfrak q_iR_i^n\ge0.                 \tag{11}
\]

令 (s_l=\log(w_{il}/b_i))，则 (\sum_l s_l=0)，(s_l) 在固定紧区间。置 (\kappa=b_i/R_i\asymp v)，

\[
F(s)=\sum_l\log(1+\kappa e^{s_l}),\qquad
\partial_{ll}F(s)=\frac{\kappa e^{s_l}}{(1+\kappa e^{s_l})^2}\ge cv.
\]

整条 (t s\)、(0\le t\le1) 仍在该紧区间；在零点梯度各分量相同，因此一阶项为零。积分 Taylor 公式给 (F(s)-F(0)\ge cv\sum_l s_l^2)。再用 (e^u-1\ge u\) 及 ((1+\kappa)^n\ge1)，得

\[
\begin{aligned}
\mathcal A_i/R_i^n
 &=e^{F(s)}-e^{F(0)}\\
 &=(1+\kappa)^n\{e^{F(s)-F(0)}-1\}
 \ge cv\sum_l s_l^2.
\end{aligned}
\]

此证明没有在逐源小偏差上作 Taylor 近似；紧区间上保留的 Hessian 下界适用于全部源，包括罕见强各向异性源。

## 4. 同一次实际账本同时保留形状与前层重复量

设 (H_h=\sum_{i<j}|X_i\cap X_j|-\sum_iJ_i)。对实际并集作二阶 Bonferroni，并代入 (10)，得到

\[
\left|\bigcup_iX_i\right|
\ge\sum_i f_i(\mu_i)-H_h+\sum_i\mathcal A_i.                \tag{12}
\]

这是在原 `eighth:ledger` 中丢弃非负 AM–GM 差之前的准确公式。原同层 joint/high 价格只上包 (H_h)，原删除成本按已证数值上界放进 (E_h)，环境扣实际重复量。故重新保留 (12) 最后项后，同一实际账本可以写成

\[
\begin{aligned}
\mathcal T\ge{}&A_vM/a+c_{\rm rec}P_I-C_{\rm b}P_\partial-2nvS_I\\
&+\tfrac12\mathcal D_h/a-b(v/a)\mathfrak d_h^0-E_h
   +\mathbf1_{\{h=8\}}\sum_i\mathcal A_i.                 \tag{13}
\end{aligned}
\]

早层的形状差只丢弃，不需另外升级。保留第八层形状没有拆走 (J_i) 的任何份额，也没有重新购买 raw、面条、删除或非正常目标价格。已有粗宽度先验在上包 (H_h) 的误差中可作为数值数据使用；那不等于再次收入它所曾用于证明该先验的正预算。

现在把正系数 (1/2) 拆成两份 (c_0=1/4)：

\[
\tfrac12\mathcal D_h/a-b(v/a)\mathfrak d_h^0
\ge \frac{c_0f_h-bvf_{h-1}}a
        +c_0\mathcal D_h/a-bvP_\partial.                  \tag{14}
\]

第一层也可用右端较弱式。令

\[
\tau=bv/c_0=4bv,\quad w_h=\tau^{8-h},\quad
\Omega=\sum_{h=1}^8w_h=1+O(v).
\]

这只改变正权平均的数值系数，不改变任何源、深度、CDF 或图。中间缺量精确消去：

\[
\sum_{h=1}^8w_h(c_0f_h-bvf_{h-1})
      =c_0f_8-bv\tau^7f_0.                                \tag{15}
\]

初端除 (a\Omega) 为 (O(v^9S))。未用于消去的另一份重复量保留为

\[
\frac{c_0}{\Omega}\sum_hw_h\mathcal D_h/a
\ge c\left(\mathcal D_8/a+v\mathcal D_7/a\right).           \tag{16}
\]

同时末层形状保留 (\Omega^{-1}\sum\mathcal A_i\ge c v\sum\mathfrak q_iR_i^n)。这里不可把 (\mathcal D_7) 全部换成 (f_7) 然后又声称保留它；拆分 (14) 正是必要的新一步。

误差也以同一个中心价平均。由 §2 的前缀价格和 (\mathcal T\asymp S)，

\[
\begin{aligned}
\bar E-Q\mathcal T
={}&\Omega^{-1}\bigl[(E_8-Q\mathcal T)
  +\tau(E_7-Q\mathcal T)
  +\sum_{h\le6}w_h(E_h-Q\mathcal T)\bigr]\\
={}&O(v^8S).                                               \tag{17}
\end{aligned}
\]

最后一和从第六层的 (O(v^2)v^6S=O(v^8S)) 开始；不需将第六层真实领先误差误说成统一 little-o。由于 (b) 固定，最后小 (v) 时 (bv\le1)，(14) 新写出的 (bvP_\partial) 落在原 (C_*-C_{\rm b}=2) 的纯维数边界预留内；原图条件仍吸收固定边界和 (2nvS_I)。同一内部收入只出现一次。

由 (13)–(17) 得

\[
(1+Q)\mathcal T-A_vM/a
\ge c\left[v\sum_i\mathfrak q_iR_i^n+
             \mathcal D_8/a+v\mathcal D_7/a\right]-Cv^8S.  \tag{18}
\]

这一侧没有将原 R049 已有正量再加作第二份收入。若需，其 (J)/面条固定保留份额也可原封保留；本稿只利用 (7) 的数值结论。

另一侧只用真实比值。令

\[
\mathcal R_Q(v)=\frac{A_v}{A_v-a(1+Q(v))}.
\]

准确恒等式

\[
(1+Q)\mathcal T-A_vM/a
=\frac{A_v-a(1+Q)}a V[\mathcal R_Q(v)-(1+M/V)]             \tag{19}
\]

及 (\mathcal R_Q=P_7+O_n(av^8))、(1)、(V\asymp S)，给公共左侧至多 (C_{n,B}v^8S)。分母在小参数时正且统一远离零。与 (18) 比较即证 (4)。此处没有把 (aS) 代为 (M)，也没有从某个特殊输入的价格反推全部输入上确界。

## 5. 等向质量偏差：须控制平方，不可重留盈余的一阶量

定义准确根归一化源缺陷

\[
\eta_i=1-\mu_i/g_i\ge0.
\]

这与旧 (\theta_i=(g_i-\mu_i)/(ar_i^n)) 只差有界根比；精确关系为

\[
\sum_i\eta_iR_i^n
=\frac{\mathcal D_8+\mathfrak d_8^0+\Delta_8^{\rm mass}}a.
                                                                  \tag{20}
\]

由 (4)、(\mathfrak d_8^0\le f_7+aP_\partial\le\mathcal D_7+aP_\partial)、(7)、(9)，得

\[
\sum_i\eta_iR_i^n=O(v^7S),\qquad
\sum_i\eta_i^2R_i^n=O(v^7S),                                \tag{21}
\]

后式只需 (0\le\eta_i\le1)。因此前层重复量不能遗漏；仅升级 (\mathfrak q_i) 不足以证明三盒 (\xi_i) 的结论。

对 (\varepsilon_i)，**不能**声称重新保留 (v\sum\varepsilon_iR_i^n=O(v^8S))：正常禁层本身已是六阶 (J) 成本的主量。正确的平方界从原真实缺集直接得出。因 (W_i\subset W_i^0)，

\[
0\le\varepsilon_i\mu_i=|W_i|-\mu_i
\le|W_i^0|-\mu_i=|\mathscr D_i|
\le C av^5R_i^n+Z_i.                                     \tag{22}
\]

令 (z_i=Z_i/(aR_i^n))。由 (2)，(0\le z_i\le C)，且 (7) 给

\[
\sum_i z_iR_i^n=\frac1a\sum_iZ_i=R_N/v=O(v^7S).
\]

于是 (\varepsilon_i\le C(v^5+z_i))，并且

\[
\sum_i\varepsilon_i^2R_i^n
\le C v^{10}\sum_iR_i^n+C\sum_i z_i^2R_i^n
=O(v^7S).                                                   \tag{23}
\]

这里保留所有 raw 高度异常；没有从加权平均跳到逐源小量。正常层逐源的 (v^5) 在平方后是 (v^{10})，而 raw 只用统一有界和全质量。没有第二次花费 (J_i)。

## 6. 三种端盒和最大宽度定义

首先中点盒的等向尺度准确为

\[
\frac{b_i}{vR_i}=[(1+\varepsilon_i)(1-\eta_i)]^{1/n}.       \tag{24}
\]

其根号内变量在固定正紧区间，故实函数 (x^{1/n}) 的导数一致有界；由 (24)

\[
\left|b_i/(vR_i)-1\right|^2
\le C(\varepsilon_i^2+\eta_i^2).                           \tag{25}
\]

又 (w_{il}=b_i e^{s_{il}})，(s_{il}) 一致有界，(b_i/(vR_i)) 亦一致有界，因而

\[
\max_l|w_{il}/(vR_i)-1|^2
\le C[\mathfrak q_i+\varepsilon_i^2+\eta_i^2].              \tag{26}
\]

对端移令 (\rho_i=\Delta_i/(vR_i))。原 CDF 给 (0\le\rho_i\le C)，(7) 给

\[
\sum_i\rho_iR_i^n=\mathfrak E/v=O(v^7S),\qquad
\sum_i\rho_i^2R_i^n=O(v^7S).                               \tag{27}
\]

三个盒对应边长差至多 (\Delta_i)，所以**按 (3) 的同一个最大值**有

\[
\xi_i^2\le C[\mathfrak q_i+\varepsilon_i^2+\eta_i^2+\rho_i^2].
                                                                  \tag{28}
\]

将 (4)、(21)、(23)、(27) 相加得到 (\sum\xi_i^2R_i^n=O(v^7S))，即 (5) 第一式。又 (R_i\le r_i) 给 (\sum_iR_i^n\le S)，Cauchy–Schwarz 给

\[
\sum_i\xi_iR_i^n
\le(\sum_i\xi_i^2R_i^n)^{1/2}(\sum_iR_i^n)^{1/2}
=O(v^{7/2}S),
\]

证明完整命题。所有指标集、包括 raw 集中的罕见异常块和非局部尺度，均在原有限和内。

## 7. 实际改善的接口与未作的升级

现正文 `eighth:exact-normal-price` 的正常名义成本宽度替换误差为

\[
Cv^6\sum_i\xi_iR_i^n.
\]

用 (5) 后它从 (O(v^{17/2}S)) 提高到

\[
\boxed{\quad O(v^{19/2}S)=o(v^9S).\quad}                   \tag{29}
\]

正常相交的有效根可比且次数有界，故逐边 \((\xi_i+\xi_j)R_e^n\) 的和也满足同一界；不需逐边归一化宽差小。该结论精确消除了所问的九阶宽度障碍。

这不消除所有九阶项。实际目标锚移的退款 (Cv\mathfrak E=O(v^9S))、真实纤维顶截损失、三近面退款、补集/内部删除及更深前缀的同阶项，仍须各自重算或保留，不得由 (29) 自动放进 (o(v^9S))。本稿不识别九阶系数，也不把本地宽度命题升格成原 (C_n) 的九阶上界。

**有限深度迁移的条件口径。** 本文无条件结论只针对已列明的同 N8。若另一个预先固定的有限 (N\ge8) 实际族已经独立供应：同一个公共图账本及 (\mathcal D_h\ge f_h)、(\mathfrak d_h^0\le f_{h-1}+aP_\partial)；末两前缀 (E_N,E_{N-1}=Q\mathcal T+O(v^8S))，其余前缀 (O(v^6S))；本族末层 (2)、(7)–(9) 和正常层逐源界 (8)，则本证明按 (8\mapsto N\)、(7\mapsto N-1) 重放，端点成为 (O(v^{N+1}S))，得到本族的 (4)–(5)。这只列迁移的充分条件，没有证明新 N9 实际构造具有这些供货，也没有把旧 N8 对象重命名为新前缀。

## 8. 实际本机核与证明边界

本轮将确切局部 Hessian 下界交 Terra 在本机核验，并实际读取成功代码与 stdout：

| 材料 | SHA256 |
|---|---|
| `r059-computation/width_amgm.wl` | `3457b2a4d32942e68522debb4aaf9284d4c44aacef451e2637ed5c1584ba60a4` |
| `r059-computation/width_amgm.log` | `4e2e73f3192d16996817b4ae076ea38b183a20ffdbe815e7c2dbf2ef1705628b` |
| `r059-computation/width_fraction.py` | `85f7c5f4f69771e1a1634b13416b6470c8faefa53cfca806a84e4c7a9d7fe25d` |
| `r059-computation/width_fraction.json` | `46f197b6df16a4742a1b3701ac87b6303b1f23b983c7ada053eafb7038d1651d` |

成功 WL 日志为 Wolfram 15.0 ARM，PID 32488。Hessian 公式残差 0。将域写成 (0\le lo\le xx\le up)，交叉乘积不等式

\[
xx(1+up)^2\ge lo(1+xx)^2
\]

的违反集为 False，其精确非负分解为

\[
(xx-lo)(1+xx)^2+xx(up-xx)(2+up+xx).
\]

代 (lo=c c_yv\)、(xx=\kappa e^s\)、(up=C C_yv_0\) 即是 §3 所用统一 Hessian 下界。二维准确 AM–GM 差 (2\kappa(\cosh s-1)) 的残差也为 0。三组独立 Fraction 输入（包括等号边界）的裕量为 (0,275/2304,359/1600)，均通过；它们只检查有理抄录，不外推一般维数。

首轮曾误用内置保护符号 C，作者保留 `first-warning.log`，最终已更换符号重跑无警告；本稿不以该首轮 warning 输出为成功证据。此计算只认证局部标量恒等与正分母；全维实际 AM–GM 差、同时保留预算、全部输入和尺度的证明是 §§1–7，而非有限计算外推。

**本轮冻结内容：一般固定 (n\ge2) 的 (4)–(5) 及正常宽度误差的 (29)。待非作者独立审查；没有修改已发表主稿，也没有声称九阶整体升级。**
