# R030：二维七层账本的条件性次阶结算

日期：2026-11-05。仓库只读；本稿只处理总装。已经使用的真实输入为 R023 六阶锐定理及其全部输入证明、R027 的二维实际半径 incidence。R028 的局部七阶价格在本稿中始终是条件，不据其当前审核状态宣布无条件七阶渐近。

## 1. 结论及唯一待供局部接口

记
\[
 n=2,\quad a=v^2=4t,\quad A_v=(1+v)^2,\quad k=\frac4{27},
 \quad L_2(t)=\frac{A_v}{A_v-a}.
\]
给定任意有限正测度输入 \(E\)，令 \(V=|E|\)，取外部严格 halo 内的紧集 \(K\)，令 \(M=|K|\)、\(\mathcal V=V+M\)。原有限严格见证贪心给 \(Q_i=Q(x_i,r_i)\)、两两不交的半盒 \(B_i\)，写 \(S=\sum_i r_i^2\)。只研究任意满足
\[
 v\downarrow0,\qquad
 1+\frac MV\ge L_2(t)+kav^6-Bav^7,
 \qquad B<\infty\text{ 固定}                                      \tag{1.1}
\]
的输入序列。下文小量均沿任意这样的序列成立，收敛模量允许依赖序列；常数不依赖输入、指标数或尺度跨度。

在 §2 重新构造的**同一七层族**第七层，记
\[
 R_i=(g_{i,7}/a)^{1/2},\qquad
 X_i=x_i+W_{i,7}+[-R_i/2,R_i/2]^2,
\]
\[
 J_i=(R_i+|W_{i,7}|^{1/2})^2-(R_i+\mu_{i,7}^{1/2})^2,
 \qquad H_7=\sum_{i<j}|X_i\cap X_j|-\sum_iJ_i.
\]
\(J_i\) 是完整实际源盈余；pair 和包含全部外部指标，特别包括边界和异常端点。对深处正交二临界边 \(e=\{i,j\}\)，取 \(R_e=\max(R_i,R_j)\)，令 \(W_7=\sum_eR_e^2\)。

**条件 \(\mathrm H_7\)。** 对每个 (1.1) 序列及上述重建第七层，有
\[
 H_7\le \frac{k}{A_v}v^6W_7+\frac{10}{9}k v^7W_7+o(v^7S).            \tag{1.2}
\]
这里要求的是该层的实际集合，不是先前已经构造的六层末层。R028 当前原稿将其陈述于 R023 六层末层；要向本稿供货，必须确认其证明对下述固定有限七深度递归也成立。所需粗先验、六阶饱和及免费面严格小量均在本稿重建，不另作局部假设。

**条件总装定理。** 若 \(\mathrm H_7\) 成立，则
\[
 \limsup_{t\downarrow0}
 \frac{C_2(1-t)-L_2(t)-\frac4{27}av^6}{av^7}
 \le-\frac{176}{243}.                                               \tag{1.3}
\]
这与 R026 的真实有限输入下界恰好相等。因此在 \(\mathrm H_7\) 成立条件下，
\[
 C_2(1-t)=L_2(t)+\frac4{27}av^6-\frac{176}{243}av^7+o(av^7).           \tag{1.4}
\]
本稿不把 (1.2) 的未完成独审替换成既成前提；已经无条件证明的仍包括 R023 六阶渐近和 R026 的这一侧七阶下界。

证明的两个精度提升是：给第六层保留已证的领先价格，使它的一阶权重与权重和的归一化相消；把小源回收的原始系数 \(a\) 留在 dyadic 估计中，再用免费面饱和排除删除量在七阶的残留。最后以实际 \(V+M\) 表示角权，避免以 \(S\) 的粗分母损失次阶系数。

## 2. 同一个图标记与七个固定深度

先固定 R014/R023 的纯维数价格
\[
 H=4,\quad c_{\rm rec}=\frac18,\quad
 C_{\rm b}=2C_2(2)+2,\quad C_*=C_{\rm b}+2.
\]
这里 \(C_2(2)\) 是原内部回收的几何常数，与 Tauberian 函数的记号区分。原完整见证图的 \(K_0,\eta_{\rm w},D\) 仍按 `global.tex` 的 `eq:global-constant-order` 选定。标记只用原 \(P,S\)：在有限图层 \(J_k\) 中取第一个满足原加权边界条件的层为 \(I\)。它与后续源盒、CDF 斜率、深度及局部角阈值均无关，并满足
\[
 C_*P_\partial\le\frac{c_{\rm rec}}4P_I,
 \qquad 4vS_I\le\frac{c_{\rm rec}}4P_I.                              \tag{2.1}
\]
此处 \(P_I=P(I)\)、\(P_\partial=P(\partial I)\)、\(S_I=S(I)\)。同一
\[
 K'=K\setminus\bigcup_{i\in I}Q_i,\quad M'=|K'|,\quad M_I=M-M'
\]
用于全部深度。

固定这些纯几何参数后，按 R014 `lem:global-exterior-union` 的递归一次执行至 \(h=7\)。已固定 \(L_{h-1}\) 后先取深目标常数，再取中心邻域常数 \(L_h\)。深处未删除源为前层中心邻域，边界未删除源仍为原私有 \(P_i\)；删除的始终是原指定严格见证半径低于该层固定阈值的点。所有大见证、源与其他同层目标不交的条件逐层保持。统一全部常数为 \(A_7\) 后才选正的强 CDF 斜率，最后令 \(v\) 充分小。常数的有限次增加不改变 (2.1) 的纯维数价格；源常数产生的 \(CvP_\partial\) 用随后的小 \(v\) 吸收。

对未删除源 \(T_{i,h}\)、有效源 \(F_{i,h}\)、目标 \(G_{i,h}\)，设
\[
 m^0_{i,h}=|T_{i,h}|,\quad \mu_{i,h}=|F_{i,h}|,\quad
 g_{i,h}=|G_{i,h}|,\quad
 \delta_{i,h}^{\rm mass}=m^0_{i,h}-\mu_{i,h},\quad
 \Delta_h^{\rm mass}=\sum_i\delta_{i,h}^{\rm mass},
\]
\[
 \mathcal D_h=\sum_i g_{i,h}-M',\qquad
 \mathfrak d_h^0=M'-\sum_i m^0_{i,h}.
\]
精确恒等式和源层递推为
\[
 \sum_i(g_{i,h}-\mu_{i,h})
 =\mathcal D_h+\mathfrak d_h^0+\Delta_h^{\rm mass},                    \tag{2.2}
\]
\[
 \mathcal D_h\ge f_h\ge0,\qquad
 \mathfrak d_1^0\le A_0vM',\qquad
 \mathfrak d_h^0\le f_{h-1}+aP_\partial\quad(2\le h\le7),             \tag{2.3}
\]
其中 \(f_h=M'-\sum_i|N_{i,h}|\)、\(N_{i,h}=K'\cap Q(x_i,L_hvr_i)\)。未删除源两两不交，所以 \(\mathfrak d_h^0\ge0\)。小源删除量不在 \(\mathfrak d_h^0\) 中。令 \(f_0=A_0vM'\)。

## 3. 同时得到第六层和第七层粗临界先验

关键是给已经一次构造的七层族分别使用长度六、长度七的前缀；不重新选 \(I\)，也不更换第六层集合。R023 `sharp-coarse-retained-recursion` 的保留正量递推在全部七层给
\[
 Z_0\ge c_0 f_h-b_0vf_{h-1}+a\mathcal J_h+\mathcal D_h/4
          -Ca(v^6S+vP_\partial),                                   \tag{3.1}
\]
\[
 Z_0=a[\mathcal V-A_vM/a-c_{\rm rec}P_I+C_{\rm b}P_\partial+4vS_I].
\]
这里 \(\mathcal J_h\ge0\) 包括固定正比例的
\[
 v\sum_i(\varepsilon_{i,h}+\mathfrak q_{i,h})r_i^2,
 \qquad \mathcal F_h,\qquad c_\tau\mathcal F_{\tau,h},
\]
分别是源过量、形状余量、原免费面量和固定稳健裕量的免费面量。\(c_0,b_0,C\) 可在所有七深度常数固定后统一选择。

分别以 \((b_0v/c_0)^{m-h}\)、\(1\le h\le m\) 加权 (3.1)，其中 \(m=6,7\)。中间 \(f_h\) 精确抵消；分开保留原重复正量的一部分，就同时保留 \(\mathcal D_m/a\)、\(v\mathcal D_{m-1}/a\) 和 \(\mathcal J_m\) 的固定正比例。初端除以 \(a\) 分别为 \(O(v^7S)\)、\(O(v^8S)\)。图项用同一个 (2.1) 吸收。

(1.1) 给
\[
 \Xi:=\mathcal V-A_vM/a
 =V-(A_v-a)M/a=O(v^6S)\quad\text{的上界}.
\]
下界由已证全输入粗余项给出，故 \(|\Xi|=O(v^6S)\)。于是，对 \(m=6,7\) 各有
\[
 P_I+\mathcal D_m/a+
 v\sum_i(\varepsilon_{i,m}+\mathfrak q_{i,m})r_i^2
 +\mathcal F_m+c_\tau\mathcal F_{\tau,m}=O(v^6S),                      \tag{3.2}
\]
\[
 \mathcal D_{m-1}/a+\mathfrak d_m^0/a+S_I+S(\partial I)=O(v^5S),
 \qquad \Delta_m^{\rm mass}/a=O(v^6S).                              \tag{3.3}
\]
盒宽紧性、端移和稳健面对根差也按 R023 `sharp-width-tightness`、`sharp-face-tightness` 逐层成立。它们只用 (3.2)、(3.3)，没有六这一特定深度的性质。

因此 R019 的原输入全尺度 LT 提取和 R017/R023 的领先 incidence，可分别用于本次第六层、第七层；所有坏根、边界及异常端点仍按原全部尺度口径处理。R023 已证完整鲜盈余价格给
\[
 H_m:=\sum_{i<j}|X_{i,m}\cap X_{j,m}|-\sum_iJ_{i,m}
 \le[k+o(1)]v^6S,\qquad m=6,7.                                    \tag{3.4}
\]
这里先固定旧证明的盈余份额和角阈值，再取 \(v\to0\)，最后令份额趋零，故 (3.4) 是任意输入序列的上极限陈述，不要求给定收敛率。它仅是已经证明的六阶价格，未使用 (1.2)。

另由原精确 \(M/a\le S\le V/(1-a)\)、(1.1) 及全输入粗上界，
\[
 M/V=a[1+O(v)],\qquad S/V=1+O(v),\qquad S/\mathcal V=1+O(v).          \tag{3.5}
\]
这些关系只在需要领先归一化时使用，没有以 \(aS\) 代替实际 \(M\)。

## 4. 末层免费面严格小量由领先饱和独立推出

先不用 (1.2)。在第七层把 R014 `lem:global-free-face-strips` 的原免费面条加回 \(\bigcup X_{i,7}\)。这些条与全部外部 \(X_{i,7}\) 及标记区域不交；有界重数给出至少 \(c\mathcal F_7\) 的新体积。共同环境的附加缺陷仍为 \(Ca(\mathcal F_7+P_\partial)\)，小源成本按 (3.3) 为 \(O(v^7S)\)，故在取小 \(v\) 后，固定正比例的 \(\mathcal F_7\) 可以保留。

保留精确的有符号 \(H_7\)，而不是先把它替换成上界。末层单层式因而可用 \(H_7-c\mathcal F_7+O(v^7S)\) 为误差；前六层仍用 \(Cv^6S\) 的粗误差。具体地，取固定足够大的公共基线常数 \(b\)，令 \(\tau=bv/c_D\)、\(c_D=1/2\)、\(\Omega=1+\tau+\cdots+\tau^6\)，在各层公共账本中使用
\[
 \widehat E_h=C_0v^6S\ (h\le6),\qquad
 \widehat E_7=H_7-c\mathcal F_7+C_1v^7S.
\]
误差 \(\widehat E_7\) 可带符号，这不影响同左侧不等式的正权相加。精确恒等式
\(\sum_{h=1}^7\tau^{7-h}(c_Df_h-bvf_{h-1})=c_Df_7-bv\tau^6f_0\)
消去中间缺量。\(\Omega-1=O(v)\)、\(M'\le aS\) 使早六层粗成本为 \(O(v^7S)\)、初端为 \(O(v^8S)\)。边界图项按 (2.1) 吸收，得到
\[
 \mathcal V\ge A_vM/a-
       \frac{H_7-c\mathcal F_7}{1+O(v)}-O(v^7S).                     \tag{4.1}
\]
此处只使用已证的粗账本、免费面排他和实际完整 \(J_i\)，未使用任何七阶局部锋利性。常数归一化后的 \(c>0\) 固定。图项仍按 (2.1) 丢弃其非负剩余。

(1.1) 给
\[
 A_vM/a-\mathcal V
 \ge[k+O(v)]v^6V.
\]
又 \(V\ge(1-a)S\)，故 (4.1) 推出
\[
 H_7-c\mathcal F_7\ge[k-o(1)]v^6S.
\]
与 (3.4) 的 \(m=7\) 情形合并，得到
\[
 \boxed{H_7/(v^6S)\to k,\qquad \mathcal F_7=o(v^6S).}               \tag{4.2}
\]
这独立核实了 R028(3) 中本次总装所需的免费面结论，且适用于重建第七层。无需从待审的 (1.2) 反推它。

## 5. 真实小源删除成本是 \(o(v^7S)\)

以下全在第七层。对深处 \(i\)，按 R014 小源回收证明构造扩展盒 \(\mathcal R_i\)，设 \(D_i=\mathcal R_i\setminus T_{i,7}\)。其边长可比于 \(vr_i\)，扩展量由入向面缝 \(\gamma_{il}^{\pm}\le2\Gamma_{ip}\) 控制。

盒外部分用原入向邻居的同一个严格见证及剩余容量，得到
\[
 \delta_{i,\mathrm{out}}^{\rm mass}/a
 \le C r_i\sum_{l,\pm}\gamma_{il}^{\pm}.                              \tag{5.1}
\]
这正是 `eq:global-incoming-small-tail` 逐半空间相加后的结论；端点不必取到，半空间质量连续性足够。

盒内部分必须保留 dyadic 引理的原始系数 \(a\)，而不是最后粗化为 \(v\)。原指定小见证满足 \(r(x)\le A_7vr_i\)，故
\[
 \delta_{i,\mathrm{in}}^{\rm mass}\le Ca|D_i|,
 \quad |D_i|\le C\varepsilon_i\mu_i
               +Cvr_i\sum_{l,\pm}\gamma_{il}^{\pm}.
\]
用 \(\mu_i\le Ca r_i^2\) 并求和，得到
\[
 \sum_{i\ \mathrm{deep}}\delta_{i,\mathrm{in}}^{\rm mass}/a
 \le C\left[a\sum_i\varepsilon_i r_i^2+v\mathcal F_7\right],
 \quad
 \sum_{i\ \mathrm{deep}}\delta_{i,\mathrm{out}}^{\rm mass}/a
 \le C\mathcal F_7.                                                  \tag{5.2}
\]
相同面对的尺度可比、重数有界，因而入向 \(\gamma\) 的和确实由原 \(\mathcal F_7\) 控制。

边界源不套深处面邻居论证。由原局部小源估计及 \(S(\partial I)=O(v^5S)\)，
\[
 \sum_{i\in\partial I}\delta_i^{\rm mass}/a
 \le CaS(\partial I)=O(v^7S).                                       \tag{5.3}
\]
等价地可直接用 `eq:global-small-source-boundary` 得 \(CvP_\partial\)。

由 (3.2)、(4.2)、\(a=v^2\)，(5.2) 的第一项为 \(O(v^7S)\)，第二项为 \(o(v^6S)\)。从而
\[
 \boxed{\Delta_7^{\rm mass}/a=o(v^6S),\qquad
        (v/a)\Delta_7^{\rm mass}=o(v^7S).}                            \tag{5.4}
\]
这里没有把真实删除量等同于 CDF 端移。即使盒内缺量为零，dyadic 引理也给盒内小源零质量，故无零测例外。

共同环境另有
\[
 \sum_i u_{i,7}\le Ca(\mathcal F_7+P_\partial)=O(v^8S)=o(v^7S).         \tag{5.5}
\]
因此环境项不需要额外七阶局部假设。边界的固定突出价格保留给图回收，不能一起称作小量。

## 6. 第七层的真实体积价格与实际半径权

取 \(U_7=\bigcup_{i\notin I}X_i\)，不加 §4 的免费面条。Bonferroni、完整 \(J_i\) 及条件 (1.2) 给
\[
 |U_7|\ge\sum_i(R_i+\mu_i^{1/2})^2
       -\frac{k}{A_v}v^6W_7-\frac{10}{9}kv^7W_7-o(v^7S).             \tag{6.1}
\]
标量基线只用一个固定粗常数：
\[
 (R_i+\sqrt{\mu_i})^2
 \ge A_vg_i/a-C(v/a)(g_i-\mu_i).
\]
代入 (2.2)，扣除物理 \(\mathcal D_7\) 和环境缺陷 (5.5)。重复量净系数 \((A_v-Cv-a)/a\ge1/(2a)\)，删除项由 (5.4) 为 \(o(v^7S)\)，故
\[
 |U_7|-\mathcal D_7-\sum_i u_{i,7}
 \ge A_vM'/a+\mathcal D_7/(2a)-b(v/a)\mathfrak d_7^0-E_7,             \tag{6.2}
\]
\[
 E_7=\frac{k}{A_v}v^6W_7+\frac{10}{9}kv^7W_7+o(v^7S).
\]
\(b>0\) 固定，可增大以统一全部七层。\(\mathfrak d_7^0/a=O(v^5S)\) 的成本仍可能是六阶，故必须保留到望远镜，不能列入七阶误差。

R027 的实际半径 incidence 适用于本层，因为 §3 已给 \(\mathcal D_7/a=O(v^6S)\)、\(|\Xi|=O(v^6S)\)、宽度紧性、稳健面根差和 R019 的全尺度 LT。其全角边和包含本稿的深处边集，所以对 \(W_7\) 同样成立；所有角权采用本次 \(R_e\)，并没有删除大角：
\[
 W_7\le\frac{\mathcal V}{A_v}+O(v^{9/2}S).                            \tag{6.3}
\]
乘以 \(v^6\) 或 \(v^7\) 后的余项均为 \(o(v^7S)\)。所以可在 (6.2) 中取非负上包
\[
 E_7=Q(v)\mathcal V+e_7,\quad e_7=o(v^7S),\qquad
 Q(v)=\frac{k}{A_v^2}v^6+\frac{10}{9}\frac{k}{A_v}v^7.                \tag{6.4}
\]
必要时增大 \(e_7\) 为其正部。第一价格和 incidence 各给一个 \(A_v^{-1}\)，这两个因子都保留到解分母。

第六层使用同一七层族的前缀及 (3.4)，原删除项只需 \(O(v^7S)\)。由 R021/R023 已证的领先末层结算，可取
\[
 E_6=[k+o(1)]v^6S.                                                  \tag{6.5}
\]
前五层取 \(E_h=C_0v^6S\)。第六层不要求局部七阶价格；(6.5) 在一阶权重下已经足够。

## 7. 公共左侧、七层精确抵消与第六层归一化

对前五层用粗增广并集，第六层用领先鲜盈余并集，第七层用 \(U_7\)。各层都加入同一标记半盒和内部突出集。R014 `eq:global-recovered-union`--`eq:global-exact-interior-identity` 的内部回收逐层给同一个公共左侧，输入体积只计一次：
\[
 \begin{split}
 \mathcal V\ge{}&A_vM/a+c_{\rm rec}P_I-C_{\rm b}P_\partial-4vS_I\\
 &+c_D\mathcal D_h/a-b(v/a)\mathfrak d_h^0-E_h,
 \qquad c_D=1/2,\quad 1\le h\le7.                                  \tag{7.1}
 \end{split}
\]
纯维数边界预留吸收各层 \(CvP_\partial\)，所需小 \(v\) 可依赖已经固定的 \(A_7\)。完整 \(J_{i,6}\)、\(J_{i,7}\) 各自用于不同单层不等式；跨层加权的是这些同左侧不等式，没有在物理上重复消费一个源盈余。

置
\[
 Z=a[\mathcal V-A_vM/a-c_{\rm rec}P_I+C_{\rm b}P_\partial+4vS_I],
 \quad \tau=bv/c_D,\quad
 w_h=\tau^{7-h},\quad \Omega=\sum_{h=1}^7w_h=1+\tau+\cdots+\tau^6.
\]
由 (2.3)、(7.1)，
\[
 Z\ge c_Df_h-bvf_{h-1}-aE_h-abvP_\partial.                            \tag{7.2}
\]
精确望远镜为
\[
 \sum_{h=1}^7w_h(c_Df_h-bvf_{h-1})
 =c_Df_7-bv\tau^6f_0.                                               \tag{7.3}
\]
特别是 \(-bvf_6\) 与第六层的一阶加权正重复量完全抵消。前端为
\[
 \frac{bv\tau^6f_0}{a\Omega}
 =\frac{A_0b^7}{c_D^6\Omega}v^8\frac{M'}a=O(v^8S),                   \tag{7.4}
\]
而非七阶项。

令 \(\bar E=\Omega^{-1}\sum_h w_hE_h\)。不能将第六层简单计作 \(O(v^7S)\) 后丢弃其系数。由 (3.5)、(6.4)、(6.5)，
\[
 E_6-Q(v)\mathcal V=o(v^6S).
\]
因而有准确的归一化比较
\[
 \begin{split}
 \bar E-Q(v)\mathcal V
 =\Omega^{-1}\bigg[&e_7+
   \tau\{E_6-Q(v)\mathcal V\}\\
 &+\sum_{h=1}^5\tau^{7-h}\{E_h-Q(v)\mathcal V\}\bigg]
 =o(v^7S).                                                         \tag{7.5}
 \end{split}
\]
最后五项的总权为 \(O(v^2)\)，其成本为 \(O(v^8S)\)。第六层的一阶领先价抵消了将第七层除以 \(\Omega=1+\tau+O(v^2)\) 时的一阶变化；不需要 \(b\) 的精确值或第六层余项的收敛率。

加权 (7.2)、用 \(f_7\ge0\)，再除以 \(a\Omega\)，得到
\[
 \mathcal V-A_vM/a\ge
 c_{\rm rec}P_I-(C_{\rm b}+bv)P_\partial-4vS_I
 -\bar E-O(v^8S).
\]
取 \(bv\le1\)，按 (2.1) 至少保留 \(c_{\rm rec}P_I/2\ge0\)。由 (7.5)，
\[
 \boxed{\mathcal V\ge A_vM/a-Q(v)\mathcal V-o(v^7S).}                \tag{7.6}
\]
外部为空时取零源、不定义空 CDF，内部计数仍给同方向；实际上 (1.1) 的小 \(v\) 正质量分支不会仅由零源造成例外。

## 8. 真实分母及 \(-176/243\)

由 \(S/V\le(1-a)^{-1}\)，(7.6) 的误差是 \(o(v^7V)\)。写 \(x=M/V\)，移去实际 \(M\) 得
\[
 (1+Q(v))(1+x)\ge (A_v/a)x-o(v^7).
\]
\(A_v-a(1+Q(v))\) 在 \(v=0\) 等于一，因此正且远离零，从而
\[
 1+x\le\frac{A_v}{A_v-a(1+Q(v))}+o(av^7).                            \tag{8.1}
\]
误差在移分母时多乘一个 \(a\)，这是七阶 halo 精度所需的因子。

精确差式为
\[
 \frac{A_v}{A_v-a(1+Q)}-\frac{A_v}{A_v-a}
 =\frac{aA_vQ}{(A_v-a)(A_v-a-aQ)}.                                  \tag{8.2}
\]
由于 \(Q=O(v^6)\)、\(a=v^2\)，替换第二分母为 \(A_v-a\) 的误差为 \(O(a^2Q^2)=o(av^7)\)。只须一阶展开：
\[
 Q(v)=kv^6-\frac{26}{9}kv^7+O(v^8),\qquad
 \frac{A_v}{(A_v-a)^2}=1-2v+O(v^2).
\]
因此
\[
 1+M/V\le L_2(t)+kav^6-
       \left(2+\frac{26}{9}\right)kav^7+o(av^7)
 =L_2(t)+\frac4{27}av^6-\frac{176}{243}av^7+o(av^7).                 \tag{8.3}
\]
两次实际半径归一化、最后的真实 \(S/E\) 分母处理和反馈退款均已保留；没有将七阶系数从旧粗 \(S/V\) 界直接读出。

## 9. 上确界、全部尺度与条件范围

设 (1.3) 失败。取固定 \(\epsilon>0\) 和 \(v_j\downarrow0\)，先按 \(C_2\) 的上确界定义选有限正测度输入，损失小于 \(\epsilon av^7\)，再对外部严格 halo 用紧集内逼近，损失也小于 \(\epsilon av^7\)。可得到
\[
 1+M_j/V_j\ge L_2(t_j)+kav_j^6+
       \left(-\frac{176}{243}+\epsilon\right)av_j^7.
\]
这些输入满足 (1.1)，可取一个固定 \(B\)，与 (8.3) 矛盾。外部 halo 的有限性、可测性及输入几乎处处进入严格 halo，均用原已证性质；不把 \(K\) 当作输入补集中的开集。有限贪心给有限族，全部原半径层被保留；R019/R027 已在任意尺度跨度下处理窗口和坏根，所以这里没有额外等尺度或周期假设。

完整次序为：纯维数回收及图参数 → 由原 \(P,S\) 固定同一个 \(I\) → 七深度递归及统一 \(A_7\) → 强 CDF 斜率 → 对六层/七层前缀各做粗先验 → 已证领先价格的固定份额、阈值和顺序极限 → 由饱和推出 (4.2) → 同一第七层条件价格 (1.2) 及 R027 → 七层抵消及上确界。源参数不反馈到图选择；第六层和第七层没有更换输入、标记或原指定见证。

**剩余最小接口。** 本总装不再需要另行假设 \(\mathcal F_7=o(v^6S)\)、小源删除严格小量或第六层七阶价；它们分别由 §4、§5、§6 得到。唯一仍待 R028 独审和供货确认的是 (1.2) 对本次重建第七层的全 pair 完整源价及其任意输入序列 \(o(v^7S)\) 量词。只证旧六层对象上的 (1.2)，尚不足以直接宣布 (1.3)。

## 10. 实跑算术证据与依赖版本

已实际读取 Terra 的本机 Mathematica 15.0 运行结果，独立目录为 `r030-computation/`：`seven_layers.wl`、`seven_layers.log`。其固定符号 \(b,c_D\) 的七层望远镜残差为零；权平均减去 \(Q(v)\) 的展开到 \(v^7\) 全为零；\(f_0=O(v^3S)\) 的初端归一化为 \(O(v^8S)\)；由精确 \(\mathcal V\) 分母算出的两个系数是 \(\{4/27,-176/243\}\)，与目标比较返回 `True`。这些是有限代数证书，不负责几何或未审局部价格。

版本：

- R028 候选 `r028-feedback-upper/planar_second_price.md`：SHA256 `615fefa82b0df856106c58fd2906f00d77b99fffbc8f3a7957500037263ee618`，本稿仅条件引用其 (2)。
- 已证 R027 `proofs/027-actual-radius-incidence.md`：SHA256 `59c99a99648bdc2bc9fa60ab6160db6693005e19b018a6056304f084eca7a00c`。
- 算术脚本 SHA256 `876c0d8191c7c100a9402e652d2aefddf31a0fd3756a5df81d4402a4bdc8c37a`；日志 SHA256 `62acf028e9467ec971eb1bb5e1b6b37cbfb80238eaaa8ad10542c84be80059db`。

实际查阅的证明接口包括 R021 修正版、正式 `global.tex` 的小源回收完整证明和图选层、正式 `sharp.tex` 的保留正量递推/领先价格/公共账本、R026 有限构造结论、R027 全文及 R028 当前全文。没有改动这些依赖源。
