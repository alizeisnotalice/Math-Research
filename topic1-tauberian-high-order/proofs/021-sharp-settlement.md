# R021：末层局部锐价格到全输入 limsup 的条件性结算

日期：2026-11-05。完整条件性 assembly，待独立审查。只写独立 scratch。

本稿不证明 R20 的真实角成本。唯一新增局部假设在 §2 精确写出。其余内容是六层粗先验、同一末层接口、内部回收、望远镜及最后移 \(M\) 的严格组装。R19/R20 独审完成前，不据此宣布真实锐定理。

## 1. 归一化及必须保持的量词

固定 \(n\ge2\)，\(H=2^n\)，记
\[
 a=v^n=2^nt,\qquad A_v=(1+v)^n,\qquad
 L_n(t)=\frac{A_v}{A_v-a},\qquad
 \mathfrak k_n={n\choose2}\frac{n^2}{27}.
\]
对有限正测度输入 \(E\) 写 \(V=|E|\)。取外部严格 halo 中紧集 \(K\)，\(M=|K|\)。原有限严格见证贪心给 \(Q_i=Q(x_i,r_i)\)，\(B_i=Q(x_i,r_i/2)\)，\(S=\sum_i r_i^n\)。

只需研究任意序列 \(v\to0\)、\(E=E_v\)、\(K=K_v\) 满足
\[
 \frac MV\ge\frac a{A_v-a}-Bav^6,\qquad B<\infty\text{ 固定}. \tag{1}
\]
下文所有 \(o_\sigma(1)\) 都必须对**任意这样的输入序列及所构造同一六层族**成立；仅对固定 \(E\) 的点态渐近不够。允许该序列的模量不同，但不允许先固定输入后再对输入取上确界。

两种小参数分开记：\(\eta_{\rm w}\) 是 R14 扩大窗口的固定参数；\(0<\sigma<1/4\) 是 R20 局部价格中保留份额的参数。先固定 \(\sigma\)，再 \(v\to0\)，最后 \(\sigma\downarrow0\)。

记标记集为 \(I\)，外边界为 \(\partial I\)，
\[
 P_I=P(I),\quad P_\partial=P(\partial I),\quad S_I=S(I),\quad
 K'=K\setminus\bigcup_{i\in I}Q_i,\quad M'=|K'|,\quad M_I=M-M'.
\]
全部六层都使用同一个 \(I,K',M'\)。

## 2. 唯一局部锐输入

六层源目标按 R14 的固定有限深度构造：未删除源 \(T_{i,h}\)，有效源 \(F_{i,h}\)，目标 \(G_{i,h}\)，质量分别为 \(m_{i,h}^0,\mu_{i,h},g_{i,h}\)。小源删除量 \(\delta_{i,h}=m_{i,h}^0-\mu_{i,h}\)。记
\[
 \mathcal D_h=\sum_i g_{i,h}-M',\qquad
 \mathfrak d_h^0=M'-\sum_i m_{i,h}^0,\qquad
 \Delta_h=\sum_i\delta_{i,h}.
\]
这里 \(\mathcal D_h,\mathfrak d_h^0\ge0\)，且有精确恒等式
\[
 \sum_i(g_{i,h}-\mu_{i,h})
       =\mathcal D_h+\mathfrak d_h^0+\Delta_h. \tag{2}
\]
边界源仍是原 \(P_i\)；深处源在 \(h\ge2\) 为前一层中心邻域。所有点沿用原指定严格见证，不重选半径。

末层写
\[
 R_i=(g_{i,6}/a)^{1/n},\quad
 X_i=x_i+W_{i,6}+[-R_i/2,R_i/2]^n,
\]
\[
 J_i=(R_i+|W_{i,6}|^{1/n})^n-(R_i+\mu_{i,6}^{1/n})^n\ge0.
\]
\(J_i\) 是完整的实际源盈余，不是粗结算消费后“剩余的一小份”。

**假设 H（末层鲜盈余价格）。** 对 §3 的六层粗先验重新选出的同一 \(I\) 和末层源盒，每个固定 \(\sigma\) 有
\[
 \sum_{i<j}|X_i\cap X_j|-\sum_iJ_i
 \le\left[\mathfrak k_{n,\sigma}+o_\sigma(1)\right]v^6S,
 \qquad \mathfrak k_{n,\sigma}=\frac{\mathfrak k_n}{(1-\sigma)^2}. \tag{H}
\]
和包含所有外部 pair，含边界、多临界及异常胞关联对；单临界对无正体积交。

R20 提供的更原始形式是右侧
\[
 \frac{n^2}{27(1-\sigma)^2}v^6\mathcal C_6+o_\sigma(v^6S),
 \quad \mathcal C_6=\sum_{e:\,q(e)=2}b_er_e^2.
\]
若将 R19 的输入桥及 R17 reviewed 的 incidence 应用于下面重新构造的第六层，得到
\(\limsup\mathcal C_6/S\le\binom n2\)，便是 (H)。本稿不重证该局部成本，也不把尚未独审的 R20 结论作为无条件定理。

## 3. 六层先验：不能复用旧第五层的 I 和 F

**引理1。** R17 的粗临界紧性可对六个固定深度重新构造。由此选出一个 \(I\)，使原图条件
\[
 C_*P_\partial\le c_{\rm rec}P_I/4,
 \qquad 2nvS_I\le c_{\rm rec}P_I/4 \tag{3}
\]
成立，且同一第六层有
\[
 P_I+\mathcal D_6/a+v\sum_i(\varepsilon_{i,6}+q_{i,6})r_i^n
      +\mathcal F_6+c_{n,\tau}\mathcal F_{\tau,6}=O_{n,B}(v^6S), \tag{4}
\]
\[
 \mathcal D_5/a=O(v^5S),\quad
 \mathfrak d_6^0/a=O(v^5S),\quad \Delta_6/a=O(v^6S),\quad
 S(I)+S(\partial I)=O(v^5S). \tag{5}
\]
对每个固定稳健面裕量 \(\tau>0\) 依原口径成立。

证明：首先把 R14 的有限深度递推完整执行到 \(h=6\)：每一步先选足够大的深目标系数及中心邻域常数，再选下一步常数。只多执行一次固定步骤；统一为 \(A_6\) 后再选强 CDF 斜率 \(\lambda_0>0\)，最后取小 \(v\)。所有完整大见证、源目标不交及边界改为 \(P_i\) 的论证逐层不变。具体取 \(N_{i,h}=K'\cap Q(x_i,L_hvr_i)\)，这些中心邻域在每层两两不交，令 \(f_h=M'-\sum_i|N_{i,h}|\ge0\)。原深目标的唯一点定位给 \(\mathcal D_h\ge f_h\)，边界仍改回私有源的质量差由 \(aP_\partial\) 支付；于是同一六层满足
\[
 \mathcal D_h\ge f_h\ge0,\qquad
 \mathfrak d_1^0\le A_0vM',\qquad
 \mathfrak d_h^0\le f_{h-1}+aP_\partial\quad(2\le h\le6). \tag{6}
\]
这里未删除源的缺量不含 \(v^nM'\) 删除地板。

先沿用 R14 预选纯维数常数 \(c_{\rm rec}=1/(2H)\)、\(C_{\rm b}=2C_2(n)+2\)、\(C_*=C_{\rm b}+2\)。令 \(Z_0=a[V+M-A_vM/a-c_{\rm rec}P_I+C_{\rm b}P_\partial+2nvS_I]\)。对全部六层重新使用 R17 的保留正量粗结算。与其五层证明相同，重复量留一半，另一半用于 \(\mathcal D_h\ge f_h\)；记源正预算为 \(\mathcal J_h\ge0\)，其包括 \(v\varepsilon\)、\(vq\)、原免费面缝及新增稳健面缝。递推为
\[
 Z_0\ge c_0 f_h-b_0v f_{h-1}+a\mathcal J_h+\mathcal D_h/4
                 -C_na(v^6S+vP_\partial),\qquad f_0=A_0vM'. \tag{7}
\]
按 \(\rho_0^{6-h}\) 加权，\(\rho_0=b_0v/c_0\)，中间五个 \(f_h\) 精确消去。初端为 \(O(v^7M')\)，除以 \(a\) 后至多 \(O(v^7S)\)，而粗余项仍为 \(O(v^6S)\)。权重和 \(1+O(v)\)，保留 \(\mathcal J_6,\mathcal D_6/a\) 以及 \(v\mathcal D_5/a\) 的固定正比例。

图选层只用原 \(P,S\)，所以仍可按原纯维数参数选出满足 (3) 的 \(I\)。由 (1) 得同一精确 \(V+M-A_vM/a=O(v^6S)\) 上界，便得到 (4) 和 \(\mathcal D_5/a=O(v^5S)\)。再由 (6)、R14 小源回收和 (2)，如 R17 原证明逐式得到 (5)。这里 \(\mathfrak d_6^0\) 只能保证为 \(O(av^5S)\)，没有声称是 \(o(av^5S)\)。

这份先验针对新六层 \(I,F_{i,6}\)；旧第五层 \(I,F_{i,5}\) 并未拿来替代。R19 的 LT 提取只用 \(S(I)=O(v^5S)\)、外部规则半盒和完整图分离；其角均匀可积性只用末层 \(\sum\varepsilon_i r_i^n=O(v^5S)\)、\(\mathcal F=O(v^6S)\)。这些条件现对新 \(I,F_{i,6}\) 成立。因此其证明和 incidence 原封应用于本第六层，无需任何一致网格模型。证毕。

粗先验不物理改变 \(W_{i,6},F_{i,6}\) 或 \(J_i\)。它是关于这些实际集合的一个预估计；最终锐结算重新使用完整定义，不能把第一次证明中的预算分份误当成集合中被删除的质量。

## 4. 从 H 到未折损的末层实际并集

末层取 \(U_6=\bigcup_{i\notin I}X_i\)，不加为粗先验提供正预算的面条。对原共同环境缺陷记
\[
 u_{i,6}=|(X_i\setminus E)\setminus G_{i,6}|.
\]
R14 的共同可见性及环境控制对这个较小并集仍成立，故
\[
 |U_6\setminus E|\le M'+\mathcal D_6+\sum_i u_{i,6},\qquad
 \sum_i u_{i,6}\le C_na(\mathcal F_6+P_\partial)=O(av^6S). \tag{8}
\]
最后使用 (3)、(4)。环境项相对 \(v^6S\) 是 \(o(1)\)，而不是把 \(P_\partial\) 的固定价格丢进误差。

二项 Bonferroni、完整 \(J_i\) 及 (H) 给
\[
 |U_6|\ge\sum_i(R_i+\mu_{i,6}^{1/n})^n
                 -[\mathfrak k_{n,\sigma}+o_\sigma(1)]v^6S.
\]
标量基线不等式的系数仅需一个固定粗上界：
\[
 (R_i+\mu_{i,6}^{1/n})^n
 \ge A_vg_{i,6}/a-C_{n,A_6}(v/a)(g_{i,6}-\mu_{i,6}). \tag{9}
\]
把 (2) 代入，并从左侧扣除物理 \(\mathcal D_6+\sum_i u_{i,6}\)。\(\mathcal D_6\) 的净系数为
\[
 (A_v-C_{n,A_6}v-a)/a\ge1/(2a).
\]
小源成本由 (5) 为 \(C(v/a)\Delta_6=O(v^7S)\)，环境项 (8) 为 \(O(av^6S)=o(v^6S)\)。因此有
\[
 \boxed{\begin{aligned}
 |U_6|-\mathcal D_6-\sum_i u_{i,6}
 \ge{}&A_vM'/a+\mathcal D_6/(2a)
       -b(v/a)\mathfrak d_6^0\\
 &-[\mathfrak k_{n,\sigma}+o_\sigma(1)]v^6S,
 \end{aligned}} \tag{10}
\]
其中 \(b>0\) 固定，可增大以统一早五层常数。(10) 与 R20 的条件式(13)完全对齐。

必须保留 \(-b(v/a)\mathfrak d_6^0\)：由 (5) 它可能正好是 \(O(v^6S)\)，不能吸收进 \(o(v^6S)\)。下一节和前层的正重复量精确抵消它。

## 5. 同一 V+M 的内部与边界结算

对早五层取 R14 粗增广并集 \(U_h^+\)，对末层取上面的 \(U_6\)。它们可以是不同实际并集，但所有深度的左侧最终都是同一 \(V+M\)。末层去掉面条不会增加内部占用；边界源仍为原 \(P_i\)，故所有层的纯维数突出价格不变。

记内部突出集 \(Z_I=(\bigcup_{i\in I}A_i^\beta)\setminus\bigcup_jB_j\)，\(\beta=1-\eta_{\rm w}\)。对每层加入 \(B_I\) 和 \(Z_I\) 后，实际同一 \(W_h\) 满足
\[
 |W_h|\ge |U_h|+S_I+H^{-1}\sum_{i\in I}e_i^\beta r_i^n
                              -2C_2(n)P_\partial, \tag{11}
\]
\[
 |W_h|\le V+M'+\mathcal D_h+\sum_i u_{i,h}+aS_I. \tag{12}
\]
\(V\) 只计一次，补集份额允许重复上界。完整 source/target 的重叠已经在 \(\mathcal D_h\) 处理。

内部质量计数给
\[
 S_I-M_I/a\ge\tfrac12\sum_{i\in I}s_ir_i^n-P_\partial/2,
 \qquad M_I/a\le S_I.
\]
相应内部余量的精确恒等式为
\[
 (1-a)S_I+M_I-A_vM_I/a
 =(1-a)(S_I-M_I/a)-(A_v-1)M_I/a. \tag{13}
\]
取纯维数
\[
 c_{\rm rec}=1/(2H),\quad C_{\rm b}=2C_2(n)+2,
 \quad C_*=C_{\rm b}+2.
\]
用 \(a<1/2\)、\(A_v-1\le2nv\)，(11)--(13) 与原单层并集结算给同一公共左侧
\[
 \begin{aligned}
 V+M\ge{}&A_vM/a+c_{\rm rec}P_I-C_{\rm b}P_\partial-2nvS_I\\
 &+c_D\mathcal D_h/a-b(v/a)\mathfrak d_h^0-E_h,
 \qquad c_D=1/2,\\
 E_h={}&C_0v^6S\quad(1\le h\le5),\\
 E_6={}&[\mathfrak k_{n,\sigma}+o_\sigma(1)]v^6S.
 \end{aligned} \tag{14}
\]
可用非负上界取代最后的 \(o\)，使 \(E_h\ge0\)。早层的 \(C_0,b\) 在六层源常数固定后统一选择。\(C_{\rm b}\) 中的纯维数预留单位吸收单层 \(C_{n,A_6}vP_\partial\)，只需随后取小 \(v\)，不倒过来令纯维数价格依赖 \(A_6\)。

若外部为空，全部源量取零，不定义空 CDF；内部计数自身给相同方向的不等式，后续仍成立。

## 6. 六层精确消去：前五层只花 O(v) 权重

置
\[
 Z=a[V+M-A_vM/a-c_{\rm rec}P_I+C_{\rm b}P_\partial+2nvS_I],
 \quad f_0=A_0vM',\quad \tau=bv/c_D.
\]
由 (6)、(14)，每层都有
\[
 Z\ge c_Df_h-bvf_{h-1}-aE_h-abvP_\partial,
 \qquad1\le h\le6. \tag{15}
\]
对第一层不需要边界附项，但保守加入同一非负成本不影响方向。

取 \(w_h=\tau^{6-h}\)、\(W=\sum_{h=1}^6w_h=1+\tau+\cdots+\tau^5\)。精确恒等式为
\[
 \sum_{h=1}^6w_h(c_Df_h-bvf_{h-1})
 =c_Df_6-bv\tau^5f_0. \tag{16}
\]
特别，末层的 \(-bvf_5\) 与第五层的 \(w_5c_Df_5=\tau c_Df_5=bvf_5\) 完全消去。这里没有假设 \(f_5=o(av^5S)\)，也没有丢掉 \(d_6^0\) 的领先成本。

将六个 (15) 加权，除以 \(aW\)，利用 \(f_6\ge0\)：
\[
 \begin{aligned}
 V+M-A_vM/a\ge{}&c_{\rm rec}P_I-(C_{\rm b}+bv)P_\partial-2nvS_I\\
 &-\frac{E_6}{W}-\frac{C_0v^6S(W-1)}{W}
  -\frac{bv\tau^5f_0}{aW}.
 \end{aligned} \tag{17}
\]
前五层总权 \(W-1=O(v)\)，所以其粗误差只有 \(O(v^7S)\)。末层价格只除以 \(W=1+O(v)\)，领先系数保持 \(\mathfrak k_{n,\sigma}\)。初端精确为
\[
 \frac{bv\tau^5f_0}{aW}
 =\frac{A_0b^6}{c_D^5W}v^7\frac{M'}a=O(v^7S), \tag{18}
\]
因为 \(M'\le aS\)。这就是从五层加到六层的实质作用。

取小 \(v\) 使 \(bv\le1\)，边界增量放进 \(C_*-C_{\rm b}=2\) 的预留。代入同一 \(I\) 已满足的 (3)，至少保留 \(c_{\rm rec}P_I/2\ge0\)。因此
\[
 \boxed{V+M\ge A_vM/a-[\mathfrak k_{n,\sigma}+o_\sigma(1)]v^6S.} \tag{19}
\]
到此为止没有将 \(aS\) 近似成 \(M\)，也没有物理删除异常胞的基线质量。

## 7. 最后移 M 的准确分母

原半盒内部不交且 \(B_i\subset Q_i\)，严格见证给
\[
 S=|\bigcup_iB_i|\le V+\sum_i|B_i\setminus E|
 \le V+\sum_i|Q_i\setminus E|\le V+aS.
\]
故对所有原输入精确有
\[
 \frac SV\le\frac1{1-a}. \tag{20}
\]
在 (19) 中直接移去实际 \(M\)。写 \(k_v=\mathfrak k_{n,\sigma}+o_\sigma(1)\ge0\)，得到
\[
 \frac MV\le\frac a{A_v-a}\left[1+k_vv^6\frac SV\right].
\]
于是
\[
 \boxed{\frac{1+M/V-L_n(t)}{av^6}
 \le\frac{k_v}{(A_v-a)(1-a)}.} \tag{21}
\]
右边的额外分母趋于一，既不需要 \(S/V=1+o(v^6)\)，也不需要 \(aS=M+o(av^6S)\)。每个固定 \(\sigma\) 先令 \(v\to0\)，再令 \(\sigma\downarrow0\)，给
\[
 \limsup_{v\to0}\frac{1+M/V-L_n(t)}{av^6}\le\mathfrak k_n. \tag{22}
\]

## 8. 条件性全部输入结论及常数顺序

**条件定理。** 若假设 H 对任意临界序列 (1)、其按 §3 重新构造的同一六层末层成立，则
\[
 \limsup_{t\downarrow0}
 \frac{C_n(1-t)-L_n(t)}{av^6}\le\mathfrak k_n.
\]

证明：若不成立，取某个 \(\delta>0\) 及 \(v_j\to0\)，使上确界比值超过 \(L_n+(\mathfrak k_n+3\delta)av^6\)。按上确界定义选择有限输入，损失小于 \(\delta av^6\)；再用外部 halo 的紧集内逼近，损失也小于 \(\delta av^6\)。密度定理给输入几乎处处进入 halo，外部 halo 可测且有限，所以得到
\[
 1+M_j/V_j\ge L_n+(\mathfrak k_n+\delta)av^6.
\]
这些输入满足 (1)，可取固定 \(B=0\)，与 (22) 矛盾。紧集不需是 \(E^c\) 中的开集；它只需包含在外部严格 halo 中，原半径严格见证和内逼近均合法。证毕。

参数次序完整如下：先固定 \(n,c_{\rm rec},C_{\rm b},C_*\) 及原全图增长常数；再选原有限图层数 \(K_0\)、窗口 \(\eta_{\rm w}\)、大固定规则阈值 \(D\)；随后一次完成六个深度、统一 \(A_6\)、正强 CDF 斜率 \(\lambda_0\)。按这些固定参数建立六层粗先验和其 \(I\)。对锐价格另固定 \(\sigma\)，再选 R20 的有限截断及辅助阈值，最后统一取小 \(v\)。先 \(v\to0\)，再 \(\sigma\downarrow0\)。\(B\) 只影响先验界的常数，不改变纯维数内部回收价格；没有把之后的 \(A_6\) 或局部误差常数反馈到 \(C_*\)。

当几何／测度输入外部为空时，以零源口径处理；不在临界分支的序列无需本定理。全局 limsup 的反证只选择高于基线的序列，已足够涵盖 Tauberian 上确界，不遗漏一般输入。

## 9. 本机符号核验及依赖范围

实际读取 Terra 在 `r021-computation/six_layer.wl`、`six_layer.log` 的 WolframScript 15.0 独立实跑结果：六层望远镜恒等返回 True；权重和及粗权重除以 \(v\) 有精确展开；端点 (18) 恒等返回 True。最终比值用两个非负 slack 验证准确恒等式，差为零，正分母与比值不等式的违反集合均为 False，\([(A_v-a)(1-a)]^{-1}\) 的极限为一。

这些核验只支持本稿的符号消去与分母。假设 H 的真实角成本归 R20 独立审查；输入到 incidence 的几何归 R19 独立审查。本稿不以软件输出替代这些条件，不将条件 assembly 写成已经确认的真实锐上界。
