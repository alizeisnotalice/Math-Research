# R062：同一实际源族中可增长测试窗的定量补集候选

2026-11-05。固定整数 \(n\ge2\)。本文是待独立审核的局部供货候选；不宣布新的全局 Tauberian 余项定理。仓库及既有冻结供货只读。

目标是把 R055／正式 `eighth:exact-complement` 的固定精度结算改成一个带显式常数增长的有限参数不等式。只改变从已有二轴边中选取的测试 packet 及其软权；原 \(I\)、八层递归、原见证半径、源／目标盒、正常集合、旧正常层、正常 clean 分配和强 CDF 均保持不变。本文不由“对每个固定 \(T,m\)”直接代入 \(T(v),m(v)\)。

## 0. 供货范围与冻结依赖

本稿把正常模型净价 \(\mathsf P_N\) 作为实际数值保留。主式先对全部二轴补集给出明确 \(O(v^{65/8}S)\)。二维高支撑为空；维数至少三时，接入同一族上已独审的 R061 高支撑纯 gross 界 \(O(v^{43/5}S)\)。因 \(43/5>65/8\)，最终得到全部外部交对的同阶定量补集接口。本文只读取 R061 的准确输入、输出及正式审查范围，不重复其几何证明。

实际读取、复核的主依赖如下。完整路径均相对 `/Users/zhengzhihao/Desktop/T/topic1-tauberian-work/`，另有正式论文路径。

| 供货 | 文件 | SHA-256 |
|---|---|---|
| 原第八层稳定性 | `scratch/r049-second-variation-upper/competitive_eighth_stability.md` | `5c4d7b371a0c33c39449afd686de599c49570aaa639d4f7db7d3c1f166a040f4` |
| 原全部半径轴帽 | `scratch/r051-refined-axial-cap/refined_axial_cap.md` | `37bbb0ec63fac1dce33fb9883e14499488baef490341df32a6fdfb0da427c979` |
| 原反单调算子及同份强 CDF | `scratch/r052-double-feedback-upper/actual_cap_antitone_iteration.md` | `d655c1e160c88914774248c81688c6bdde3f27614387c4eb7f0b9c1a0a0a5296` |
| 原整源补集 | `scratch/r055-eighth-complement/whole_height_complement.md` | `0684d01f3665a5a064b4786b6eb1d78ea5780ccdae64f96df3ad01647224a1a1` |
| 二维原位桥／同族接口 | `scratch/r057-eighth-global-assembly/global_eighth_assembly.md` | `e3ed7614d8775cd57fe6792ab6235bfdb8ec0029ec94532564052ba67d78277c` |
| 实际 adaptive clean 分配 | `r036-weighted-slots/adaptive_soft_slots_r036.md` | `c6b761427b26c1ab88d2ecf17f47ef856cb7c9565d53dc980bcb78fbd6e75c3d` |
| 同族高支撑纯 gross | `scratch/r061-high-support-rate/high_support_rate.md` | `a037675512160f4717f49f37f466f0579fb386e8cdc67078fa1a4a34b8b9c663` |
| R061 正式独审 | `scratch/r061-review/audit.md` | `f2a74b5ba9bb0d07c508ef42d5ae568a451016f74a195b3917c9be99aea6cb02` |
| 当前正式合本 | `/Users/zhengzhihao/Desktop/T/Math-Research-topic1-tauberian/topic1-tauberian-high-order/paper/sections/eighth.tex` | `175b314a9f076c39583a8926921cdd8f99921b3579f8c2c9860df4221006f151` |

原全域标量强制性、固定正常小管长条以及边界分解，使用上述正式合本 `eighth:global-coercivity`、`eighth:long-cap-fraction` 和 `eighth:exact-complement` 中已签原实际证明。本文追踪其参数常数，不变更原几何定义。

## 1. 固定输入及常数口径

固定竞争常数 \(B<\infty\)。沿任意真实有限正测度输入及严格 halo 外部紧集的序列，假设 \(1+M/V\ge P_{7,n}(t)-Bav^8\)，其中 \(P_{7,n}\) 是已证七阶表达式。取该序列的同一 R049 八阶竞争族，记
\[
a=v^n,\quad S=\sum_i r_i^n,\quad s_v=(1+v)^{n-1},\quad
\lambda_v=s_vv/a,\quad A_6(v)=\frac{n^2}{27(1+v)^2}v^6.
\]
全部指标都是实际已有源。令
\[
\mathscr D_i=B_i^0\setminus F_i,\qquad
\eta_i=|\mathscr D_i\setminus Y_{i,N}^{\rm old}|,\qquad
\mathcal R_*=(v/a)\sum_i\eta_i.
\]
原正常集合 \(\mathcal N\)、正常模型迭代 \(Y_{2,N}\)、槽容量 \(\mathcal C\)、槽欠量 \(D_N\) 及全部正常占用均先固定。使用
\[
\mathcal R_*+\mathcal R_{\rm all}^{\rm old}+\mathcal F+
 \mathcal F_{1/128}+\mathfrak E\le C_*v^8S,
\qquad D_N\le C_*v^2S,\qquad S(\partial I)\le C_*v^5S,
\tag{1}
\]
其中 \(\mathcal R_{\rm all}^{\rm old}\) 及其原固定大管都是旧稳定构造的对象，不是随后选择的 \(B_A^-\)。添加这两个旧正量是为了准确接入 R061；本稿二轴证明只用 \(\mathcal R_*,\mathcal F,\mathfrak E\)。
\[
X_\xi:=\sum_i\xi_iR_i^n\le C_*v^{5/2}S,
\qquad \sum_i\varepsilon_iR_i^n\le C_*v^5S,
\qquad W_N=\mathcal C-D_N+O(v^3S).
\tag{2}
\]
实际盒各边属于 \([c_wvR_i,C_wvR_i]\)。所有原邻接度和逆度有固定上界，相交根比为 \(1+O(v)\)，真正面邻居及其切向整宽裕量均是原构造的固定对象。每个原源点只用它已有的严格见证
\[
G_i\subset Q(x,r(x)),\quad R_i<r(x)\le(1+Av)r_i,\quad R_i/r_i=1+O(v).
\tag{3}
\]
各 \(F_j\) 与自身 \(G_i\) 的排他关系不变。

以下 \(C,c\) 只依固定维数、原构造参数、固定正常窗和 \(C_*\)，不依输入数目、尺度跨度、\(v,T,m\)。可以在证明各处增大 \(C\)，再一次取足够小的 \(c_0\)。局部匹配只要求原好根完整清单达到一个固定精度；这个精度与 \(T,m\) 及软权常数无关，详见第 3 节。

记正常实际模型净价
\[
\mathsf P_N=\sum_{e\in\mathcal N}I_e-\lambda_v\sum_i|Y_{2,N,i}|.
\tag{4}
\]
它没有扣全体 \(J_i\)。任取原先固定的
\[
q_N(v)=A_6(v)+O(v^7),\qquad
\mathsf E_N=\mathsf P_N-q_N(v)W_N.
\tag{5}
\]
本文不把 \(\mathsf E_N=o(v^8S)\) 擅自改成某个幂次，也不调用正常价格来再次消费 \(J\)。

## 2. 统一参数域和候选结论

正常两交长位于固定区间 \([c_Nv^3R_e,T_Nv^3R_e]\)。取
\[
T\ge T_0\ge\max(T_N,1),\quad
0<m\le m_0<\min(c_N/2,1),
\tag{6}
\]
并要求
\[
Tv\le c_0,\qquad Tm^{-3}v\le c_0,\qquad T^5v^2\le c_0.
\tag{7}
\]
这不是固定 \(T,m\) 的渐近约定，而是每个实际有限输入都能检查的共同参数域。特别它保证 \(Tv^3\ll v\)、\(T^2v^4\ll v\)、\(T^2v^3\ll1\) 以及长条中的 \(T\le c_0/v\)。

对两端深处、两轴正交长的已有边定义
\[
\mathcal B=\{e\notin\mathcal N:m\le p_e,q_e\le T\},\qquad
\mathscr P=\mathcal N\cup\mathcal B,
\]
其中 \(p_e=L_1/(v^3R_e),q_e=L_2/(v^3R_e),R_e=\max(R_i,R_j)\)。临界轴判据和边集合仍是原来的；只选择其中一个子集。

令 \(\mathsf M\) 为该 packet 的完整非正常名义高度成本，\(\mathcal R_A^-\) 为共同真实保证层之外的 raw 成本，精确定义见第 5 节。取
\[
H_2=\sum_{q(e)=2}I_e-\sum_iJ_i,
\]
包括二轴边界边和任意非好端点。候选结论为存在固定正 \(c_D,c_B,c_R,c_E\) 使
\[
\boxed{
H_2+c_Dv^6D_N+c_B\mathsf M+c_R\mathcal R_A^-+c_E\mathfrak E
\le q_N\mathcal C+\mathsf E_N+C\,\mathcal E(v,T,m)S,
}
\tag{8}
\]
其中
\[
\boxed{
\mathcal E(v,T,m)=
\frac{v^8}{T}+Tm\,v^8+Tm^{-3}v^9+T^5v^{10}
 +T^2v^{11}+T^3v^{17/2}+m^{-1}v^9.
}
\tag{9}
\]
全部外部交对版本的右端先加 \(\sum_{q(e)\ge3}I_e\)，再用 R061 的 \(O(v^{43/5}S)\) 纯 gross 界；第 10 节给最终式。零／单临界正体积交已由原真实面容量引理排除。此处没有把高支撑的固定精度小量冒充幂次。

## 3. 同一 clean 几何上的软常数线性增长

需要的加强不是引用“每个固定 \(C_s\)”的结论，而是重放原有限证明，得到对所有有限 \(C_s\ge1\) 同时成立的
\[
W_{C_s}:=\sum_e b_eR_e^2\min(1,C_s\min(p_e,q_e))
\le\mathcal C+C(1+C_s)v^3S.
\tag{10}
\]

在 \(n\ge3\) 时，固定原 probe 大小 \(\epsilon\)、窗口清单 \(\mathcal I_i(20)\)，及匹配精度 \(\epsilon/100\)。每边的 clean 阈值仍取
\[
\ell_e=\min L_l,\qquad \lambda_e=\ell_e/(1000r_i)>0.
\]
所有临界交长原来就满足 \(\ell_e\le Cv r_i\)，所以 \(\lambda_e<\epsilon^2\) 的起效阈值不依 \(C_s,T,m\)。probe 捕获、实际 helper 的唯一性、跨锚清单识别及不同 \(\lambda_e\) 的两对角冲突都没有使用软权数值。尤其共享一个实际角面点的两条不同边，由原精确端点恒等式迫使至少一个相邻面侵入量大于第一边的 \(\ell_e/2\)，触发它自己的 \(\lambda_e\) 碰撞。无需两阈值相同，也无需阈值有正下界。

收费中唯一含 \(C_s\) 的步骤是
\[
\frac{\omega_e}{\lambda_e}
\le\frac{C_s\ell_e/(v^3R_e)}{\ell_e/(1000r_i)}
\le C C_s/v^3.
\]
原 pair gross 总和是 \(O(v^6S)\)，固定窗口反向重数不变，故碰撞费至多 \(C C_s v^3S\)。probe 空隙费 \(O((v^{n+1}+v^5)S)\)、坏端点费 \(O(v^{9/2}S)\)、根权传输费 \(O(v^6S)\) 和宽度费全都只用 \(\omega_e\le1\)，其常数不含 \(C_s\)。同一 helper 先对片求并再收费，未引入随阈值数目增长的重数。角面单射遂给 (10)。

二维用 R057 已重放的实际四角／四胞离散计数及原根桥。其无权角和已有统一界，乘任意 \(0\le\omega_e\le1\) 只减小它；(10) 因而同样成立，且不需高维 transverse helper。

取 \(C_s=C_0/m\)，其中固定 \(C_0\) 足够大。正常边及 \(\mathcal B\) 的权均为一，正常分配和 (2) 不变。对任意其余边子集 \(\mathcal D\)，有
\[
W_B,W_{\mathcal D}\le D_N+C m^{-1}v^3S.
\tag{11}
\]
这是沿原分片仅改变权值，不是新选正常集合或清单。

## 4. 实际全部半径必要性的显式 \(T\) 账本

使用原 effective-box global-min-face chart：最近面符号先固定，全局最浅轴为 \(q\)，\(y_q=0\)，\(s\) 为该轴深度，真实纤维上界
\[
H_i(y)=\min_h(w'_{ih}/2-y_h).
\]
对目标两临界轴 \(l,k\)，保留原截幅
\[
c_{ij}=\operatorname{clip}((R_i-R_j)/2;-\max L_h,\max L_h),
\quad P_{ij,h}=(L_h+c_{ij})_+.
\]
上截幅只缩小真实可保证帽，下截幅给零帽。由 (6)–(7) 及相交根可比性，所有帽宽至多 \(CTv^3R_i\)，严格小于源半宽和目标原始全宽。仍在定义中保留真实 target-width clip，不先删掉它再论证。

同源实际目标先相加，定义 \(h_A^0\) 和实际高度 \(h_A(\mathscr D)\)。有限度及原整横宽给
\[
0\le h_A\le h_A^0\le CT^2v^4R_i,
\qquad D_{i,T}\le CTv^{n+1}R_i^{n-1},
\tag{12}
\]
其中 \(D_{i,T}\) 是总角帽边带质量斜率。在 \(Tv\le c_0\) 下它小于容量斜率 \(naR_i^{n-1}\)。因此原 \(\ell=r(x)-R_i\le s\) 分支的斜率比较一致成立，不需另选见证半径。

当 \(\ell\ge s\)，原真正浅面邻居的整横截面被同一见证捕获。用固定旧正常层定义
\[
\beta_h=(1-\theta_N(v))\prod_{l\ne q}w^0_{hl}
\asymp v^{n-1}R_h^{n-1},\qquad \theta_N(v)=O(v^2).
\]
实际轴帽下界仍是 \([\beta_h d-\eta_h]_+\)，对任意小 \(d\) 一致。全目标进入时 \(|F_h|\ge caR_i^n\) 与原剩余容量 \(O(avR_i^n)\) 矛盾，故排除饱和分支。转折点
\[
z=s+g_{i\to h}+\Delta_i+\eta_h/\beta_h
\]
满足 \(z\le CvR_i\)：\(s\le CT^2v^4R_i\le CvR_i\)，其余三项各自原本不超过 \(CvR_i\)。因而比较容量 \(C_i(z)\) 与 \(C_i(s)\) 时的导数常数不含 \(T\)。

容量曲率则有明确增长
\[
0\le C_i(s)/(naR_i^{n-1})-s\le CT^4v^8R_i.
\]
原 R051 两分支逐式给同一真实必要性
\[
s\ge h_A(\mathscr D;y)-e_i(T),\qquad
e_i(T)=C\left[T^4v^8R_i+\Delta_i+
 \max_h(g_{i\to h}+\eta_h/\beta_h)\right].
\tag{13}
\]
关键是后面三项无需再乘任意未知 \(C_T\)。有限面逆度和 (1) 给
\[
\sum_iR_i^{n-1}e_i(T)\le C(1+T^4)v^8S.
\tag{14}
\]

正帽 base 的测度，在最浅轴属于目标 pair 时至多
\(CTv^{n+1}R_i^{n-1}\)，否则至多 \(CT^2v^{n+3}R_i^{n-1}\)。由 \(Tv^2\le c_0\) 可合成前者。于是对所有实际 chart 同时侵蚀 \(e_i(T)\) 的物理损失及鲜成本满足
\[
e_T^{\rm phys}\le CT^5av^9S,\qquad
\lambda_ve_T^{\rm phys}\le CT^5v^{10}S.
\tag{15}
\]
第二类 chart 的名义主量仍保留；这里仅估计侵蚀支撑。

同样打开 Tonelli，固定一个目标物理点，可覆盖它的源 base 测度除以 \(naR_i^{n-1}\)，分别为 \(CTv\) 和 \(CT^2v^3\)。全 packet 算子的物理 \(L^1\) 常数至多
\[
\Lambda_T\le C(Tv+T^2v^3)\le CTv.
\tag{16}
\]
名义总物理量至多 \(CT^3av^5S\)；pair 外最浅轴的附加项是 \(CT^4av^7S\)，由 \(Tv^2\le c_0\) 吸入。乘鲜导数后，退款核为 \(CTv^2/a\)。

## 5. 全高度、真实 clip 与新窗成本的统一八阶界

写原同一 chart 上的非负高度
\[
h_A=h_N+h_B,\quad h_B=h_B^0-\mathcal K_B(\mathscr D)\ge0,
\qquad \mathsf M=\lambda_v\sum_{\rm charts}\int h_B^0.
\]
名义积分保留目标原始整横宽、源完整有效横宽、实际根分母，以及最浅轴不在目标 pair 的所有 chart。令
\[
B_A=\{0\le s<\min(H_i,h_N+h_B)\},\quad
B_N=\{0\le s<\min(H_i,h_N)\}.
\]
两者用同一 \(e_i(T)\) 侵蚀得到 \(B_N^-\subset B_A^-\subset D'_i\)。原旧正常层精确包含于 \(B_N\)，所以
\[
\mathcal R_A^-\le\mathcal R_N^-\le\mathcal R_*+CT^5v^{10}S,
\qquad
\mathcal R_N^-\le\mathcal R_A^-+s_v^{-1}\mathsf M.
\tag{17}
\]
末式常数是一，来自同一正部和上截断的单调 1-Lipschitz 性，不依 \(T,m\)。

真实增量的精确等式为
\[
\Delta_A:=\lambda_v|B_A\setminus B_N|_1
=\mathsf M-\mathfrak F_B-\mathfrak L_B.
\tag{18}
\]
在每个名义目标项为正之处，临界底变量只占 \(CTv^3R_i\) 窗。全高度至多 \(CT^2v^4R_i\)，而最浅轴和两个临界轴剩余半宽均为至少 \(cvR_i\)。故上截断只能发生于非临界自由底变量的中分面 \(CT^2v^4R_i\) 带。该名义项独立于这一自由变量，其全长至少 \(cvR_i\)。有限变量并给
\[
0\le\mathfrak L_B\le CT^2v^3\mathsf M.
\tag{19}
\]
二维没有该自由变量，\(\mathfrak L_B=0\)。这不是删除三近面主量。

原必要性、(15) 和旧层锚定给
\[
\Delta_A\le s_v\mathcal R_*+CT^5v^{10}S.
\tag{20}
\]
坏退款先把目标缺集分成旧正常层及其原 raw。核 (16) 给余量退款 \(CTv\mathcal R_*\)。旧正常层逐源体积是 \(Cav^5R_i^n\)，故逐坏边退款至多 \(CTv^7R_e^n\)。另一方面，两个真实有向端的名义成本由原宽下界和正部凸性至少为
\[
cv^6R_e^n p_eq_e(p_e+q_e)\ge cm^3v^6R_e^n.
\]
即使实际根差发生截幅，两个方向同时保留的凸性下界仍成立。因此明确有
\[
\mathfrak F_B\le CTm^{-3}v\mathsf M+CTv\mathcal R_*.
\tag{21}
\]
不将此式用于近零边。

合并 (18)–(21)，
\[
[1-CTm^{-3}v-CT^2v^3]\mathsf M
\le(s_v+CTv)\mathcal R_*+CT^5v^{10}S.
\tag{22}
\]
一次把 (7) 的 \(c_0\) 取得足够小，左系数至少 \(1/2\)，且 \(T^5v^{10}\le c_0v^8\)。得到真正统一、非 fixed-window 上极限式
\[
\mathsf M\le Cv^8S,
\quad \mathfrak F_B\le C(Tm^{-3}v^9+Tv^9)S,
\quad \mathfrak L_B\le CT^2v^{11}S.
\tag{23}
\]
这里是所需的 \(C_{T,m}\le CTm^{-3}\) 明确增长和吸收过程。

## 6. 只对联合真实层用一次完整 CDF

仅对固定正常算子使用原第二次迭代估计：
\[
\lambda_v|B_N(\mathscr D)|_1
=\lambda_v|Y_{2,N}|_1+O(v^9S).
\tag{24}
\]
其算子窗始终固定，所以此处无新 \(T\) 次幂。没有对增长的全 packet 再做两次迭代。

联合增量由 (18) 和 (23) 给出。对真实 \(B_A^-\) 使用原强 CDF 一次，保留 raw 和端移的固定份额。鲜导数 epsilon 损失至多
\[
C\lambda_v\sum_i\varepsilon_i|B_A(0)_i|
\le CT^3v^6\sum_i\varepsilon_iR_i^n
\le CT^3v^{11}S.
\]
于是
\[
\begin{aligned}
\sum_iJ_i\ge{}&\lambda_v|Y_{2,N}|_1+\mathsf M
 +c_R\mathcal R_A^-+c_E\mathfrak E\\
&-C[Tm^{-3}v^9+Tv^9+T^2v^{11}
 +T^5v^{10}+T^3v^{11}+v^9]S.
\end{aligned}
\tag{25}
\]
此处 \(c_R,c_E\) 来自原同一强 CDF，不随测试窗缩小。旧 \(\mathcal R_*,\mathcal F,\mathfrak E\) 的稳定性只作为已证数值上界引用，并未再抽取一份 \(J\) 来支付尾部。

## 7. 有界非正常净价的统一严格折扣

沿用固定正常窗及全域强制性
\[
G(b,x,y,\zeta)=b/27-bxy+D(x,y,\zeta)
\ge c_\nu(b+D)
\tag{26}
\]
（在正常窗外），其中 \(D\) 是两个截幅有向端的完整三次正部成本。其常数在 \(T,m\) 之前固定。正常集合的判据是原参数窗，没有额外要求端点属于好根，故所有深处非正常边均可用此式。

将 \(\mathsf M=\mathsf M_2+\mathsf M_3\) 按最浅轴在／不在目标 pair 中分组。对前者，每个有向端的实际排序积分是
\[
\frac{\lambda_v S_i^\perp T_j^\perp}
 {2naR_i^{n-1}}P_lP_k(P_l+P_k).
\tag{27}
\]
这里 \(S_i^\perp\) 是源完整有效横宽，\(T_j^\perp\) 是目标原始整横宽。\(\mathsf M_3\ge0\) 保留。取 \(c_B<1\) 后，
\[
I-\mathsf M+c_B\mathsf M
=I-(1-c_B)\mathsf M_2-(1-c_B)\mathsf M_3
\le I-(1-c_B)\mathsf M_2.
\]
因此不需要用一个随 \(T\) 增长的比较把 \(\mathsf M_3\) 偷换成零；它自动有利于所需上界，并仍保留在左边 \(c_B\mathsf M\) 中。

将 (27) 的宽度替换成 \(vR_i,vR_j\) 时，正部三次式在 \(|p|,|q|,|z|\le CT\) 上至多 \(CT^3\)。归一化横宽属于固定紧区间，其有限乘积的导数一致有界。固定度数遂给绝对宽度误差
\[
CT^3v^6X_\xi\le CT^3v^{17/2}S.
\tag{28}
\]
实际相交根比造成的误差为 \(Cv\mathsf M\le Cv^9S\)，而不是把 \(CT^3v^7S\) 当误差。接触横权的截幅超额
\((\widehat b-1)_+\le Cv(\xi_i+\xi_j)\)
给 \(CT^2v^7X_\xi\le CT^2v^{19/2}S\)，被 (28) 支配。

分取 (26) 的两个固定小份额支付角权折扣和完整 \(\mathsf M_2\) 预留，得到与 \(T,m\) 无关的 \(\gamma_\nu,c_B>0\)，使
\[
\sum_{e\in\mathcal B}I_e-\mathsf M+c_B\mathsf M
\le(1-\gamma_\nu)A_6W_B
 +C[T^3v^{17/2}+v^9]S.
\tag{29}
\]
整个 \(b+D\) 亏损没有使用两遍；两个固定份额之和小于原强制性常数。

## 8. 近零、长角、边界及二维桥

### 8.1 近零直接 gross

设深处 \(\max(p,q)\le T\)、\(\min(p,q)<m\)。对 \(C_s=C_0/m\) 有逐边
\(pq\le Tm\,\omega_e\)，包括权饱和与未饱和两支。因此由 (11) 和 (1)
\[
\sum_{\mathcal D_{T,m}}I_e
\le CTm\,v^6(D_N+m^{-1}v^3S)
\le C[Tm\,v^8+Tv^9]S.
\tag{30}
\]
没有为极薄目标申请任何源名义成本下界或第二份 \(J\)。

### 8.2 长角的真正统一尾

只用固定旧正常层和其固定小管 \(K_N\)。取大根所有者，在原尺度写
\(P=\max L/r_i,Q=\min L/r_i,H=\min(Q,PQ/v^2)\)。原任意短 \(Q\) 的帽比例界是
\[
|Y_N\cap\mathcal P|/|\mathcal P|\le C(v^3/P+v^2).
\]
当 \(p_e>T\ge T_0\) 时，前项至多 \(C/T\)，固定大 \(T_0\) 后可用原实际条。目标原 raw 例外总 gross 为 \(Cv\mathcal R_*\)，面对／端移例外为 \(Cv(\mathcal F+\mathfrak E)\)，常数不随 \(T\) 增长。

原实际缺源条的成本满足
\[
(v/a)|\mathcal R_e|\ge cr_i^nPQ\min(1/v,P/v^3).
\tag{31}
\]
这里必须使用 (7)：\(P/v^3\ge cT\) 且 \(1/v\ge T/c_0\)，故两个分支均至少 \(cT\)。若遗漏 \(Tv\ll1\)，不能直接由 (31) 声称任意增长 \(T\) 的 \(C/T\) 尾。

条沿长轴的深度至少 \(cP\)，避开固定正常小管；另轴可任意薄。全部目标原整宽截断保留固定比例，原同一见证捕获完整自身目标和异源目标，固定源度数控制重复。无需避开新大 packet 的管，得到
\[
\sum_{\mathcal L_T}I_e
\le\frac C T\mathcal R_*+Cv(\mathcal R_*+\mathcal F+\mathfrak E)
\le C[v^8/T+v^9]S.
\tag{32}
\]

### 8.3 边界与二维

沿用实际 \(\varepsilon_i|F_i|\le|\mathscr D_i|\)、\(|F_i|\asymp aR_i^n\) 及旧正常层体积 \(Cav^5R_i^n\)，边界邻端满足
\[
\sum_{\operatorname{Adj}(\partial I)}\varepsilon_iR_i^n
\le Cv^{10}S+C\mathcal R_*/v=O(v^7S).
\]
原二轴单源 gross 因此为 \(O(v^9S)\)，完全不含 \(T,m\)。未按非好根权粗删深处目标。

二维唯一轴对为 \(\{1,2\}\)，所有最浅轴都在 pair 中；\(\mathsf M_3=\mathfrak L_B=0\)，原始横积为空积一。两方向根差截幅、全部原半径轴帽以及固定正常管长条证明仍逐式成立。槽收费使用第 3 节的二维实际桥。高支撑为空；这里不涉及二维盒内小半径源的独立八阶价格。

## 9. 一次槽转换与总误差

把 (25)、(29)、(30)、(32) 和边界式相加，只减一次完整 \(J\)，得
\[
H_2+c_B\mathsf M+c_R\mathcal R_A^-+c_E\mathfrak E
\le\mathsf P_N+(1-\gamma_\nu)A_6W_B
 +C\mathcal E(v,T,m)S.
\]
在 \(T\ge1,m\le1,v\le1\) 下，未单列的 \(Tv^9,v^9,T^3v^{11},T^2v^{19/2}\) 分别被 (9) 的 \(Tm^{-3}v^9,m^{-1}v^9,T^5v^{10},T^3v^{17/2}\) 支配。

最后由 (11)、(5)、(2)，
\[
\begin{aligned}
\mathsf P_N+(1-\gamma_\nu)A_6W_B
\le{}&q_N\mathcal C+\mathsf E_N\\
&-[q_N-(1-\gamma_\nu)A_6]D_N
 +Cm^{-1}v^9S.
\end{aligned}
\]
方括号至少为固定 \(c_Dv^6\)，即得 (8)。这是唯一一次槽欠量转换；没有再加正常模型已转换后的第二份 \(D_N\) 收入，也没有从本次面空隙收入重新购买 \(\mathcal F\)。

## 10. 显式幂次选择

令
\[
T(v)=v^{-1/8},\qquad m(v)=v^{1/4},\qquad C_s(v)=C_0v^{-1/4}.
\tag{33}
\]
对足够小 \(v\)，固定正常窗包含条件和 (7) 同时成立。七项误差的绝对指数如下。

| 误差 | 代入后的指数 |
|---|---:|
| \(v^8/T\) | \(65/8\) |
| \(Tm v^8\) | \(65/8\) |
| \(Tm^{-3}v^9\) | \(65/8\) |
| \(T^5v^{10}\) | \(75/8\) |
| \(T^2v^{11}\) | \(43/4\) |
| \(T^3v^{17/2}\) | \(65/8\) |
| \(m^{-1}v^9\) | \(35/4\) |

吸收参数 \(Tv,Tm^{-3}v,T^5v^2\) 的指数依次为 \(7/8,1/8,11/8\)，均严格为正。所以候选式 (8) 给
\[
H_2+c_Dv^6D_N+c_B\mathsf M+c_R\mathcal R_A^-+c_E\mathfrak E
\le q_N\mathcal C+\mathsf E_N+O(v^{65/8}S).
\tag{34}
\]
令 \(H=\sum_{i<j}I_e-\sum_iJ_i\)。R061 的固定旧大管、原正常集合及稳定常数正是 (1) 的对象，其结论是数值 gross 上界，不需从 (34) 再消费 \(Z,J,D_N\) 或面预算。二维空和为零，维数至少三时 \(43/5-65/8=19/40>0\)，故得到本稿的全补集候选接口
\[
\boxed{
H+c_Dv^6D_N+c_B\mathsf M+c_R\mathcal R_A^-+c_E\mathfrak E
\le q_N\mathcal C+\mathsf E_N+O(v^{65/8}S).
}
\tag{35}
\]

本账本中 \(\delta=1/8\) 同时平衡长尾、近零、原正常层退款及粗宽度误差。它不是关于实际最优余项的断言。事实上即使只把宽度指数 \(17/2\) 提升到 \(19/2\)，其余三项仍迫使账本最优 \(\delta\le1/8\)：若 \(T=v^{-\alpha},m=v^\beta\)，就须 \(\delta\le\alpha\)、\(\delta\le\beta-\alpha\)、\(\delta\le1-\alpha-3\beta\)，从而 \(1\ge8\delta\)。后续要改进本账本的正幂，需要同时处理这些具体费用，而非仅改善宽度供货。

## 11. 一致性、可测性和未宣称的全局内容

每个 \(v\) 上的 packet 由原有限边清单的实数不等式选取，chart 和目标帽是原 Borel 盒。所有求和均有限，Tonelli 及同一实际物理点的固定逆度给出各积分界。\(T,m\) 不改变根、原半径、目标位置、原分离、CDF 或源／目标排他性。

几何中唯一来自局部匹配的要求，是第 3 节已有固定窗口达到一个固定精度，及原固定稳健面路径的匹配裕量。没有出现匹配误差需要 \(o(mv^3)\)、\(o(T^{-1})\) 或任何随 \(v\) 变化的新精度。原竞争输入序列的 LT 接口因而仍足够；它的未知收敛速度不作为数值误差进入 (9)。所有数值常数独立于该序列，所需固定匹配阈值最终成立。若以后证明全输入的统一 \(v_0\)，可沿失败输入序列使用原 LT 接口；不能在本稿中更换该序列的 \(I\) 或正常模型来优化误差。

本稿目前的完成范围是：在已签同族稳定条件下，构建一个可审的增长测试窗不等式，并接入已签 R061 得到全部交对的明确正幂。它没有改善正常模型 \(\mathsf E_N\)、二维 inside 项或八层前缀的误差。因此 (35) 不能单独宣称 \(C_n\) 的 \(O(av^{65/8})\) 总余项，更不能宣称九阶系数。

## 12. 计算核验与状态

已实际读取 Terra 本机 Mathematica 15.0 ARM 成功脚本及日志，PID 33307，未重跑同一批次：

* `/Users/zhengzhihao/Desktop/T/topic1-tauberian-work/r062-computation/ledger.wl`，SHA `ab3ea93b133d11d9e697f3eeb46f35e2f6e346b96d73e92860ff6ce4689e340c`；
* 同目录 `ledger.log`，SHA `c6cf7067fc55c1ed531ee987bf12427ea6b91be3c264bc5407efd87b22f4fb05`；
* 同目录独立有理核 `ledger_fraction.json`，SHA `62cc0ab757b3e8e0d3c3a5f041cc77bd0edb3daa31855ded6fbe556b2da0f88a`。

八项账本包括 (9) 七项和被支配的 \(Tv^9\)，实际输出指数为
\(\{65/8,65/8,65/8,75/8,43/4,65/8,35/4,71/8\}\)，最小值 \(65/8\)。精确 `Maximize` 返回 \(\delta=1/8,\alpha=1/8,\beta=1/4\)，`Reduce` 对 \(\delta>1/8\) 返回 `False`；只改善宽度一项后的结果仍相同。代码中的一般 guard 测试不代替 (7)，其三个实际吸收指数已在第 10 节逐项给出，严格为正。

计算只核有限指数账本；同一原见证几何、常数增长、正预算及全补集范围由第 3–11 节提供。冻结后交独立审核；正式独立审核前状态为 **PROOF CANDIDATE**。
