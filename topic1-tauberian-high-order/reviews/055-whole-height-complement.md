# R055 整源补集价格：独立数学审查

2026-11-05。本报告仅写入指定审查目录；作者源、前轮冻结稿及仓库均未修改。

## 1. 冻结对象和裁决

全文独立审查的源为 scratch/r055-eighth-complement/whole_height_complement.md，443 行、18643 字节，SHA256：

**0684d01f3665a5a064b4786b6eb1d78ea5780ccdae64f96df3ad01647224a1a1**。

**裁决：PROOF-VALID，无必要数学补丁。** 该源在所列固定 \(n\ge3\) 的同层八阶竞争族内证明 (4) 和中心化形式 (6)。整源高度所引入的新几何、固定窗成本上极限、补集尾估计和一次完整 \(J\) 分配均闭合。正常模型价格仍是明确定义的实际泛函 \(\mathsf P_N\)，本报告不使用 R053 的待审系数或 R054 的 \(P_I\) 改进，也不签署 \(C_n\) 的精确八阶极限。

主要实际读取并核合的依赖为：

* R049 源，SHA **5c4d7b371a0c33c39449afd686de599c49570aaa639d4f7db7d3c1f166a040f4**；
* R050 全文源，SHA **22f848879b3c9ba839edb5d307a205ae3245e5ba3a89a99f7ec2ce2974010748**，及其正式全文审报 scratch/r050-review/audit.md，SHA **68e9532a8b80926806b0ed11219c7630a46fe7aa4a9d72c0a5e94a73b520377e**；
* 已由本审查者完整独审的 R052 源 **d655c1e160c88914774248c81688c6bdde3f27614387c4eb7f0b9c1a0a0a5296**，并复核 R051 已审原见证高度接口；
* R045 的共同实际层、全域非正常强制性、任意长短比的实际帽、两管长条和槽占用公式；本次读取源 SHA **8f0ad8b01bb443789b33f55188c7a04d000a3ceaee3fd61d0a83eb458c9bf8ea** 的相关全部证明；
* R036 全文 adaptive-clean 分配，SHA **c6b761427b26c1ab88d2ecf17f47ef856cb7c9565d53dc980bcb78fbd6e75c3d**，及 R042 强 CDF、实际宽度误差的原定义。

## 2. 稳定常数与新窗的次序

正常集合、旧正常实际层、正常 clean 分配、\(D_N\) 和 R049 的一次稳定预算份额均先固定。R049 的原固定构造供应
\[
 \mathcal R_*+\mathfrak E+\mathcal F+\mathcal F_{1/128}
 \le C_*v^8S,\qquad D_N\le C_*v^2S,\qquad
 S(\partial I)=O(v^5S),
\]
其中 \(\mathcal R_*=(v/a)\sum|\mathscr D_i\setminus Y_{i,N}^{\rm old}|\)。\(C_*\) 不通过后选 \(T,m\) 的新价格重新证明。因此用此常数控制新窗、长尾和近零角没有阈值循环。

新 packet 的成员是同一实际外部族中两端均深处的固定有界二轴目标；实际根、原见证、原始／有效／中点盒不改变。R051 的侵蚀和 R052 的小算子常数可依新固定 \(T\)，但最终小 \(v\) 才使用这些常数。未使用的源仍保留其全部 raw 缺量，未从外部族删除。

## 3. same-eroder 的常数一 raw 转移

在整个有效源盒的几乎处处不交 global-min chart 中，实际高度均非负：
\[
 h_A=h_N+h_B,\qquad h_N,h_B\ge0.
\]
非负性来自帽内真实 \(F_j\) 质量，而非名义高度减一个可能任意符号的误差。目标扣量仍在原始目标盒，源始终在有效盒。取全 packet 的同一个非负侵蚀 \(\epsilon_i(y)\)，同时用于两个高度。R051 对 \(h_A\) 的必要性自动给较小 \(h_N\) 的同一保证，所以 \(B_N^-\subset B_A^-\subset D'\)。

对任意 \(H,u,w,\epsilon\ge0\)，
\[
 0\le\min(H,(u+w-\epsilon)_+)-\min(H,(u-\epsilon)_+)\le w.
\]
这在零高度和饱和分支也成立。逐 chart 积分后，
\[
 \mathcal R_N^--\mathcal R_A^-
 =\frac va|B_A^-\setminus B_N^-|_1
 \le\frac va\sum\int h_B^0
 =(1+v)^{-(n-1)}\mathsf M_{T,m}.
\]
这是精确系数一，不需要正常／非正常层分开，也不依赖侵蚀常数比较。

旧 \(Y_N^{\rm old}\subset B_N(\mathscr D)\) 是 R051/R052 同有效源坐标下的精确锚定；旧核心和旧侵蚀没有重新选择。侵蚀物理损失 \(O_T(av^9S)\) 乘 \(v/a\) 为 \(O_T(v^{10}S)\)，所以
\[
 \mathcal R_A^-\le\mathcal R_N^-\le
 \mathcal R_*+O_T(v^{10}S).
\]
所有集合包含按零测高度图形的通用约定理解，不影响任何体积式。

## 4. 真实上截断只损失相对 \(O_T(v^3)\)

在一条纤维上，
\[
 \min(H,h_N+h_B)-\min(H,h_N)=h_B-\ell,\quad
 0\le\ell\le h_B\,\mathbf1_{\{H<h_N+h_B\}}.
\]
因此 \(\Delta_A=\mathsf M_{T,m}-\mathfrak F_B-\mathfrak L_B\) 是精确恒等式，且可用 \(h_B^0\) 上包损失。

对名义坏目标的固定项，正帽只约束它的两个临界 base 变量为 \(y_l,y_m=O_T(v^3R_i)\)。在此支撑内，临界轴的 \(w'_{ih}/2-y_h\) 与最浅轴的 \(w'_{iq}/2\) 都是 \(vR_i\) 的固定正倍数，大于总高度 \(C_Tv^4R_i\)。故饱和只能由某个不在该目标轴对内、且不等于 \(q\) 的自由 base 变量靠近其区间上端造成。

该自由变量范围为 \(w'_{ih}/2\asymp vR_i\)，危险带宽至多 \(C_Tv^4R_i\)。固定目标的名义项独立于它，故 Fubini 给危险带相对积分至多 \(C_Tv^3\)。对有限个非临界自由变量、目标和 chart 取并，得到
\[
 0\le\mathfrak L_B\le C_Tv^3\mathsf M_{T,m}.
\]
若没有这样的自由变量，根本没有该项的饱和损失。此论证使用真实 base 和真实宽度，不把整个横向核心裁掉。特别最浅轴不在目标轴对的图仍完整存在；其三近面主体未以误差名义删除。

## 5. 新窗名义坏成本的八阶界与统一上极限

旧正常层精确包含于 \(B_N\)。因此实际坏增量中位于真实缺集内的部分属于 \(\mathscr D\setminus Y_N^{\rm old}\)，位于真实缺集外的部分由 R051 必要性控制：
\[
 \Delta_A\le\lambda_v r_*+O_T(v^{10}S)
 =(1+v)^{n-1}\mathcal R_*+O_T(v^{10}S).
\]
该估计先于任何关于新 \(\mathsf M_{T,m}\) 的八阶断言，未借用旧窗的坏成本来替代新窗。

坏退款的全族实际帽核范数为 \(C_Tv^2/a\)。目标缺量按固定旧正常层和其余量分解，后者贡献至多 \(C_Tv\mathcal R_*\)。旧正常层逐目标体积 \(O(av^5R_j^n)\)，所以每条固定有界坏边的退款为 \(O_T(v^7R_e^n)\)。

另一方面，当 \(p,q\ge m>0\)，两个方向的实际完整名义成本之和至少为 \(c_m v^6R_e^n\)。原始目标整横宽、有效源整横宽均有统一正下界；根差截幅的两方向凸性保证至少这个成本，未要求坏根差未截断。有限度求和给
\[
 \mathfrak F_B\le C_{T,m}v\mathsf M_{T,m}+C_Tv\mathcal R_*.
\]
此式没有用于近零角。

结合精确增量式及上一节，
\[
 [1-C_{T,m}v-C_Tv^3]\mathsf M_{T,m}
 \le[(1+v)^{n-1}+C_Tv]\mathcal R_*+O_T(v^{10}S).
\]
每个 \(T,m\) 固定后先取小 \(v\)，左系数趋一。因此
\[
 \mathsf M_{T,m}=O_{T,m}(v^8S),\quad
 \mathfrak F_B,\mathfrak L_B=o_{T,m}(v^8S),\quad
 \limsup_{v\to0}\frac{\mathsf M_{T,m}}{v^8S}\le C_*.
\]
最后这个常数确与后选固定 \(T,m\) 无关；并未声称其小参数阈值或误差常数一致。作者 (17) 的量词正确。

## 6. 一次完整 CDF 的可加性

R052 对正常 packet 单独给
\[
 \lambda_v|B_N(\mathscr D)|_1
 =\lambda_v|Y_{2,N}|_1+O(v^9S).
\]
这只是模型层替换，没有减去任何 \(J\)。加上实际增量及上一节小量，得到
\[
 \lambda_v|B_A(\mathscr D)|_1
 =\lambda_v|Y_{2,N}|_1+\mathsf M_{T,m}+o_{T,m}(v^8S).
\]
对全 packet 的同一个真实保证层 \(B_A^-\) 使用强 CDF 一次，侵蚀成本 \(O_T(v^{10}S)\)，鲜导数相对 epsilon 损失 \(O_T(v^{11}S)\)，于是
\[
 \sum_iJ_i\ge
 \lambda_v|Y_{2,N}|_1+\mathsf M_{T,m}
 +c_R\mathcal R_A^-+c_E\mathfrak E-o_{T,m}(v^8S).
\]
原强 CDF 同时分配层成本、raw 和端移；常数可在新 \(T,m\) 之前固定。没有把正常减全 \(J\) 的价格与另一条补集源价格相加。最终保留的正常项确为
\(\mathsf P_N=\sum_{\mathcal N}I_e-\lambda_v|Y_{2,N}|_1\)，不是再次减全体 \(J\) 的正常净价。

## 7. 全横宽积分、三近面项与固定非正常折扣

固定目标轴对，最浅轴为这两轴之一的两个排序积分分别为 \(P_lP_m^2/2\)、\(P_mP_l^2/2\)。对其他最近面符号求和后，源 base 的半宽区间恰拼成完整有效源横宽 \(S_i^\perp\)。故作者 (21) 的每个有向端点因子
\[
 \frac{\lambda_vS_i^\perp T_j^\perp}{2naR_i^{n-1}}
 P_lP_m(P_l+P_m)
\]
正确，未漏 \(n\) 维横截面或排序因子。

若某个额外轴 \(q\) 最浅，两个 \(q\) 面与两个临界 base 的积分相对上述量的比率为
\[
 \frac{P_lP_m}{(P_l+P_m)w'_{iq}}\le C_Tv^2.
\]
零帽时两边直接为零，不除以零。因此 \(\mathsf M_3\ge0\)，且 \(\mathsf M_3\le C_Tv^2\mathsf M_2\)。这是真实纳入名义成本的三近面项，不是把高支撑 gross 小量误用于源层的三近面交叠。

实际整宽相对 \(vR_i\) 的误差仍用原 \(\xi_i\) 量。固定窗下总宽度误差 \(O_T(v^6\sum\xi_iR_i^n)=O_T(v^{17/2}S)\)；原实际相交根比 \(1+O(v)\) 只产生 \(Cv\) 倍完整名义成本，现已为 \(o(v^8S)\)。接触权截至一的误差也更晚。旧核心所致的相对 \(O_T(v^2)\) 名义损失已不在这里；不能把它留着再改称 little-o。

在归一化变量 \(x=(1+v)p/n,\ y=(1+v)q/n\) 下，理想交量／成本共同因子是 \(n^2(1+v)^{n-4}v^6R_e^n\)。R045 的全域
\[
 G=b/27-bxy+D\ge c_\nu(b+D)
\]
在正常窗外成立，\(c_\nu\) 不依 \(T,m\)。分出固定角权折扣及固定完整成本份额，结合上述 little-o 宽度误差，确得
\[
 \sum_{\mathcal B_{T,m}}I_e-\mathsf M_{T,m}
 +c_B\mathsf M_{T,m}
 \le(1-\gamma_\nu)A_6(v)W_B+o_{T,m}(v^8S).
\]
\(c_B,\gamma_\nu\) 可先固定；\(C_Tv^2\) 的三近面比较只改变最终小 \(v\) 阈值，不使这些保留系数依赖 \(T,m\)。

## 8. 近零 gross、长条和全边界

对于固定窗内 \(\min(p,q)<m\)，\(\omega=\min(1,C_{\rm soft}\min(p,q))\)、\(C_{\rm soft}\ge1/m\) 给 \(pq\le Tm\omega\)：未饱和支用 \(1/C_{\rm soft}\le m\)，饱和支用 \(\min(p,q)<m\)。正常边及有界坏边软权均为一。

R036 的 clean 阈值为 \(\ell_e/(1000r_i)\)，仅依真实交长，与后选软权系数无关。probe、helper、Borel 清单片和正常占用并可保持原样；改变 \(C_{\rm soft}\) 只改变碰撞收费常数。故同一个正常 \(D_N\) 满足
\[
 W_{\mathcal D},\,W_B\le\mathcal C-W_N+O_{C_{\rm soft}}(v^3S)
 =D_N+O_{C_{\rm soft}}(v^3S).
\]
其中 \(W_{\mathcal D}\) 使用其软权。于是近零 gross 的归一化上极限至多 \(TmC_*\)。没有对极薄帽使用坏退款中依赖 \(m>0\) 的估计。

长角使用固定旧正常小管。原实际目标帽中的旧正常层比例 \(C(v^3/P+v^2)\) 统一于任意小短边；目标 raw 例外只花其实际数值 \(\mathcal R_*\)，面对／端移例外只用其实际数值 \(\mathcal F+\mathfrak E\)。其余实际缺源条避开固定旧正常层，成本／交量比至少 \(cT\)，实际截宽保留固定比例，同一原指定见证及统一度数处理所有条。因此
\[
 \sum_{\mathcal L_T}I_e\le (C_L/T)\mathcal R_*
   +Cv(\mathcal R_*+\mathcal F+\mathfrak E).
\]
\(C_L\) 只依旧正常小管和原实际常数，独立于新 \(T,m\)。这些是已稳定的几何量给出的 gross 数值界，不从本轮 (19) 再抽一份 \(J\) 或 raw 正收入。

边界 epsilon 的正确方向来自实际 \(\varepsilon_i|F_i|\le|\mathscr D_i|\)，再按旧正常层分解，给
\[
 \sum_{\operatorname{Adj}(\partial I)}\varepsilon_iR_i^n
 \le Cv^5S(\partial I)+(C/v)\mathcal R_*=O(v^7S).
\]
原边界二轴单源界因此为 \(O(v^9S)\)。R050 已审源／正式审报明确在同一固定竞争构造上供应全部 \(q\ge3\) gross \(o(v^8S)\)，包括外部边界；其使用独立的尾份额，不再扣一份 \(J\)。原 global 双向面容量排除单临界正体积交，零临界由中心分离不存在。故分组覆盖完整 \(H\)，不按非好根权丢弃任何实际交对。

## 9. 合并、中心化和最终量词

给定 \(\delta>0\)，先用固定 \(C_*,C_L\) 取大而有限的 \(T\)，再取固定小 \(m\) 使 \(TmC_*\) 足够小，随后固定 \(C_{\rm soft}\)，最后取 \(v\to0\)。长尾和近零 gross 共同置于 \(\delta v^8S\)；其余误差是对这些固定参数的 little-o。与第 6、7 节相加便得 (4)，完整 \(J\) 只出现一次。

正常模型、\(W_N,D_N\) 和任给 \(q_N=A_6+O(v^7)\) 在此过程不改变。由 \(W_B\le D_N+O(v^3S)\) 及 \(W_N=\mathcal C-D_N+O(v^3S)\)，
\[
 \mathsf P_N+(1-\gamma_\nu)A_6W_B
 \le q_N\mathcal C+\mathsf E_N
 -[q_N-(1-\gamma_\nu)A_6]D_N+O(v^9S).
\]
括号为 \(\gamma_\nu k v^6+O(v^7)\)，固定正下界 \(c_Dv^6\) 成立，得到 (6)。真实 \(D_N=O(v^2S)\) 的八阶领先损失明确保留，未塞入未知误差。

对当前 \(n\ge3\)，原小源盒内扩张／环境绝对余项 \(v^{n+6}S\) 最早九阶；盒外 \(Cv\mathcal F\) 亦为九阶。此判断只核局部接口能接原总装的误差阶，不代替同族前缀望远镜或完整全输入上确界证明。二维未被纳入这一升级。

## 10. 实际本机计算与最终签署

实际读取 r055-computation/clip_increment.wl、clip_increment.log，SHA256 分别为 **274bba84c1a184907a6e118ad9193ab2b5c05e725f38302767d775da572462ea**、**6b95cb0c57a112a139e44d6632ac64854b24bd66271506e611ab693237a373d4**。成功日志的两个违反集为 False；首次无效路径日志不作证据。这仅核截断标量。

本次另实际调用本机 WolframScript 15，核两个排序积分及额外最浅轴比率。实际输出为：

    two-order integral residual: 0
    extra-axis ratio residual: 0
    ratio upper bound violation: False
    Null

具体核验的是
\(\int_0^q p(q-y)\,dy+\int_0^p q(p-x)\,dx=pq(p+q)/2\)，以及两侧额外浅面积分与此量之比 \(pq/[(p+q)w]\le\min(p,q)/w\)。这些代数不证明真实 chart 的覆盖、非临界中分面危险带、原指定见证或强 CDF 分配；相应解析审查已在上文给出。

**最终签署：作者冻结 SHA 0684d01f3665a5a064b4786b6eb1d78ea5780ccdae64f96df3ad01647224a1a1 的 (4)、(6) 在明确固定 \(n\ge3\) 的同层竞争输入范围内 PROOF-VALID，无必要修补。** 本报告不依赖 R053 的具体正常八阶系数或 R054 的 \(P_I\) 改进，不宣布 \(J_n\) 上系数或真实八阶极限已闭合。
