# R055：八阶竞争族中全部非正常补集的整源共同结算

2026-11-05。固定 \(n\ge3\)。新局部证明候选，供独立审核。仓库、论文及前轮冻结源只读。

目标是把补集的旧 \(O(v^8S)\) 未定常数去掉，同时把正常价格保留为一个明确的实际模型泛函。本文不计算 R053 的正常八阶积分，也不宣布精确八阶系数。全部源成本仅由一个完整 packet、同一份强 CDF 盈余 \(J\) 结算一次。

## 1. 固定实际族、正常模型和稳定常数

在同一八阶竞争输入序列、同一新固定 \(N=8\) 族上，使用已审 R049/R050/R051/R052。记
\[
 a=v^n,\quad S=\sum r_i^n,\quad T_{\rm amb}=V+M,\quad
 k=n^2/27,\quad A_6(v)=\frac{k}{(1+v)^2}v^6,\quad
 \lambda_v=(1+v)^{n-1}v/a .
\]
所有实际外部源、目标根、CDF、正常集合 \(\mathcal N\)、旧正常实际层 \(Y_{i,N}^{\rm old}\)、正常 clean 分配和槽欠量 \(D_N\) **先固定**，不随后选的补集阈值改变。

记
\[
 \mathscr D_i=(x_i+W_i^0)\setminus F_i,\qquad
 r_*=\sum_i|\mathscr D_i\setminus Y_{i,N}^{\rm old}|,\qquad
 \mathcal R_*=(v/a)r_*.
\]
已审 R049 给一个固定 \(C_*<\infty\)，使
\[
 \mathcal R_*+\mathfrak E+\mathcal F+\mathcal F_{1/128}
       \le C_*v^8S,\qquad
 D_N\le C_*v^2S,\qquad
 S(\partial I)=O(v^5S),
\tag{1}
\]
\[
 \mathfrak E=\sum_iR_i^{n-1}\Delta_i^{\rm end}.
\]
\(C_*\) 来自该次先固定构造和其稳定预算份额，不随后面的 \(T,m\) 重新选择。可用给定八阶竞争常数 \(B\) 下的 R049 统一上界；也可沿一条竞争序列取固定有限上界。其他粗先验及有限度／相交根比不变。

在 R052 的整个有效源盒 global-min-face 坐标中，仅取正常 packet，定义模型算子 \(B_N\) 和
\[
 Y_{0,N}=B_N(0),\quad Y_{1,N}=B_N(Y_{0,N}),\quad
 Y_{2,N}=B_N(Y_{1,N}).
\]
源始终用有效盒 chart，目标缺量始终用原始盒；全部几何、实际根、clip 和 chart 不随测试缺集改变。定义正常**模型**净价格
\[
 \boxed{\mathsf P_N(v)=\sum_{e\in\mathcal N}I_e
                         -\lambda_v\sum_i|Y_{2,N,i}|.}
\tag{2}
\]
它没有减去全体 \(J_i\)。R053 的正常模型积分可以给 (2) 的上界，本稿只把它作为实际数值保留。

令
\[
 W_N=\sum_{\mathcal N}b_eR_e^2,\qquad
 \mathcal C=\frac14\sum_{\text{实际槽 }s}R_i^2|F_s|.
\]
同一正常 clean 占用给
\[
 W_N=\mathcal C-D_N+O(v^3S).
\tag{3}
\]
正常两交长在一个固定 \([c_Nv^3R_e,T_Nv^3R_e]\) 中；它们在下面所有软权中均取权一。

## 2. 本稿结论及抽象正常接口

存在固定 \(\gamma_\nu,c_B,c_R,c_E>0\)，满足下述命题。给任意固定 \(\delta>0\)，先取有限大 \(T>T_N\)，再取有限小 \(m<c_N/2\)，在两端均深的二轴边中定义
\[
 \mathcal B_{T,m}=
 \{e\notin\mathcal N:m\le p_e,q_e\le T\},\qquad
 p_e=L_1/(v^3R_e),\quad q_e=L_2/(v^3R_e).
\]
用 \(\mathscr P=\mathcal N\cup\mathcal B_{T,m}\) 的**同一个整源联合高度**作强 CDF。令 \(\mathsf M_{T,m}\) 为第 3 节的完整非正常名义高度成本，\(\mathcal R_A^-\) 为这个 packet 的真实保证层之外的 raw 成本。写
\[
 W_B=\sum_{\mathcal B_{T,m}}b_eR_e^2,\qquad
 H=\sum_{i<j}I_e-\sum_iJ_i
\]
（全部外部交对，包括边界和高支撑）。则
\[
 \boxed{\begin{aligned}
 H+c_B\mathsf M_{T,m}+c_R\mathcal R_A^-+c_E\mathfrak E
 \le \mathsf P_N(v)
       +(1-\gamma_\nu)A_6(v)W_B
       +\delta v^8S+o_{\delta}(v^8S).
 \end{aligned}}
\tag{4}
\]
常数 \(\gamma_\nu,c_B,c_R,c_E\) 只用固定正常窗及原实际源常数，可在 \(T,m\) 之前固定；起效阈值和小量允许依赖固定 \(T,m,\delta\)。

更方便的容量形式如下。任取不随后面阈值改变的标量
\[
 q_N(v)=A_6(v)+O(v^7),\qquad
 \mathsf E_N(v)=\mathsf P_N(v)-q_N(v)W_N .
\tag{5}
\]
对某个固定 \(c_D>0\)，(4) 推出
\[
 \boxed{\begin{aligned}
 H+c_Dv^6D_N+c_B\mathsf M_{T,m}
          +c_R\mathcal R_A^-+c_E\mathfrak E
 \le q_N(v)\mathcal C+\mathsf E_N(v)
       +\delta v^8S+o_\delta(v^8S).
 \end{aligned}}
\tag{6}
\]
因此若 R053 独立交付一个**同一正常模型**上界
\(\mathsf P_N\le q_NW_N-\mathsf L_N+o(v^8S)\)，
其中 \(\mathsf L_N\ge0\)，则 (6) 保留 \(\mathsf L_N\)，并在
先 \(v\to0\)、再 \(\delta\downarrow0\) 后没有新增八阶常数。
该推论不假设 R053 的任何待审系数已经成立。

## 3. 整源共同高度与常数一的 raw 转移

对固定 \(T,m\)，采用 R051/R052 的完全相同有效源 chart 和真实纤维上截断 \(H_i(y)\)。将 packet 的名义与实际高度分别写成
\[
 h_A^0=h_N^0+h_B^0,\qquad
 h_A(\mathscr D)=h_N(\mathscr D)+h_B(\mathscr D),
\tag{7}
\]
\[
 h_N=h_N^0-\mathcal K_N(\mathscr D)\ge0,\qquad
 h_B=h_B^0-\mathcal K_B(\mathscr D)\ge0.
\]
这里每项都是实际目标帽中 \(F_j\) 质量除以 \(naR_i^{n-1}\)，所有目标先相加；图包括最浅轴不属于该目标轴对的部分。定义
\[
 B_A(\mathscr D)=\text{height }\min(H_i,h_N+h_B),\qquad
 B_N(\mathscr D)=\text{height }\min(H_i,h_N).
\]
未作旧核心剪裁。

完整非正常名义成本定义为
\[
 \boxed{\displaystyle
 \mathsf M_{T,m}=\lambda_v
       \sum_{\rm all\ charts}\int h_B^0(y)\,dy.}
\tag{8}
\]
积分保留每个目标原始整横宽、每个源实际有效宽和 \(R_i^{n-1}\) 分母，不含接触小权 \(b_e\)。

取 R051 对整个 packet 有效的共同非负侵蚀 \(\epsilon_i(y)\)，同时用于 \(h_N\) 和 \(h_A\)，得到
\[
 B_N^-:\ \min(H_i,[h_N-\epsilon_i]_+),\qquad
 B_A^-:\ \min(H_i,[h_N+h_B-\epsilon_i]_+).
\]
两者都是真实 \(D'_i\) 的子集，且
\[
 B_N^-\subset B_A^-,\qquad
 \sum_i|B_N\setminus B_N^-|+
 \sum_i|B_A\setminus B_A^-|=O_T(av^9S).
\tag{9}
\]
对近零高度也用同一个正部与同一个上截断。令
\[
 \mathcal R_N^-=\frac va\sum_i|\mathscr D_i\setminus B_{N,i}^-|,
 \qquad
 \mathcal R_A^-=\frac va\sum_i|\mathscr D_i\setminus B_{A,i}^-|.
\]
逐 chart 的 1-Lipschitz 单调性给
\[
 0\le |B_A^-\setminus B_N^-|
       \le\sum_{\rm charts}\int h_B\,dy
       \le\sum_{\rm charts}\int h_B^0\,dy.
\]
于是有**精确常数一**的转移
\[
 \boxed{\displaystyle
 \mathcal R_N^-\le
       \mathcal R_A^-+(1+v)^{-(n-1)}\mathsf M_{T,m}
       \le\mathcal R_A^-+\mathsf M_{T,m}.}
\tag{10}
\]
这不要求正常与非正常层物理分离。R051 的精确锚定
\(Y_N^{\rm old}\subset B_N(\mathscr D)\) 又给
\[
 \mathcal R_A^-\le\mathcal R_N^-
       \le\mathcal R_*+O_T(v^{10}S).
\tag{11}
\]
因此整源／真实 clip 没有破坏 R045 需要的 raw 方向。

## 4. 真实纤维 clip 损失是相对 \(O_T(v^3)\)

定义共同高度的非正常实际增量成本
\[
 \Delta_A=\lambda_v\sum_i|B_A(\mathscr D)_i
                                \setminus B_N(\mathscr D)_i|.
\]
因为 \(h_N,h_B\ge0\)，
\[
 \Delta_A=\mathsf M_{T,m}-\mathfrak F_B-\mathfrak L_B,\qquad
 \mathfrak F_B=\lambda_v\sum_{\rm charts}\int
                                   \mathcal K_B(\mathscr D)(y)\,dy,
\tag{12}
\]
其中
\[
 0\le\mathfrak L_B\le
 \lambda_v\sum_{\rm charts}\int
 h_B^0(y)\,\mathbf1_{\{H_i(y)<h_N+h_B\}}\,dy.
\]
各有界高度均不超过 \(C_Tv^4R_i\)。固定一个非正常目标的名义项，它只依赖该目标两个临界底变量，且在这些变量超过 \(C_Tv^3R_i\) 后为零。

在该项为正的地方，源最浅轴半宽以及两个临界轴的
\(w'_{ih}/2-y_h\) 均大于固定倍数 \(vR_i\)，不可能小于 \(C_Tv^4R_i\)。所以 \(H_i<C_Tv^4R_i\) 只能由某个**非临界自由底变量**落在靠近源中分面的宽 \(C_Tv^4R_i\) 薄带造成。该名义项独立于这个自由变量，其全范围为固定倍数 \(vR_i\)，故带占比至多 \(C_Tv^3\)。逐目标、逐图有限求和得到
\[
 \boxed{0\le\mathfrak L_B\le C_Tv^3\mathsf M_{T,m}.}
\tag{13}
\]
若没有非临界自由底变量，该项没有此 clip 损失。这是源盒中分面截断的误差，不是删除三近面主体；最浅轴不在目标轴对的整张图始终在 (8)、(12) 中。

## 5. 新固定窗的非正常完整成本也有八阶界

先不借用 R049 中旧固定窗的 \(B_{\rm bad}\) 来声称新
\(\mathsf M_{T,m}\) 已为八阶。它由下述实际集合估计重新证明。

R051 真实必要性给
\(\sum|B_A(\mathscr D)\setminus\mathscr D|\le C_Tav^9S\)，
而 \(Y_N^{\rm old}\subset B_N(\mathscr D)\)。因此
\[
 \Delta_A\le\lambda_v(r_*+C_Tav^9S)
       \le(1+v)^{n-1}\mathcal R_*+O_T(v^{10}S).
\tag{14}
\]

对坏退款 \(\mathfrak F_B\)，实际帽核的全族范数至多
\(C_Tv^2/a\)。目标缺集按固定的旧正常层及其余量分解。余量部分给
\[
 \mathfrak F_{B,\rm raw}\le C_Tv\mathcal R_*.
\]
正常层逐源体积至多 \(Ca v^5R_j^n\)，故每个有界坏端的退款至多 \(C_Tv^7R_e^n\)。另一方面，\(p,q\ge m\) 及原宽的统一正下界，使每条无向非正常边在 (8) 中的两个方向名义成本至少为
\[
 c_m v^6R_e^n.
\]
根差可截幅：两方向正部凸性仍给
\(D\ge pq(p+q)/n\ge2m^3/n\)，不需要坏根也接近正常。于是
\[
 \boxed{\mathfrak F_B\le
       C_{T,m}v\mathsf M_{T,m}+C_Tv\mathcal R_*.}
\tag{15}
\]
没有对近零角使用这个需要 \(m>0\) 的估计。

结合 (12)–(15)，
\[
 [1-C_{T,m}v-C_Tv^3]\mathsf M_{T,m}
 \le[(1+v)^{n-1}+C_Tv]\mathcal R_*+O_T(v^{10}S).
\]
对每个固定 \(T,m\)，可得
\[
 \mathsf M_{T,m}=O_{T,m}(v^8S),\quad
 \mathfrak F_B=o_{T,m}(v^8S),\quad
 \mathfrak L_B=o_{T,m}(v^8S).
\tag{16}
\]
更重要的是
\[
 \boxed{\displaystyle
 \limsup_{v\to0}\frac{\mathsf M_{T,m}}{v^8S}
 \le\limsup_{v\to0}\frac{\mathcal R_*}{v^8S}\le C_*.}
\tag{17}
\]
这个上极限常数不随后选固定 \(T,m\) 增长；增长的误差常数只乘消失的 \(v\) 或更高幂。这不是对 \(T=T(v)\) 的一致断言。

## 6. 同一份强 CDF 的模型成本可加性

由 (12)、(16)，以及正常 R052 的第二次迭代替换，
\[
 \lambda_v\sum_i|B_A(\mathscr D)_i|
 =\lambda_v\sum_i|Y_{2,N,i}|+\mathsf M_{T,m}
                         +o_{T,m}(v^8S).
\tag{18}
\]
R051 的侵蚀成本 \(O_T(v^{10}S)\)、强鲜导数 epsilon 误差
\(O_T(v^{11}S)\) 均已晚于八阶。对真实保证层 \(B_A^-\subset D'\) 使用强 CDF **一次**，得
\[
 \boxed{\displaystyle
 \sum_iJ_i\ge
   \lambda_v\sum_i|Y_{2,N,i}|+\mathsf M_{T,m}
        +c_R\mathcal R_A^-+c_E\mathfrak E
        -o_{T,m}(v^8S).}
\tag{19}
\]
若先对正常 packet 写一个减全部 \(J\) 的价格，再加补集源成本，将不能推出 (19)；这里的依据是整个 \(h_A=h_N+h_B\) 的真实层。

## 7. 有界非正常净价的严格折扣不再含未知八阶误差

正常窗外的全域标量亏损仍为
\[
 G=b/27-bxy+D\ge c_\nu(b+D),
\tag{20}
\]
其中 \(c_\nu\) 在 \(T,m\) 前固定。将 (8) 按最浅轴是否属于该 target pair 分成 \(\mathsf M_2+\mathsf M_3\)。前者的每个有向端点，带实际源、目标整横宽的积分恰为
\[
 \frac{\lambda_v S_i^\perp T_j^\perp}
           {2naR_i^{n-1}}P_{ij,l}P_{ij,m}(P_{ij,l}+P_{ij,m}),
\tag{21}
\]
此处 \(S_i^\perp\) 是**完整有效源横宽**，没有旧核心剪裁。
后者始终非负，并由一维帽积分直接给
\[
 0\le\mathsf M_3\le C_Tv^2\mathsf M_2.
\tag{22}
\]
具体比率每个额外最浅轴至多是
\(C\,P_lP_m/[(P_l+P_m)w'_{iq}]\le C_Tv^2\)；
零帽自动给零，不除以可能为零的交长。

将 (21) 的宽度换成 \(vR_i,vR_j\) 的平均误差是
\(O_T(v^6\sum\xi_iR_i^n)=O_T(v^{17/2}S)=o(v^8S)\)。
实际相交根比 \(1+O(v)\) 的误差是完整名义成本的 \(Cv\) 倍；由 (16) 亦为 \(o(v^8S)\)。接触横权超过一的截断误差为
\(O_T(v^7\sum\xi_iR_i^n)=o(v^8S)\)。
容量／高度必要性误差已由 R051 给到十阶。

因此可以用 (20) 保留固定角权折扣及一份完整名义成本，得到
\[
 \boxed{\displaystyle
 \sum_{\mathcal B_{T,m}} I_e-\mathsf M_{T,m}
       +c_B\mathsf M_{T,m}
 \le(1-\gamma_\nu)A_6(v)W_B+o_{T,m}(v^8S).}
\tag{23}
\]
\(c_B,\gamma_\nu>0\) 可在 \(T,m\) 前固定；(22) 只需最终小 \(v\) 使其比较因子小于二。这里原旧核心造成的 \(C_Tv^8S\) 常数已经真正消失，不能仅把它改写成 little-o。

## 8. 近零、长角、边界和高支撑

### 8.1 近零有界角用软容量直接控制 gross

设 \(\mathcal D_{T,m}\) 是深处、\(\max(p,q)\le T\)、\(\min(p,q)<m\) 的边。取固定
\[
 C_{\rm soft}\ge\max(1/m,1,12/n),\qquad
 \omega_e=\min(1,C_{\rm soft}\min(p,q)).
\]
直接有 \(pq\le Tm\,\omega_e\)，包括软权饱和和未饱和两支。R036 的全软容量及 (3) 给
\[
 W_{\mathcal D}\le\mathcal C-W_N+O_{C_{\rm soft}}(v^3S)
                  =D_N+O_{C_{\rm soft}}(v^3S).
\]
故
\[
 \boxed{\displaystyle
 \limsup_{v\to0}\frac{\sum_{\mathcal D_{T,m}}I_e}{v^8S}
 \le Tm\,C_*.}
\tag{24}
\]
没有把极薄角纳入 (15)，也没有在其退款上使用一个不存在的统一正高度。

### 8.2 长二轴 gross 的固定窗尾

对深处 \(\max(p,q)>T\) 的边，沿 R045 两管证明，但只用**固定旧正常层**及其固定小管。实际正常层在任意目标长帽中的占比仍为
\[
 C(v^3/P+v^2),
\]
统一于短边 \(Q\)，常数不依 \(T\)。目标例外只用
\(\mathscr D_j\setminus Y_{j,N}^{\rm old}\)，总 gross 为
\(Cv\mathcal R_*\)；面对／端移例外总 gross 为
\(Cv(\mathcal F+\mathfrak E)\)。其余长短实际源条避开固定正常小管，成本与交量之比至少 \(cT\)。条的固定重数给
\[
 \sum_{\mathcal L_T}I_e
 \le\frac{C_L}{T}\mathcal R_*
          +Cv(\mathcal R_*+\mathcal F+\mathfrak E).
\tag{25}
\]
所有穿透均截于真实目标宽而保留固定比例，源点使用同一个原半径；这些是已审 R045 的原实际长条，不要求条避开新的整源大管。

由 (1)，
\[
 \boxed{\displaystyle
 \limsup_{v\to0}\frac{\sum_{\mathcal L_T}I_e}{v^8S}
              \le C_LC_*/T.}
\tag{26}
\]
这里是由已证稳定性获得的数值 gross 界，不再从 (19) 的 \(J\) 抽一份成本来支付同一长条。\(C_L\) 只依固定正常小管和原源常数，独立于新 \(T,m\)。

### 8.3 边界及高支撑

实际 epsilon 分解
\[
 \varepsilon_iR_i^n\le Cv^5R_i^n+
             C|\mathscr D_i\setminus Y_{i,N}^{\rm old}|/a
\]
来自 \(\varepsilon_i|F_i|\le|\mathscr D_i|\) 与
\(|F_i|\asymp aR_i^n\)，不是倒置粗 epsilon 上界。
对边界邻端用有限度、\(S(\partial I)=O(v^5S)\)，有
\[
 \sum_{\operatorname{Adj}(\partial I)}
                \varepsilon_iR_i^n
 \le Cv^{10}S+C\mathcal R_*/v=O(v^7S).
\]
因此原二轴单源界给
\[
 \sum_{\text{边界二轴}}I_e
 \le Cv^2O(v^7S)+Cv^4S(\partial I)=O(v^9S).
\tag{27}
\]
不按非好根权乘 gross 剔除任何端点。全部 \(q\ge3\) gross 已由已审 R050 给
\(o(v^8S)\)。零／一道临界 pair 依原已证强面容量引理无正体积交，故上述分组覆盖完整 \(H\)，没有漏掉一个单临界补集。这两个步骤均不消耗另一份 \(J\)。

## 9. 合并、容量归一及极限次序

先给定 \(\delta>0\)。用固定 \(C_*\) 选有限 \(T\) 足够大，使
\(C_LC_*/T<\delta/4\)；再取固定 \(m<c_N/2\) 足够小，使
\(TmC_*<\delta/4\)。最后固定 \(C_{\rm soft}\)，随后令 \(v\to0\)。由 (24)、(26) 可把小角及长尾上极限共同放进任意预定的 \(\delta v^8S\)，加上 (27)、R050 和 (19)、(23)，即得 (4)。

实际正常集合、正常模型 \(Y_{2,N}\)、\(D_N\) 和 \(\mathsf P_N\) 均与这些新阈值无关。软权常数只改变全软容量误差 \(O_{C_{\rm soft}}(v^3S)\)，乘六阶主价后为 \(o(v^8S)\)。同一 clean helper／碰撞分片可保持不变，正常软权仍为一。

因为 \(\mathcal B_{T,m}\) 的软权为一，
\[
 W_B\le\mathcal C-W_N+O_{C_{\rm soft}}(v^3S)
       =D_N+O_{C_{\rm soft}}(v^3S).
\]
代入 (4)、(5)，有
\[
 \begin{aligned}
 \mathsf P_N+(1-\gamma_\nu)A_6W_B
 &\le q_N\mathcal C+\mathsf E_N
    -[q_N-(1-\gamma_\nu)A_6]D_N+o(v^8S).
 \end{aligned}
\]
由 \(q_N=A_6+O(v^7)\)，括号至少为固定 \(c_Dv^6\)，证明 (6)。这里 \(D_N\) 可能真为 \(O(v^2S)\)，负的领先槽损失被明确保留，而非放入八阶误差。

等价极限口径是：每个辅助 \(T,m,C_{\rm soft}\) 先固定，取 \(v\to0\)；随后让 \(T\to\infty\)，并让 \(Tm\to0\)。这可以通过任意固定 \(\delta\) 的上述选择实现。没有一个阈值随 \(v\) 变化，也不需要新窗代价常数对 \(T\) 一致。只有用于尾部的 \(C_*,C_L\) 和全域非正常折扣在 \(T,m\) 前固定；第 5 节已逐项证明这一点。

## 10. 全局接头的范围

本稿只给同一实际末层的补集结算，并以明确正常模型泛函保留 R053 的全部未审内容。它不把 normal-only 减全 \(J\) 的价格与 (4) 相加，不声称一个模型系数就是 \(C_n\) 的上系数。

对于当前 \(n\ge3\)，原精 dyadic 小源、盒内扩张余项及环境的绝对阶
\(v^{n+6}S\) 已为 \(o(v^8S)\)；盒外的 \(Cv\mathcal F\) 由 (1) 也是九阶，故不再产生一个八阶常数。前缀及全局图仍须在最终总装中保持同族；本稿不替其省略原始源缺量望远镜。

**二维不在本稿的全局升级范围。** 当 \(n=2\)，盒内小源首项可能恰为八阶，不能因局部 source 算子成立就忽略它或否定新的两尺度输入。本稿没有使用 R054 尚在审查的 \(P_I=O(v^8S)\) 接口；上述边界处理只需原 \(S(\partial I)=O(v^5S)\) 和真实缺量分解。

## 11. 计算证据与当前状态

新有限代数只需真实 clip 的增量关系：
\(0\le\min(H,u+w)-\min(H,u)\le w\)，
以及 \(H\ge u+w\) 时其等于 \(w\)，\(H,u,w\ge0\)。
已实际读取 Terra 单核成功的 r055-computation/clip_increment.wl 及 clip_increment.log：两项违反集均为 False。首次路径笔误日志不作证据。不以程序代替第 4 节中分面薄带的物理比例证明。全域非正常强制性复用已审 R045 的解析证明，不重跑相同大消元。

R051/R052 已审整源高度、实际反单调迭代和 R049/R050 的固定竞争族稳定性是本稿的数学依赖。实际读取的核心供货为：

* R049 源 SHA 5c4d7b371a0c33c39449afd686de599c49570aaa639d4f7db7d3c1f166a040f4；
* R050 源 SHA 22f848879b3c9ba839edb5d307a205ae3245e5ba3a89a99f7ec2ce2974010748；
* R051 源 SHA 37bbb0ec63fac1dce33fb9883e14499488baef490341df32a6fdfb0da427c979；
* R052 源 SHA d655c1e160c88914774248c81688c6bdde3f27614387c4eb7f0b9c1a0a0a5296，其正式独审明确供应 (25) 的完整强 CDF 正余量，解除 R051 历史条件；
* R045 的全域 \(c_\nu\)、正常小管长角估计及 R036 同一 normal clean 分配；R053 的具体待审模型系数未在本文使用。

状态：本稿提供 (4)、(6) 的完整补集证明候选，待独立审核。没有给出或预设 R053 正常八阶系数，没有宣布精确八阶 Tauberian 定理。
