# R043　正常角槽空缺与实际高支撑交量的联合预算

日期 2026-11-05。新候选，供独立审查。只处理同一末层的局部几何，不据此宣布全输入八阶上界。原仓库和前轮冻结稿只读。

固定 \(n\ge3,\ a=v^n\)，采用 R040 同一末层的原有限族、实际源 \(F_i\)、目标根 \(R_i\)、中点盒 \(X_i=x_i+W_i+[-R_i/2,R_i/2]^n\)。写 \(S=\sum r_i^n\)。本稿使用已审 R017 的固定裕量面对预算、R029 的完整固定窗口及实际根传输、R035/R036 的实际 adaptive-clean 四份角槽分配；所有量均在本末层重取。具体先验为
\[
 \sum\varepsilon_i r_i^n=O(v^5S),\qquad
 \sum\Delta_i^{\rm end}r_i^{n-1}=O(v^6S),\qquad
 \mathcal F_\tau:=\sum_{\mathscr E_\tau}\Gamma_{cd}r_c^{n-1}=O_\tau(v^6S).
 \tag{1}
\]
这里最后一式对任意预先固定的稳健切向裕量 \(\tau>0\) 成立；\(\Gamma_{cd}\ge0\) 是真正单临界面对的实际中点端隙。另有 R019/R029 坏根 \(Z_v\)，权 \(S(Z_v)=O(v^{9/2}S)\)，其补集的任意固定物理窗口有完整清单、统一可比尺度及全部深处源。相交度、固定窗口清单数和窗口重复次数均有固定上界。无需给 LT 收敛速度。

## 1　实际正常槽缺量

固定足够小的 \(\nu>0\)，取 R037 意义下的正常深处二临界边集 \(\mathscr N\)。本稿几何只用其以下后果
\[
 c_n v^3R_e\le L_{e,1},L_{e,2}\le T_0v^3R_e,\qquad
 b_e\ge(1-\nu)R_e^{n-2},\qquad R_e=\max(R_i,R_j),
 \tag{2}
\]
其中 \(T_0<\infty\) 固定。正常判据中的根差限制可以保留，但下面不依赖它。取 \(\nu\) 使正常两端在各非临界方向的实际中心差小于 \(r_i/100\)。这由 (2)、各 \(X\) 边长为 \(R+w\)、\(w=O(vr)\)、相交尺度比 \(1+O(v)\) 直接得到。

一个实际角槽记为 \(s=(i,q,l,\sigma_q,\sigma_l)\)，其横向容量是
\[
 F_s=\pi_{\widehat{ql}}X_i .
\]
对 \(\mathscr N\) 中两端都不在 \(Z_v\) 的边，作 R036 的 adaptive-clean 四份分配，阈值为
\(\lambda_e=\min(L_{e,1},L_{e,2})/(1000r_{\rm anchor})\)。
令 \(O_s\subset F_s\) 为送入该槽的全部正常 clean 片之并。R036 的单射保证片在每槽几乎处处不重复。定义真正非负的槽空缺
\[
 \boxed{\displaystyle
 D_{\mathscr N}=\frac14\sum_s R_i^2|F_s\setminus O_s|.}
 \tag{3}
\]
占用可在多个 helper 指标之间分片；这里没有假设一条正常边占满一个 helper 的整个槽。

写 \(W_{\mathscr N}=\sum_{e\in\mathscr N}b_eR_e^2\)、
\(\mathrm{Cap}=\sum_sR_i^2|F_s|\)。已审 clean 漏片收费和实际根传输给
\[
 \boxed{\displaystyle
 D_{\mathscr N}=\frac14\mathrm{Cap}-W_{\mathscr N}+O(v^3S).}
 \tag{4}
\]
确实，正常边软权为一，全部漏片的 \(R_e^2\) 权是 \(O(v^3S)\)；坏端点边权 \(O(v^{9/2}S)\) 更小。四份分配中 \(R_e^2\) 换为分配胞的 \(R_c^2\)，由 R029 固定窗口根传输总计仅 \(O(v^6S)\)。对分片先按同一 helper 取并后求面积，不能将误差再乘分片数。

固定
\[
 0<\delta<\frac1{100}\left(\frac13\right)^{n-2}.
 \tag{5}
\]
若某外部胞有至少一个槽满足 \(|F_s\setminus O_s|\ge\delta r_i^{n-2}\)，称其为缺槽胞。由 \(R_i\asymp r_i\)，
\[
 \sum_{\text{缺槽胞}}r_i^n\le C_{n,\delta}D_{\mathscr N}.
 \tag{6}
\]
再固定足够大的物理窗口 \(R_*=100(n+1)\)。称根 \(i\) 坏，若 \(i\in Z_v\)，或 \(i\notin Z_v\) 而其完整清单 \(\mathcal I_i(R_*)\) 中有缺槽胞。记坏根集为 \(\mathcal B\)。窗口尺度比较和 packing 给
\[
 \boxed{\displaystyle S(\mathcal B)\le C_{n,\delta,R_*}
          \big(D_{\mathscr N}+v^{9/2}S\big).}
 \tag{7}
\]
对 \(i\notin Z_v\) 的清单，所有相关源都深处；因此无需把一个未量化 LT 误差作为额外坏权。每个缺槽胞被可比尺度的固定窗口计入至多固定次，证明 (7)。

## 2　正常占用迫使真正面邻居的切向端点对齐

以下是新增的实际几何步骤。先按 R014 在某深处胞 \(A\) 的正 \(q\) 方向选一个真正面邻居 \(B\)。它满足
\[
 x_{B,q}-x_{A,q}=r_A+O(vr_A),\qquad
 |x_{B,h}-x_{A,h}|\le(1/2+1/128)r_A\quad(h\ne q),
 \tag{8}
\]
以及 \(r_B/r_A=1+O(v)\)。假设 \(A\) 不是缺槽胞，所需清单在一个好根固定窗口内。

固定 \(l\ne q\)。考察 \(A\) 的 \((q+,l+)\) 槽。对其余 \(n-2\) 轴，\(X_A,X_B\) 的投影交各至少为 \(r_A/3\)，故公共横向块
\[
 T_{AB}^{ql}=\pi_{\widehat{ql}}X_A\cap\pi_{\widehat{ql}}X_B
\]
的体积至少为 \(3^{-(n-2)}r_A^{n-2}\)。由 (5)，可选一个避开全部零测边界的 \(y\in T_{AB}^{ql}\cap O_s\)。该点属于一条实际正常 clean 边送入 \(A\) 的片。

按 R035 的实际四胞布局，把这条片的四个胞重新记成
\[
 \begin{matrix} C&D\\ A&\widetilde B\end{matrix}.
 \tag{9}
\]
因为共享槽的符号已固定，\(A\) 在此布局中必为左下。正常边为 \(AD\) 或 \(\widetilde BC\)。四胞平面中心相对于 \(A\) 的近单位列误差不超过 R035 固定 probe 常数的一个固定倍数，例如 \(100\epsilon\)，\(\epsilon=10^{-6}\)。四胞都含同一个实际横向点 \(y\)。

**识别 \(\widetilde B=B\)。** 若二者不同，则由 (8) 和 (9)，它们在 \(q,l\) 平面的交面积至少为 \(r_A^2/4\)，而它们的横向投影都含 \(y\)。两胞均在该正常边锚点的完整固定清单内：经 \(A\) 比较，平面距离有固定上界，横向距离因共同 \(y\) 亦有固定上界，尺度可比。于是它们触发该边自己的 collision 判据，因为 \(\lambda_e=O(v^3)<1/100\)。这与该分配片 clean 矛盾。故确为同一个实际胞。这个论证不把缺失一个指定 helper 当作 void，也不要求各槽选择同一条边。

现在用实际端点写 \(d=B_l^+-A_l^+\)，上下标 \(\pm\) 均指 \(X\) 的实际端点。

若正常对角为 \(AD\)，则 \(L_{AD,l}=A_l^+-D_l^->0\)。正常两端 \(A,D\) 在所有横向轴的中心差 \(<r_A/100\)，而 \(A,B\) 满足 (8)。因此 \(B,D\) 只有 \(l\) 方向临界，其他方向都有固定切向裕量；在 \(q\) 方向还由 (9) 有固定大重叠。它是 R017 的稳健面对，故
\[
 \Gamma_{BD}=D_l^- -B_l^+\ge0,\qquad
 d=-L_{AD,l}-\Gamma_{BD}.
 \tag{10}
\]
若正常对角为 \(BC\)，则同理 \(A,C\) 是稳健面对，且
\[
 \Gamma_{AC}=C_l^- -A_l^+\ge0,\qquad
 d=L_{BC,l}+\Gamma_{AC}.
 \tag{11}
\]
这些辅助面对的原半盒法向缝非负，因其他方向无法实现原强中心分离；由 (9) 小于 \(r_A/8\)。其切向交长至少为固定倍数的 \(r_A\)，可统一取 \(\tau=1/16\)。故 (1) 对它们适用，而不仅是对最初选定的面邻居适用。

反射 \(l\) 轴并使用 \((q+,l-)\) 槽，得到左端同样的结论。所有法向符号反射亦然。于是，对每个真正选定面对 \(A\to B\) 和每个切向轴 \(l\)，
\[
 \boxed{\displaystyle
 |B_l^\pm-A_l^\pm|
 \le C_nT_0v^3r_A+
       \text{一个固定窗口内的稳健实际面对端隙}.}
 \tag{12}
\]
这里没有以 \(vR\) 代换任何源宽；左右端点分别处理，宽度差自动包含在同一精确恒等式中。

## 3　有限面对路径识别真实高支撑目标

对好根 \(i\notin\mathcal B\)，在其 \(R_*\) 窗口内定义长度量
\[
 \mathcal E_i=
 \sum_{\substack{c,d\in\mathcal I_i(R_*)\\
                 \{c,d\}\text{ 稳健面对},\,\tau=1/16}}\Gamma_{cd}
 +\sum_{c\in\mathcal I_i(R_*)}\Delta_c^{\rm end}.
 \tag{13}
\]
无向面对仅计一次。由 (1) 和固定窗口重复次数，
\[
 \sum_{i\notin\mathcal B}\mathcal E_i r_i^{n-1}=O(v^6S).
 \tag{14}
\]

取任意深处高支撑相交对 \(\{i,j\}\)，令所有者满足 \(R_i\ge R_j\)，且先设 \(i\notin\mathcal B\)。反射所有临界轴使 \(j\) 均在正侧，临界轴数为 \(q\ge3\)，交长 \(0<L_l\le C_Avr_i\)。

若 \(\mathcal E_i\ge c_*vr_i\)，直接有
\[
 |X_i\cap X_j|\le C_Av^qr_i^n
 \le C_{A,c_*}v^{q-1}\mathcal E_i r_i^{n-1}.
 \tag{15}
\]
以下设 \(\mathcal E_i<c_*vr_i\)，其中 \(c_*>0\) 在维数及原选择常数之后足够小地固定。

从 \(i\) 出发，依次沿每个临界轴正方向走一次 R014 的选定真正面对邻居，得 \(q\) 步路径，终点记 \(k\)。其长度至多 \(n\)，所有节点及 (12) 所需的辅助四胞均在 \(R_*\) 窗口内。每个节点不缺槽，所以 (12) 可逐步使用。每步根差由 R017 精确式
\[
 |R_c-R_d|\le2\Gamma_{cd}+\Delta_c^{\rm end}+\Delta_d^{\rm end}
\]
支付，故 \(|R_k-R_i|\le C_n\mathcal E_i\)。

在某个临界轴 \(l\)，走到该轴之前的各步均为切向步，保持右端点到 \(O_n(T_0v^3r_i+\mathcal E_i)\)；该轴一步把新左端放在旧右端加非负实际面隙的位置；其余各步保持这个左端到同样误差。因此
\[
 |X_{k,l}^- -X_{i,l}^+|
       \le C_n(T_0v^3r_i+\mathcal E_i).
 \tag{16}
\]
在非临界轴，路径每步都是切向步，故左右两端均与 \(X_i\) 的对应端相差同样的量。

**识别 \(k=j\)。** 在临界轴，\(X_{j,l}^-=X_{i,l}^+-L_l\)，故 \(x_k,x_j\) 在这些轴的差为 \(O_{n,A}(vr_i)+C_n\mathcal E_i\)。在非临界轴，使用原临界判据 \(A_1=5(1+A)\) 以及所有者 \(R_i\ge R_j\)，有
\[
 |x_{j,l}-x_{i,l}|
 <(1/2-A_1v)(r_i+r_j)
 \le \frac{1-2A_1v}{1-Av}R_i.
 \tag{17}
\]
这里用原 \(R_h\ge(1-Av)r_h\)。实际 \(X\) 中心与 \(x_h\) 的差每轴至多 \(Avr_h\)。由非临界轴的左右端传输，
\[
 |x_{k,l}-x_{i,l}|
 \le (2A+o(1))vR_i+C_n(T_0v^3r_i+\mathcal E_i).
 \tag{18}
\]
而 \(r_k\ge R_k\ge R_i-C_n\mathcal E_i\)。先取 \(c_*\) 小、再取 \(v\) 小，(17)、(18) 的严格剩余裕量给
\(|x_{k,l}-x_{j,l}|<r_k\)；临界轴的 \(O(vr_i)\) 差亦小于 \(r_k\)。若 \(k\ne j\)，遂违反原 \(\|x_k-x_j\|_\infty>\max(r_k,r_j)\)。故 \(k=j\)。

于是 (16) 给所有临界交长同时满足
\[
 L_l\le C_n(T_0v^3r_i+\mathcal E_i).
 \tag{19}
\]
结合原 \(L_l\le C_Avr_i\)，按 \(\mathcal E_i/r_i\le v^3\) 或 \(>v^3\) 分支，得到
\[
 \boxed{\displaystyle
 |X_i\cap X_j|
 \le C_{n,A,T_0}\big(v^{3q}r_i^n+
                         v^{q-1}\mathcal E_i r_i^{n-1}\big).}
 \tag{20}
\]
非临界横向交长只用了原 \(O(r_i)\) 上界。由有限度及 (14)，
\[
 \boxed{\displaystyle
 \sum_{\substack{q(e)\ge3\\\text{owner}(e)\notin\mathcal B}}I_e
        =O(v^8S).}
 \tag{21}
\]
尤其这里没有把三临界片的薄横宽当成槽空缺的权；先用整胞各槽的实际覆盖识别端点关系，再估计交量。

## 4　坏根净价接口及联合结论

本节接 R044 的独立真实高支撑净价候选引理，冻结源为
“r044-high-support-net/high_support_net.md”，SHA256
ef3a6df6ee4581e990fd55fbe76c2b03f43b2e5bdbd0cb27f7eba18231cf0cc4。
已全文读取其 (27)、(28) 及假设，以下只作这条独立接口的合成；其新证明仍须独审。
它使用同一个源、目标、指定半径以及同一个固定二轴管
\[
 \mathcal C_i(K)=
 \{x\in x_i+W_i':\text{至少两个轴到有效盒某面的距离}\le Kv^3r_i\},
\quad
 Z(K)=\frac va\sum_i|[(x_i+W_i')\setminus F_i]\setminus\mathcal C_i(K)|.
 \tag{22}
\]
固定 \(K\) 后取 \(T(K)\) 足够大。R044 (27)、(28) 给，对任意固定 \(\eta>0\)，
\[
 \sum_{\substack{q(e)=3\\\text{owner}(e)\in\mathcal B}} I_e
 \le \eta Z(K)+C_{\eta,K}v^7S(\mathcal B)+C_{\eta,K}v^8S,
 \tag{23}
\]
且全部 \(q\ge4\) 亦有 \(\eta Z(K)+C_{\eta,K}v^8S\) 的界。该引理不使用正常角价格；管外实际矩形与所有管内正常缺层物理分离。其长尾实际净价更强，为 \(C_\eta v^{15/2}S(\mathcal B)\)；(23) 的 \(v^7\) 来自 \(p_2\le T(K)v^3\) 支的 gross 交量。

由 (7)、(21)、(23)，对任意预先固定的 \(\eta>0\)，重新选 R044 中的系数并令 \(v\) 足够小，得到
\[
 \boxed{\displaystyle
 \sum_{q(e)\ge3}I_e
 \le \eta Z(K)+\eta v^6D_{\mathscr N}
                      +C_{\eta,K,\nu}v^8S.}
 \tag{24}
\]
边界 \(q\ge3\) 交量由原单源界、\(S(\partial I)=O(v^5S)\) 及 (1) 为 \(O(v^8S)\)，并入末项。全部源基线保持。没有从完整 \(J_i\) 抽一个会污染六阶价格的份额。

若正常 packet 的实际反馈禁层 \(Y_i^{\rm fb}\subset\mathcal C_i(K)\)，则
\[
 Z(K)\le\frac va\sum_i|\mathscr D_i\setminus Y_i^{\rm fb}|,
 \qquad \mathscr D_i=(x_i+W_i^0)\setminus F_i .
 \tag{25}
\]
因此 (24) 与 R042 保留的 raw CDF residual 是同一物理预算方向。对任意给定的正 residual 系数，可先把 \(\eta\) 取小。

还须在全二轴结算中真正保留 \(v^6D_{\mathscr N}\) 的正系数；不能凭正常 packet 一项就假定它已存在。式 (4) 说明所需量精确等于正常边未使用的实际四份容量，误差乘 \(v^6\) 仅 \(O(v^9S)\)。若另有非正常二轴的定量标量正罚，可将其与总软容量的空缺联合为该项。此为 R045 二轴整体结算的明确接口，不在本稿冒充已完成的全输入八阶上界。

## 5　范围与本机算术

第 1–3 节是本稿新增的实际槽几何，准确结论为 (3)–(7)、(12) 及 (20)–(21)，只依赖已有已审输入。第 4 节的 (24) 另依赖上述确切 SHA 的 R044 新净价候选，不把其尚未完成的独审当作已完成。审查尤其应核 (9) 的同 fibre helper 识别、(10)–(11) 辅助面对的完整切向裕量、(17)–(18) 原中心严格分离，以及所有窗口重复次数。

Terra 已在本机实际运行 r043-computation/branches.wl 和 branches.log。好窗乘积界的两个分支、AD/BC 端点绝对值界、三步路径三角界的违反集均为 False，望远镜残差为零。一般 \(q\) 长尾的幂和 Young 分支归入 R044 的 r044-computation/young_branches.wl/log/json 与 Fraction 证书，不在此重复计算。软件只核所列代数，不代替实际清单、可测分配和源缺集的几何证明。
