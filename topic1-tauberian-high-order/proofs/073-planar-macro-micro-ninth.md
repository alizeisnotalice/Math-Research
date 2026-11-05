# 二维宏微两尺度的九阶实际下界：完整 host 费用与有限恢复

2026-11-05。固定维数 n=2。本文只给一个明确实际输入族的九阶展开及相应下界；不证明九阶归一化极限存在，也不声称该族在所有输入、微孔分布或宏观参数中最优。旧冻结源、论文和其他代理文件未改。

## 1. 命题、参数及常数次序

令
\[
 u=\sqrt t,\quad v=2u,\quad a=4u^2,\quad
 g=(1+2u)^{-1},\qquad L_2=\frac1{1-ag^2},
\]
\[
 P_7=L_2+\frac4{27}av^6-\frac{176}{243}av^7,\qquad
 \widehat J_2=\frac{17638}{6561}.
 \tag{1}
\]
这里 P₇ 保留完整 L₂，不是把 L₂ 自身截到七阶的多项式。

**实际有限构造命题。** 存在 u₀>0，使每个实 0<u<u₀ 都有有限测度输入 Eᵤ，满足
\[
 \frac{|\mathcal H_2(1-u^2;E_u)|}{|E_u|}
 =P_7+\widehat J_2av^8
       -\frac{23295344}{2657205}av^9+O(av^{19/2}).
 \tag{2}
\]
因此同式去掉等号并把最后项改为 −Cav¹⁹ᐟ² 是 C₂ 的实际下界。余项为 o(av⁹)。本文不从单侧 (2) 推出 C₂ 的九阶系数。

构造取已明确给定的解析剪切
\[
 \kappa(u)=\frac{16}{3}g^2+\frac{640}{81}g^3u,\qquad
 \kappa_0=\frac{16}{3},\qquad e=\kappa u^3,
\]
\[
 T=\begin{pmatrix}1&e\\e&1\end{pmatrix},\quad d=1-e^2,
 \quad w=ug+\Lambda u^6,\quad
 \rho=1-2w=g-2\Lambda u^6,\quad M=a\rho^2.
 \tag{3}
\]
M 从此是宏孔加微孔的总设计质量，不由某次删层后的剩余宏孔质量重算。

先定义显式纯宏参数
\[
 \Lambda_{\rm p}
 =\frac{\kappa^3}{32g}
 -\frac{5\kappa^4}{768g^2}u
 +\frac{13\kappa^5}{8192g^3}u^2
 -\frac{33\kappa^6}{81920g^4}u^3,
 \qquad b=4u^2-16u^3.
 \tag{4}
\]
本次实际参数为
\[
 \boxed{\Lambda=\Lambda_{\rm p}-b\Lambda_{\rm p}
                +\frac{5\kappa_0^4}{192}u^3+u^{7/2}.}
 \tag{5}
\]
u⁷ᐟ² 是主动留下的正安全偏移；它对比值只造成 O(av¹⁹ᐟ²) 损失，并非对自然下一阶指数的判断。

所有几何和积分常数先在包含 (3)–(5) 的固定正紧 κ 区间、固定有界 Λ 区间上取得，再缩小 u₀。可同时要求
\[
 e<w/4,\quad \rho>2w,\quad
 \delta:=\frac e{2a\rho}<1/8,\quad
 F^{-1}(2e^2)<w-e,\qquad
 F(s)=a(2\rho s+s^2).
 \tag{6}
\]
因为 e=O(u³)、w∼u、δ∼u，这些守卫一致成立。后文所有 O 常数只依上述固定区间。每次使用的 w、ρ、F、K、H 都按 (5) 的新参数重新定义。

## 2. 完整物理算子与有限宏观 host

写 Cw=[−w,w]²，s=w−‖x‖∞，取长规则半径 r(x)=ρ+s。对任意可测 V⊂Cw 定义
\[
 R_V(x)=\sum_{k\in\mathbb Z^2\setminus\{0\}}
 |(Tk+(C_w\setminus V))\cap Q(x,\rho+s)|,
 \qquad
 \mathcal B(V)=\partial C_w\cup\{R_V(x)\ge F(s)\}.
 \tag{7}
\]
它是实际外国孔交量，不是只扣一个名义帽的松弛。r≥ρ>2w 包含同周期全部实际孔。

若正交量来自 Tk+Cw，则
‖Tk‖∞≤‖x‖∞+r+w=1。由 ‖Tk‖∞≥(1−e)‖k‖∞，模至少二的指标不能进入。剩余指标中轴邻孔仅接触，同号双支撑不进入。只有 k=(1,−1)、(−1,1) 可能产生正面积，每个帽面积≤e²。因此对所有 V，坏点均在
\[
 s\le h_H:=F^{-1}(2e^2)\le2e\delta<w-e.
 \tag{8}
\]
深部可同时见两目标，但 (8) 已先将其排除。

浅层只有最近面符号 (+,−)、(−,+) 两个活动角，每角一个目标。其两面深为 (s,s+y) 或交换，0≤y<e，令 q=1−y/e。目标帽从目标面对端计为
[0,e]×[0,eq] 或交换，独立于 s；四张图的 Jacobian 都是一。输出高度≤eδq，而原纤维顶端 w−y 满足
w−y−eδq≥w−e>0，故此处真实 clip 不活跃。对后加任意 H，只要目标仍是 Cw 的子集，同一守卫仍成立。

帽质量随矩形端点连续，故 B(V) 为闭集。它反单调。取固定有限次数
\[
 U_0=\mathcal B(\varnothing),\quad U_{j+1}=\mathcal B(U_j),
 \qquad K=U_4.
 \tag{9}
\]
有 B(K)=U₅⊂K；这就是本次全部微孔必须保留的 host。没有调用未知无限固定点生成实际孔。

为记录必要精度，令
\[
 \Psi(z)=\frac{2\delta z}{1+\sqrt{1+4a\delta^2z}},\quad
 \Phi(h)=h/\delta+ah^2,\quad \Psi=\Phi^{-1}.
 \tag{10}
\]
全局中心反射和交换两轴保持本二维晶格；不使用分别翻转单轴的对称。故 Uⱼ 的两个活动角有相同高度 hⱼ(q)，
\[
 h_0(q)=\Psi(q),\qquad h_{j+1}(q)=
 \Psi(q-\mathcal M_{h_j}(q)),\qquad |U_j|=4e^2\int_0^1h_j(q)\,dq,
 \tag{11}
\]
\[
 \mathcal M_h(q)=
 \int_0^q\min\{h(1-z),q-z\}\,dz+
 \int_0^1\min\{h(1-z),q\}\,dz.
 \tag{12}
\]
这两个 min 包括近零帽。hⱼ(0)=0、0≤hⱼ′≤δ，所以 hⱼ(q)≤δq；由帽的截线长度小于一，hⱼ 在内部严格递增。

由 min 的 Lipschitz 性和 Ψ′≤δ，
‖hⱼ₊₁−hⱼ‖∞≤2δ‖hⱼ−hⱼ₋₁‖∞。首差 O(δ²)，从而
\[
 |K\setminus\mathcal B(K)|=O(e^2\delta^6)=O(u^{12}).
 \tag{13}
\]
同时固定目标点打开 Tonelli，四图中可见底变量 y 的总长≤4e，给完整集合算子的
\[
 |\mathcal B(V)\triangle\mathcal B(W)|
 \le4\delta|V\triangle W|.
 \tag{14}
\]
这里允许任意可测 V,W、任意薄帽及跨角目标，不要求其为对称高度图。

## 3. 所需宏观积分及到九阶的显式参数

本节只使用已有真实有限层积分的四项，不重做常容量全阶递归。为使新预算可从定义核验，先写出有限积分规格。

把 Ψ 改为 δz 的线性层记为 hⱼˡ。因
Ψ(z)=δz+O(aδ³) 于 0≤z≤1 一致成立，(11) 的有限收缩归纳给
‖hⱼ−hⱼˡ‖∞≤Caδ³（j≤5）。其对真实面积的误差 O(ae²δ³)=O(u¹¹)，已晚于本次所需面积 u¹⁰。

在线性层中，置
\[
 M_0(q)=
 \begin{cases}
 q-(1-\delta)q^2/(2\delta),&0\le q\le\delta,\\
 \delta(1-2\delta+2q-q^2)/(2(1-\delta)),&\delta\le q\le1,
 \end{cases}
 \quad h_1^{\rm l}(r)=\delta[r-M_0(r)].
 \tag{15}
\]
对满足 0≤h(r)≤r 的高度，Tonelli 的精确恒等式为
\[
 \int_0^1\mathcal M_h(q)\,dq
 =\int_0^1[(1+r)h(r)-h(r)^2]\,dr.
 \tag{16}
\]
设 P(r)=∫₀ʳh₁ˡ、A=P(1)，则 M(h₁ˡ) 的分支为
\[
 q=h_1^{\rm l}(r):\quad m_L=q^2/2+q(1-r)+P(r),
 \qquad
 q=1-r+h_1^{\rm l}(r):\quad m_H=2A-P(r)+h_1^{\rm l}(r)^2/2.
 \tag{17}
\]
像区间分别为 [0,h₁ˡ(1)]、[h₁ˡ(1),1]，没有两个独立反根同时出现。
令 Zδ(q,m)=δ(1+q)(q−m)−δ²(q−m)²，则
\[
 I_2:=2\int_0^1\mathcal M_{h_2^{\rm l}}(q)\,dq
 =2\int_0^1 Z_\delta(h_1^{\rm l},m_L)(h_1^{\rm l})'\,dr
 +2\int_0^1 Z_\delta(1-r+h_1^{\rm l},m_H)
                   [1-(h_1^{\rm l})']\,dr.
 \tag{18}
\]
只需把两个积分在 r=δ 处分段；(15)、P 都是分段多项式，(18) 为有限有理积分。逐项积分给
I₂=5δ/3−13δ²/4+33δ³/5+O(δ⁴)。
由 ∫2h₃ˡ=δ−δI₂，再用 h₄ˡ−h₃ˡ=O(δ⁵)，得到
\[
 |K|=2e^2\left[\delta-\frac53\delta^2+
                  \frac{13}{4}\delta^3-\frac{33}{5}\delta^4\right]
          +O(u^{11}).
 \tag{19}
\]
(18) 是四个定限多项式积分的核验式，所用系数的原始精确计算保存在已冻结的平面积分材料；它们不承担新 host 几何证明。

由 δ=κu/(8ρ)、ρ=g+O(u⁶)，(19) 等价于
\[
 \frac{|K|}{a}
 =\frac{\kappa^3}{16g}u^5-\frac{5\kappa^4}{384g^2}u^6+
   \frac{13\kappa^5}{4096g^3}u^7
   -\frac{33\kappa^6}{40960g^4}u^8+O(u^9).
 \tag{20}
\]
在本次全部有界 Λ 上统一。特别 (4) 给
\[
 |K|=8u^7\Lambda_{\rm p}+O(u^{11}),\qquad
 D:=(2w)^2-M
   =8\Lambda u^7+4\Lambda^2u^{12}(1-4u^2).
 \tag{21}
\]
(20) 是按新 w 的真 K 计算的；把 ρ 换成 g 的误差已晚于所写余项，不是把旧孔宽的层直接搬来。

还将用到较粗且与后面新增费用有关的式子
\[
 h_4(q)=\delta q+O(\delta^2),\qquad
 \int_0^1\mathcal M_{h_4}(q)\,dq
       =\frac56\delta+O(\delta^2).
 \tag{22}
\]
第一式来自 (11) 的 M_h≤2δ；第二式直接由 (16)，因为 ∫₀¹(1+r)r=5/6。近零帽没有从积分中删除。

## 4. 密集真微孔：全部矩形帽的一致误差

必须区别 K 的零测代表和有面积的两角。K 按 (7) 人为含整个 ∂Cw，其字面边界长度可为 O(w)。令
Kreg=closure(int K)，则 |K△Kreg|=0；Kreg 正是两个活动角的四张连续图。每张边界由参数 q∈[0,1] 的 e-Lipschitz 曲线及有限线段组成，故长度≤Ce。
于是对 0<h≪e，
\[
 |\{x:\operatorname{dist}_\infty(x,\partial K_{\rm reg})\le 2h\}|
 \le Ceh.
 \tag{23}
\]
可把每条参数曲线按步长 h/e 分成 O(e/h) 段，每段的 h 邻域落在边长 Ch 的方盒内，直接证明 (23)。不以整个 ∂Cw 代替该有面积边界。

取
\[
 t'=u^2(1-u^2)<t,\quad v'=2\sqrt{t'},\quad a'=4t',
 \quad \beta=\frac{a'}{(1+v')^2}
       =4u^2-16u^3+O(u^4)=b+O(u^4).
 \tag{24}
\]
令 h=2⁻ᵐ，其中 m 为满足 h≤u²⁰ 的最小非负整数。只保留 q∈hℤ² 中
Q(q,h)⊂int Kreg 的中心，基本胞为 Cq=q+[−h/2,h/2)²。置
\[
 w_\mu=\frac{hv'}{2(1+v')},\qquad
 H=\bigcup_q(q+[-w_\mu,w_\mu]^2),\qquad S_h=\bigcup_q C_q.
 \tag{25}
\]
这是有限个不交紧微孔之并；2wμ<h，且 H⊂Sₕ⊂int K⊂K。遗漏的 Kreg\Sₕ 在 (23) 的边界邻域内：若点距边界>2h，它所属胞中心的 Q(q,h) 必在 int Kreg。因此
|Kreg\Sₕ|≤Ceh。

每个完整基本胞中微孔质量精确为 βh²。任一实际目标矩形 Π 的两边≤e，其边界所遇基本胞数≤C(e/h+1)。在这些边界胞上，微孔质量及 β 倍的胞面积均≤βh²。其余完整胞精确抵消。再付 Kreg\Sₕ 的 β 倍面积，得到
\[
 \boxed{\quad
 |\,|H\cap\Pi|-\beta|K\cap\Pi|\,|\le C\beta eh,\qquad
 |H|=\beta|K|+O(\beta eh).\quad}
 \tag{26}
\]
第一个式子对四张图的所有 q∈[0,1] 同时成立，包括宽度趋零的帽；不需要正最小帽宽。这里是由细格全胞误差得来的真实分布结论，不是把 H 按平均密度未经证明地代入每个帽。

对任意 z 属于以 q 为中心的微孔，保留短半径
\[
 r_\mu(z)=h-w_\mu-\|z-q\|_\infty.
 \tag{27}
\]
其见证在 Q(q,h) 内，包含整个自身微孔；与任意其他同格微孔至多边界接触。因为 rμ≥h−2wμ，
\[
 \frac{|H\cap Q(z,r_\mu)|}{|Q(z,r_\mu)|}
 \le\left(\frac{w_\mu}{h-2w_\mu}\right)^2=t'<t.
 \tag{28}
\]
这是所有微孔点（含边界点）的严格原阈值见证。不同宏胞的 Cw 分离，所以该 host 内没有外国微孔或宏孔。

## 5. 全部 host 强制保留后的实际新增保护费

固定已经构造的 H，不再改变它。只作一次实际集合运算
\[
 W=K\cup\mathcal B(K\setminus H).
 \tag{29}
\]
因 W\H⊃K\H，反单调性给
\[
 \mathcal B(W\setminus H)\subset\mathcal B(K\setminus H)\subset W.
 \tag{30}
\]
W 含全部 K，因此以后切出的宏孔不会侵入任何短见证 host。
由于 K\H⊂Cw，B(K\H)⊂U₀，故本节新坏点仍受 (8) 和真实 clip 守卫控制；没有遗漏因额外微孔产生的深源点。

记 E=|W\K|。在四张图上，B(K) 高度为
s₅=F⁻¹(RK)，K 高度为 s₄≥s₅，新增外国帽质量为 mH(y)=|H∩Π(y)|。故新增高度为
[F⁻¹(RK+mH)−s₄]₊。对 A≥s₅ 以及 s₄≥s₅，
0≤(A−s₅)−(A−s₄)₊≤s₄−s₅。积分后由 (13) 得
\[
 E=\sum_{\text{四图}}\int_0^e
       [F^{-1}(R_K+m_H)-F^{-1}(R_K)]\,dy+O(u^{12}).
 \tag{31}
\]
这一步显式处理有限偶层 K 非精确固定层造成的正部；不能把 B(K)=K 当成等式使用。

在 s≤hH=O(u⁴) 上，F′(s)/f=1+s/ρ=1+O(u⁴)，f=2aρ。
所以 (31) 的首项为 (1+O(u⁴))/f 乘四图 ∫mH。
(26) 的统一误差在该积分里至多 Cβe²h/f，远小于 u¹¹。K 对称而 H 不必对称，因此只有将 H 以 (26) 比较之后才能把四图写成同一个 M。由 (22)，
\[
 \begin{aligned}
 E
 &=\frac{4\beta e^3}{f}
           \int_0^1\mathcal M_{h_4}(q)\,dq
       +O(u^{12})+O(\beta e^2h/f)+O(u^{14})\\
 &=\frac{10}{3}\beta e^2\delta^2+O(u^{11})
   =\frac{5\kappa_0^4}{24}u^{10}+O(u^{11}).
 \end{aligned}
 \tag{32}
\]
面积 E∼u¹⁰，经过质量根与真实 halo 比值才首次付到 av⁹。既有 host 面积中的首反馈负项和 (32) 是两笔不同费用：由 (19)，
\[
 |H|-E
 =\beta\left[2e^2\delta-\frac{20}{3}e^2\delta^2\right]
       +O(u^{11}).
 \tag{33}
\]
其中第一笔 −(10/3)βe²δ² 来自原 K 已经变薄，第二笔来自新增保护。不能只扣一次。

若继续使用保留 host 的完整反单调算子
Tₖ,ₕ(V)=K∪B(V\H)，其收缩率为 4δ。它在 V⊃K 上有唯一固定类 Vₖ,ₕ，且 K⊂Vₖ,ₕ⊂W。由 (14)，
\[
 |W\setminus V_{K,H}|
 \le4\delta|V_{K,H}\setminus K|
 \le4\delta E=O(u^{11}).
 \tag{34}
\]
因此对本次给定细格微孔，全 host 后续反馈影响比值只从 av¹⁰ 起；一次保护已给准确九阶费用。此为固定 H、固定 K 的类内结论，不优化其他 host。

## 6. 新孔宽下的正预算与全部严格见证

由 (21)、(24)、(26) 和 h≤u²⁰，
\[
 |H|=8b\Lambda_{\rm p}u^7+O(u^{11}),\qquad
 |K|=8\Lambda_{\rm p}u^7+O(u^{11}).
 \tag{35}
\]
这些估计直接对新 Λ 的 K 成立；不需要借用旧参数的微孔或把旧孔平移。结合 (5)、(21)、(32)，
\[
 \begin{aligned}
 D+|H|-|W|
 &=8u^7(\Lambda-\Lambda_{\rm p})
   +8b\Lambda_{\rm p}u^7
   -\frac{5\kappa_0^4}{24}u^{10}+O(u^{11})\\
 &=8u^{21/2}+O(u^{11})>0.
 \end{aligned}
 \tag{36}
\]
安全偏移的常数固定为一；由于 21/2<11，先取统一余项界再小 u 即得严格正预算，不存在安全常数循环。

W 闭且含 ∂Cw，Cw\W 是开集。由 (36)，其面积严格大于 M−|H|>0；后一个正性来自 |H|=O(u⁹)、M∼4u²。先取足量有限内部紧盒，再对最后一个盒连续切片，可选紧宏孔 A⊂Cw\W，|A|=M−|H|。令
\[
 P=A\cup H,\qquad |P|=M.
 \tag{37}
\]
它是实际总孔，而不是在旧宏孔质量 M 上额外加 H。A 与 H 不交。

对 x∈A，r=ρ+s 包含本胞整个 P。外国实际孔由
(Cw\W)∪H=Cw\(W\H) 上包。由 x∉W 和 (30)，
\[
 M+\sum_{k\ne0}|(Tk+P)\cap Q(x,r)|
 \le M+R_{W\setminus H}(x)
 <M+F(s)=ar^2.
 \tag{38}
\]
严格不等式来自坏集补集，适用于选定紧宏孔的每一个点，不只几乎处处。

对 z∈H，其短见证 (27) 全在 int K⊂W 中，故不遇宏孔 A；其余微孔至多接触。(28) 因而仍是总实际孔 P 的严格见证。由于外包宏孔的周期平移分离，没有外国周期的污染。至此宏、微全部孔点都已有原阈值的实际严格见证。

## 7. 有限外盒、全部半径与精确九阶分数

取整数 N=ceil(u⁻¹²)≥4，ΩN=T([−N,N]²)，定义
\[
 E_u=\Omega_N\setminus\bigcup_{k\in\mathbb Z^2}(Tk+P).
 \tag{39}
\]
保留所有边界切孔。经 T⁻¹ 变换是两坐标 1 周期函数在整数边长区间上的积分，故
|Eᵤ|=(2N)²(d−M) 精确成立；没有把只挖内层孔的输入误写成这个分母。
‖T⁻¹‖∞<2、长半径<1，距周期外盒边界至少三层的全部孔点保留前述见证。边界丢失孔质量 O(NM)，所得比值误差 O(a/N)=O(au¹²)。

对 ΩN 外任意中心及任意正半径，分离半空间与中心方体的对称性给输入占比≤1/2，故不在阈值 1−u²>1/2 的严格 halo 内。输入密度点几乎处处进入 halo。因此完整极大算子（全部正半径）的真实比值满足
\[
 \frac d{d-M}-C\frac aN
 \le\frac{|\mathcal H_2(1-u^2;E_u)|}{|E_u|}
 \le\frac d{d-M}.
 \tag{40}
\]
不是只估计规定半径子族的大小。

对 (3)–(5) 的实际质量，有安全偏移造成 O(au¹⁹ᐟ²)，其余有限有理式可精确展开。Λp 与不含安全偏移的 Λμ 的前三阶 jet 分别为
\[
 \begin{aligned}
 \Lambda_{\rm p}
 &=\frac{128}{27}-\frac{2560}{81}u+
       \frac{24704}{243}u^2-\frac{287952128}{2657205}u^3+O(u^4),\\
 \Lambda_\mu
 &=\frac{128}{27}-\frac{2560}{81}u+
       \frac{20096}{243}u^2+\frac{305512192}{2657205}u^3+O(u^4).
 \end{aligned}
 \tag{41}
\]
代入真实 d/(d−a[g−2Λu⁶]²)，并保留完整 L₂，得
\[
 \begin{aligned}
 Q_{\rm p}
 &=P_7+\frac{15694}{6561}av^8
       -\frac{17871584}{2657205}av^9+O(av^{10}),\\
 Q_\mu
 &=P_7+\frac{17638}{6561}av^8
       -\frac{23295344}{2657205}av^9+O(av^{10}).
 \end{aligned}
 \tag{42}
\]
安全偏移把第二式的余项改为 O(av¹⁹ᐟ²)，再用 (40) 即得 (2)。

同剪切的纯宏 canonical 分支与 Qp 的差为 O(av¹⁰)：其预算根与 Λp 相差 O(u⁴)，质量导数 ∂ΛM=−4aρu⁶。该陈述也可直接由 (19)–(21) 与原正斜率预算根证明。因此两尺度机制相对该纯宏实际分支的改变量为
\[
 \frac8{27}av^8-\frac{496}{243}av^9+o(av^9).
 \tag{43}
\]
新增保护 (32) 在 halo 中单独付 −40av⁹/243；其余九阶项包括 β 的 −4u、原 host 已有反馈、移动 κ 和真实质量/分母，均已在 (41)–(42) 中合并。

这些是选定实际族的绝对展开。没有用未知 O(av⁹) 拟合系数，也没有把二维纯宏第三反馈直接加到旧微孔族而漏掉同阶 host 保护费。本文未给全输入九阶上界；微孔更高密度形状、不同 host 或更一般宏观孔可能改善 (2)。

## 8. 实际计算、依赖与可复核材料

所有新算核在本目录 computation/。本机实际运行 WolframScript 15.0.0 ARM 的 ninth_budget.wl；最终原始日志 ninth_budget.log 自带实际 PID 81523。独立 Python 3.15 标准库 Fraction 的 ninth_fraction.py 实际运行并通过九条精确断言，原始 stdout 和 JSON 均保留，PID 81866。

两套运算都从 (3)–(5) 的有限有理式独立展开到 u¹¹；结果给 (41)、(42)、(43)，低于 u¹⁰ 的残差为零，八阶 \(\widehat J_2\) 回归残差为零。Wolfram 另核两种 higher-shear 参数任意 u²/u³ 系数在本九阶结果中消失，及 cap 矩 ∫(1+r)r=5/6。没有用计算替代邻孔清单、(26) 的真实矩形误差、正部 (31)、严格阈值或有限化。

首次运行的最后两项回归有输入错误：Coefficient 未先展开 P₇ 的有理部分，及 c2u/c3u 被读成单一标识符。原代码与 stdout 保存在 ninth_budget.first-partial.wl/log，只作失败历史，不作证据；修正后实际重跑，最终两项分别为 0 和 −23295344/2657205。没有补造版本、时间或退出信息。

已全文读取的冻结基础：

- R054 原微孔几何 small_source_eighth.md，SHA bc74ac7f20b63be1387bf1f465316fb61e00513b6567002246ac32095b4f00d1；
- R060 第三反馈及宏微未合并边界 third_feedback_and_micro.md，SHA c7bcd28338661fa782029c87c883467d7ea6dd812d3b16301a02c74598cc9ccc；
- R067 v2 真实收缩与 host 算子 iterated_protection.v2.md，SHA 7e148e1fb837e3c09f595bfcbc5d4d5ec2c27f4d8f8d5486dff1ae29d91ff429；
- R070 真实二维非线性有限层 physical_feedback.md，SHA a83818e4e092977566fc516a2edab5f20639183fa22e908732dbdab3ff499014。

新证明的关键增量为：所有 cap 的显式细格 discrepancy；包含整个短见证 host 的新增保护 leading；有限偶层的正部误差；同一新参数下的预算抵消。既有纯宏有限积分只用到 (18)–(20)，非线性容量曲率相对线性首影响实际比值 av¹⁰，晚于本轮。原完整八阶上界不参与本次单侧构造证明。
