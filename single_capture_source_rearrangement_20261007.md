# 真实单捕获输出的来源二阶重排：多峰常数账与合作交叉

2026-10-07。独立解析审计；只新增本稿，不改主账，不运行toy或旧守卫。沿用已读L01/L02的完整来源/真实winner/资格审计合同。查重source_partition_dominance、greedy witness、gated_occupancy_bridge与deep_gram_bridge：旧dominance付交通与Γ能量，单原子hardband付短壳；没有据此自动得到本稿的来源column二阶量。本稿不对同一actual事件再领旧receipt。

主结论是root所提重排完全正确，且可以对任意多峰的**真实only-one-captured输出**求和，常数不依n或峰数。它只支付spike来源的isolated profile；共同背景来源、cooperative profile及iso/coop交叉均须单列，不宣称全geom已付。

## 1. 完整输入与精确单捕获门

设 \(\rho=\sum_i\alpha_i\delta_{z_i}\)、\(\alpha_i>0\)、\(\sum_i\alpha_i<\infty\)，
\[
 \mu=\rho+\beta\mathbf1_Ddx,\quad\beta>0,\quad
 M(x)=h_{R(x)}*\mu(x),\quad E=\{M>\beta\}.
\tag{1}
\]
R为同一输入的真实finite winner（或已经给定的合法最大selector），不重选原门。假定所研究的输出上selected \(Q(x,R(x))\subset D\)，使背景平均**实际恒等β**。这个条件只要求D含selected query并集U，不强加E+Q_b。用原closedcube/可测约定，尺度有限时所有捕获集是Borel；边界Lebesgue零测可沿原tie规范保留。

令
\[
 E_i=\{x\in E:\rho(Q(x,R(x)))=\alpha_i,\
              Q(x,R(x))\text{捕获的正质量spike仅为 }z_i\}.
\tag{2}
\]
重复位置先合并成一个正质量原子；否则“仅一spike”标签不唯一。E_i不交，\(E_{\rm iso}=\bigcup_iE_i\)，\(v_i=|E_i|\)。来源支撑有限与R≤b使E有限；更一般直接假定当前E有限，符合原审计窗口。

在E_i上，真实M=\(\beta+\alpha_i/R^n\)，故原threshold normalization下来源profile为
\[
 S_i^{\rm iso}=\int_{E_i}
  \frac{\beta}{M(x)}h_{R(x)}(x-z_i)dx
 =\int_{E_i}\frac{\beta}{\beta R(x)^n+\alpha_i}dx.
\tag{3}
\]
其它spike在E_i的source kernel为0。因此**整个isolated receiver部分**的spike-source square，恰为
\(\int(S^{\rm iso})^2d\rho=\sum_i\alpha_i(S_i^{\rm iso})^2\)；不是把每个孤立单峰输入的profile免费叠加。

## 2. 径向重排的自含证明

捕获门保证 \(R(x)\ge2\|x-z_i\|_\infty\)。定义
\[
 V_i(x)=(2\|x-z_i\|_\infty)^n,\qquad
 \phi_i(u)=\frac{\beta}{\beta u+\alpha_i}.
\]
精确 \(|\{V_i\le u\}|=u\)，φ_i递减。对任意Lebesgue集合A、|A|=v，层饼给
\[
 \int_A\phi_i(V_i(x))dx
 =\int_0^{\beta/\alpha_i}|A\cap\{\phi_i(V_i)>t\}|dt
 \le\int_0^{\beta/\alpha_i}
       \min(v,|\{\phi_i(V_i)>t\}|)dt
 =\int_0^v\phi_i(u)du.
\tag{4}
\]
最后等号由相同层饼换序，或直接取中心cube \(\{V_i\le v\}\)。不需要球对称定理、来源uniform、a→0或假设E_i本身是cube。由(3)/(4)，
\[
 \boxed{S_i^{\rm iso}\le
      \log(1+\beta v_i/\alpha_i).}
\tag{5}
\]
v_i=0给0；α_i>0避免除零，零质量项忽略。有限v无额外核可积问题。\(1+s\ge2\sqrt s\)在s>0积分给
\[
 \log(1+t)=\int_0^t\frac{ds}{1+s}\le\sqrt t.
\]
因此
\[
 \boxed{\sum_i\alpha_i(S_i^{\rm iso})^2
 \le\beta\sum_i v_i=\beta|E_{\rm iso}|.}
\tag{6}
\]
非负Tonelli使(6)适用于countable原子集合，无n、atom-count或来源深度费用。严格阈M>β保证(3)的真实归一化；常数证明本身不依winner极大性，只依捕获、完整背景和原M。

## 3. L1 full-packet扩展与准确支付范围

原子不是必要条件。预先固定来源分解 \(\rho=\sum_i\rho_i\)、\(M_i=\|\rho_i\|>0\)，无输出依赖的新归一化。若E_i的selected query **完整捕获ρ_i，且不捕获其它正packet**，则真实M=\(\beta+M_i/R^n\)。对ρ_i-a.e.原来源y，捕获给R≥2||x−y||∞，逐y重排得
\[
 S^{\rm iso}(y)\le\log(1+\beta v_i/M_i),\quad
 \int(S^{\rm iso})^2d\rho_i\le\beta v_i.
\tag{7}
\]
可用真正L1窄packet，无需先令宽度为0。全捕获是实质门，不能把partial capture自动替成完整M_i。若只知道partial捕获≥ηM_i，可对**预先不交receiver分配及相应来源子交通**得βv_i/η，但未分配的交叉packet source交通仍保留；本稿不借此再报一个旧dominance的主账费用。

全μ的square还含
\[
 \beta\int_D(S^{\rm iso}(y))^2dy,
\]
这个共同background-source项没有由(6)/(7)支付。每个spike在iso与coop输出都可贡献，故整体只能用
\[
 \int(S^{\rm iso}+S^{\rm coop})^2d\rho
 \le2\int(S^{\rm iso})^2d\rho+
       2\int(S^{\rm coop})^2d\rho.
\tag{8}
\]
不允许把来源square按receiver不交性直接相加。若后续能独立支付coop，则(8)只引入常数2；当前coop仍未知。本稿亦没有将β|Eiso|写成来源W的一次traffic费用：它是针对source-tail合同I2/I1的有偿receiver项，I1=\(\beta|E|\)，原weak体积尚须其它论证。

## 4. 按真实捕获集合/总质量的合作block：diagonal常数，交叉未付

有限原子输入时，按真实selected capture集合 \(\mathcal A(x)=\{i:z_i\in Q(x,R(x))\}\)定义不交E_A，非空A；设
\[
 m_A=\sum_{i\in A}\alpha_i,\quad v_A=|E_A|,\qquad
 s_{i,A}=\int_{E_A}\frac{\beta}
                  {\beta R(x)^n+m_A}dx\quad(i\in A).
\tag{9}
\]
每个s_{i,A}保留原R与所有实际捕获，且核支持捕获i。逐i的(4)给
\[
 s_{i,A}\le\ell_A:=\log(1+\beta v_A/m_A),\qquad
 \sum_{i\in A}\alpha_i s_{i,A}^2
 \le m_A\ell_A^2\le\beta v_A.
\tag{10}
\]
于是所有真实capture-block **diagonal** 之和≤β|E|。它按总捕获质量m_A归一化，不按|A|分层，也不要求每source质量相等。

但是实际 \(S(z_i)=\sum_{A\ni i}s_{i,A}\)，所以精确square仍为
\[
 \int S^2d\rho
 =\sum_{A,B}\sum_{i\in A\cap B}
                  \alpha_i s_{i,A}s_{i,B}.
\tag{11}
\]
为明确缺口而不删交叉，记共同来源单位向量
\[
 e_A(i)=\sqrt{\alpha_i/m_A}\mathbf1_{i\in A},\quad
 \Gamma_{AB}=\langle e_A,e_B\rangle
       =\frac{m_{A\cap B}}{\sqrt{m_Am_B}},\quad
 b_A=\sqrt{m_A}\ell_A\le\sqrt{\beta v_A}.
\]
正性仅给
\[
 \int S^2d\rho
 \le\sum_{A,B}b_Ab_B\Gamma_{AB}
 =\left\|\sum_A b_Ae_A\right\|_{\ell^2}^2.
\tag{12}
\]
Γ是真正同源capture-overlap Gram；它不同于deep_gram_bridge中带LCA/远对门的有符号核，本文没有把后者PSD化。receiver的E_A不交不保证e_A正交，同一原source可被许多真实输出集合捕获。只逐block使用(10)后套Cauchy会引入capture-set数量，违反无atom-count要求；只按m_A的dyadic档合并也不删除同source跨档交叉。

因此本次合作推导给出真实有偿diagonal及精确源重叠缺陷；未得到Γ作用于实际b向量的dimension预算。需要利用中心cube实际可实现的capture集合及winner，证明weighted overlap费用，而不能从任意非负incidence matrix的谱性免费推出。没有将“无法证明”当反例，也不制造未认证actual数据。

### 4.1 保留log非线性，近重复block不能假收个数费

(12)保留精确 \(b_A=\sqrt{m_A}\log(1+\beta v_A/m_A)\)，不是将全部b_A替成√βv_A后声称没有细分代价。当v_A≪m_A/β时，b_A≤βv_A/√m_A，receiver细分近重复block的贡献是volume线性。比如一个真实receiver子集B的所有捕获质量≥m0，Minkowski直接给
\[
 \int(S_B)^2d\rho\le\beta^2|B|^2/m_0.
\tag{13}
\]
可在源空间直接证明：每行的L2(ρ)norm为
\(\beta\sqrt{m(x)}/(\beta R^n+m(x))\le\beta/\sqrt{m_0}\)，再积分dx。它不依block数、atom数，说明任意复制near-identical capture标签的fakeJ收费是错误的。但当m0≪β|B|时该式仍大，不能据此支付所有合作交通；也不删除同source跨真实质量层的cross。

### 4.2 真实大份额source incidence，square费用≤β|E|/θ

固定θ∈(0,1]，只保留source-receiver incidences
\[
 z_i\in Q(x,R(x)),\qquad\alpha_i\ge\theta m(x),
 \quad m(x)=\rho(Q(x,R(x))).
\]
记相应receiver集合 \(E_i^\theta\)及原restricted profile
\[
 S_i^\theta=\int_{E_i^\theta}\frac{\beta}
                  {\beta R(x)^n+m(x)}dx.
\]
因m(x)≥α_i，逐i重排仍给
\[
 \alpha_i(S_i^\theta)^2\le\beta|E_i^\theta|.
\]
每个receiver的这些正spike满足 \(\sum_{\mathrm{dom}\ i}\alpha_i\le m(x)\)，所以个数≤1/θ。准确地
\[
 \boxed{\sum_i\alpha_i(S_i^\theta)^2
 \le(\beta/\theta)|E|.}
\tag{14}
\]
θ=1就是single-capture；θ=n^-1/2给spike-source square的sqrt n receiver预算。费用无需先按capture set再付count。剩余profile仅保留每次捕获中α_i<θm(x)的source incidence；同source大/小份额profile之和的square要用(8)类常数2，而非无交叉相加。

这是一份source-column二阶合同，原source-dominance/greedy主要是完整μ hardband的输出traffic和Γ能量合同。尤其此θ分母是spike捕获ρ，不是完整μ捕获；背景可占大多数，所以不能把它换成旧完整μ dominance子事件并额外领actual费用。两者思路相关，但量与资格不同。

原子α_i阈值会随人为细分改变；把一个原子拆为多个同位置标签后，各α_i小，并不说明输入几何分散。因此(14)的small incidence补集不可冒称几何nonconcentration。若推广，必须预先固定真正来源packet并保留完整捕获，或同时要求捕获量≥ηM_i。后者对不重复来源分解 \(\rho=\sum\rho_i\)及其真实捕获量 \(m_i(x)\ge\theta m(x), m_i(x)\ge\eta M_i\)，以同一restricted source kernel可得β|E|/(ηθ)；理由是m(x)≥ηM_i，逐来源y重排后每packet square≤β|E_i|/η，每row packet数≤1/θ。未满足packet内部捕获比例时不能把partial m_i升级为总M_i。任意输出重组packet或重领来源均不合法。

若来源分解是fractional而非支撑不交，写 \(\rho_i=a_i\rho,\sum_i a_i\le1\)，把未分配份额作为零接受标签。原unlabelled restricted profile为 \(\sum_i a_i(y)S_i(y)\)，其square由Jensen≤\(\sum_i a_i(y)S_i(y)^2\)；积分才可用各ρ_i预算。不能把fractional labelled profiles直接相加后平方却仍只付sum labelled square。

## 5. 结论与查重范围

(6)/(7)是多峰/packet真实single-capture来源二阶的常数预算；与旧source-dominance交通receipt是不同量，不能各自给同一actual支加一次费用。它可作为source-tail证明中的一项互补控制，尚不支付共同背景或完整actual R_angle。

合作部分保留(11)/(12)的真实cross项及log非线性；(14)能source-once控制真实大质量份额incidence的square，却不把剩余任意fine标签解释成几何分散。小份额source重叠仍须中心cube约束，不是每个capture-block再加一个行条件。本文没有运行toy或重复旧数值；(14)只是直接解析一般质量预算，不声称已测actual资格。全FIRST、fullfuture、LCA/history若要迁移仍须保留原来源端点和原kernel；仅名字相同的“isolated”不够。
