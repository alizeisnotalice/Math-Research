# 尖峰加完整背景：点态高来源占用与来源加权二阶压力

2026-10-07。只新增本稿及同前缀证书，不改主账。使用已读L01/L02的真实输入、总质量、strict阈值及反例范围合同；本稿不将容量LP的放松参数当实际中心cube来源。读取现有source-tail registration与cross_shell_pair_geometry；source-tail正文尚未落盘时按其注册定义固定
\[
 S(y)=\int_E\frac{\tau}{M(x)}h_{R(x)}^{(n)}(x-y)dx,\quad
 I_1=\int S\,d\mu=\tau|E|,\quad I_2=\int S^2d\mu.
\tag{1}
\]
R为完整同一输入的真实finite最大赢家，E={M>τ}。这是soft-free hard来源压力，未认证实际FIRST/fullfuture/CP-GP/history。

待审强合同为 \(I_2\le\sqrt n\,\mathrm{polylog}(n)I_1\)。单峰可以令pointwise S≈n，但下面完整计费与角向二阶几何表明，该候选不否定强合同；本稿亦不证明强合同一般成立。单峰上界只作反例资格审计，不作一般新paid成果。

## 1. 同一输入下的唯一winner，严格完整E

取a=1、J≥n、\(q=1+1/(2J)\)、\(L_j=q^j\)、b=L_J<2，
\[
 \mu=\alpha\delta_0+\beta\mathbf1_{Q(0,2b)}dx,\quad
 \tau=\beta>0,\quad \alpha=\eta\beta,\quad0<\eta\le1.
\tag{2}
\]
背景支撑完整side2b，质量 \(\beta(2b)^n\)，总W=\(\alpha+\beta(2b)^n\)。对于x∈Q(0,b)，所有candidate cube完整处于背景，响应
\[
 h_{L_j}*\mu(x)=\beta+\alpha L_j^{-n}\mathbf1_{Q(0,L_j)}(x).
\]
因此winner唯一为**最小捕获尖峰尺度**：A0=Q1，Aj=QLj\QLj−1，j≥1；边界a.e.约定不影响积分。在Qb外没有candidate捕获尖峰，背景响应≤β，严格阈值不接受。故完整E恰Qb，且
\[
 I_1=\beta b^n,\qquad
 g(x)=\frac{\beta}{\beta+\alpha R(x)^{-n}}.
\tag{3}
\]
背景响应的ties在E外，不用人为tie规则把它们搬进E；E内由α>0严格分出winner。

尖峰的来源占用准确为
\[
 S(0)=\frac1{1+\eta}
 +\sum_{j=1}^J\frac1{1+\eta L_j^{-n}}(1-q^{-n}).
\tag{4}
\]
于是
\[
 \frac{1+J(1-q^{-n})}{1+\eta}\le S(0)
 \le1+J(1-q^{-n})\le1+n\log b\le1+n/2.
\tag{5}
\]
由 \(q^n\ge1+n/(2J)\) 与J≥n，下界≥\((1+n/3)/(1+\eta)\)。确有真实pointwise Θ(n)，没有删门的自由selector。

但尖峰的first占比至多 \(\eta b^{-n}(1+n/2)\)，second占比至多 \(\eta b^{-n}(1+n/2)^2\)。它们已显示不能只取S(0)而忘记背景来源；下一节也控制背景的完整second项。

## 2. 自含cube-cone二阶体积，背景square没有隐藏密度常数

定义cone概率ω（均匀面/符号，其余坐标均匀[-1,1]），对Q≥1
\[
 P_Q(s)=\Pr(|s-U|\le Q),\quad U\sim\mathrm{Unif}[-1,1],\qquad
 B_Q(s)=\tfrac12[\mathbf1_{|s-1|\le Q}+\mathbf1_{|s+1|\le Q}],
\]
\[
 \mathcal L_Q(z)=\frac1n\sum_i B_Q(z_i)\prod_{k\ne i}P_Q(z_k).
\]
这是 \(\Pr_\omega(\|z-\omega\|_\infty\le Q)\)。P_Q在|s|≤Q−1为1，在[Q−1,Q+1]线性下降为0。直接一维积分给
\[
 A_Q:=\int P_Q^2=2Q-\tfrac23,\quad
 C_Q:=\int B_Q^2=\int B_QP_Q=2Q-1.
\]
展开平方，独立坐标积分逐项相乘，n≥2时
\[
 \int\mathcal L_Q^2
 =\frac{C_QA_Q^{n-1}}n+
    (1-\frac1n)C_Q^2A_Q^{n-2}.
\tag{6}
\]
在Q=1，变元z=2y/r得
\[
 \int\mathcal L_1(2y/r)^2dy
 =c_n(2r/3)^n,\qquad c_n=\tfrac9{16}(1+\tfrac1{3n})\le\tfrac34.
\tag{7}
\]
这里是对Lebesgue dy的exact geometry；将它用于背景是因为背景密度**实际等于β**，不替任意μ假设L∞上界。

本finite grid有q−1≤1/(2n)。用 \(C_q/C_1=1+2(q-1)\)、\(A_q/A_1=1+\tfrac32(q-1)\)，(6)的两项分别控制，得
\[
 \int\mathcal L_q^2\le16\int\mathcal L_1^2.
\tag{8}
\]
自含粗常数：比值≤\((1+1/n)^2(1+3/(4n))^n\le4e^{3/4}<16\)。只需e<3，未引用未读的高维cone定理。

将g≤1先用于资格审计。core term为
\(\prod_i(1-|y_i|/a)_+\)，其L2 norm²=(2a/3)^n。对shell \(r\in(L_{j-1},L_j)\)，x=rω/2，Jacobian为 \(nr^{n-1}dr\,d\varsigma\)，而R=Lj，\(L_j/r\le q\)。因此
\[
 S(y)\le\prod_i(1-|y_i|/a)_+
      +n\int_a^b\mathcal L_q(2y/r)\frac{dr}r.
\tag{9}
\]
Minkowski、(7)/(8)及\(\sqrt{16c_n}<4\)给
\[
 \|S\|_{L^2(dy)}
 \le(2/3)^{n/2}\{a^{n/2}+8(b^{n/2}-a^{n/2})\}
 \le8(2b/3)^{n/2}.
\]
所以完整背景项（先扩到全dy，正性合法）
\[
 \beta\int_{Q(0,2b)}S(y)^2dy\le64\beta(2b/3)^n.
\tag{10}
\]
没有对sources逐坐标uniform后验、没有把薄Lebesgue壳当μ概率，也没有暗收J。

## 3. 这族并未否定来源加权强合同

由(3)/(5)/(10)，
\[
 \boxed{\frac{I_2}{I_1}
 \le64(2/3)^n+\eta b^{-n}(1+n/2)^2.}
\tag{11}
\]
J≥n而b=q^J随J趋e^1/2；此族的点态尖峰S~n与来源平均second之间存在指数体积补偿。无需将背景当免费输入。甚至 \(I_1/W=\beta b^n/[\eta\beta+\beta(2b)^n]\) 本身约2^-n。它不能作为一般强合同反例；单峰失败也不能反向证明强合同。

增大单峰质量也不能免费保留(5)的Θ(n)profile。对任意η>0，(4)的shell和是对递减函数1/(T+η)的right Riemann sum，故
\[
 S(0)\le\frac1{1+\eta}+\log\frac{b^n+\eta}{1+\eta}.
\]
记t=η/b^n，则右边≤\(1/(1+t)+\log(1+1/t)\le3/(2\sqrt t)\)。最后使用 \(1+t\ge2\sqrt t\)，以及由 \(1+s\ge2\sqrt s\)积分得到的 \(\log(1+x)\le\sqrt x\)。因此单峰贡献 \(\alpha S(0)^2/(\beta b^n)=tS(0)^2\le9/4\)。这只是参数质量tradeoff的解析审计，未包含在82项有限守卫中，也不推广为多峰或一般来源界。

可转成真L1窄峰：把αδ0改为均匀α质量于Q(0,2ε)，背景不变，按全部新响应重判winner。固定n,J后ε→0，除finite shell面零测外所有响应、strict E与uniquewinner稳定；来源y=εV的Sε(y)→S(0)，背景source的Sε(y)→S(y) a.e.。Eε在Q(0,b+2ε)，column有原sup-hard包络1+nlogb，故正dominated convergence证明I1/I2收敛。该定维转移不声称ε固定对所有n一致，三轮guard也不冒称执行连续L1积分。

## 4. 多尖峰共享背景：哪些容量不能免费加

令 \(\rho=\sum_k\alpha_k\delta_{z_k}\)、背景β1_D。完整winner必须由真实background捕获分数一同重算：
\[
 M(x)=\max_j L_j^{-n}
       [\beta|Q(x,L_j)\cap D|+\rho(Q(x,L_j))],\quad
 g(x)=\frac{\beta}{M(x)}
\tag{12}
\]
不能把(4)孤立源profile逐峰相加：新增点会改变ρ捕获量、真实winner和g。若要求每个selected query的background平均恒β，必要support为 \(D\supset U:=\bigcup_{x\in E}Q(x,R(x))\)，故背景费用至少β|U|。有 \(E+Q_a\subset U\subset E+Q_b\)，所以至少β|E+Q_a|；**不**能据此强加β|E+Q_b|。D⊃E+Q_b是方便但更强的充分选择，它保证所有candidate背景平均β，(12)此时才可简化为β+max spike平均。使用规范有限网/可测selector时U的体积可按可测外包理解；本文不依赖未证明的一般selector并集Borel性。

这一必要费用与充分选择的区分是多峰审计的关键；若原LP没有β|U|及候选query的真实截断背景，就是外放松，但不能因它没有更强的b膨胀费而判其不可行。

例如想维持g≥1/(1+η)，需在所有selected query上有真正共同capacity
\[
 \rho(Q(x,R(x)))\le\eta\beta R(x)^n.
\tag{13}
\]
它仍只给spike的first贡献≤ηI1、second粗界≤η(1+nlogb)I1，尚不够sqrt n。故不把(13)当新行条件进展。

更具体的几何检查是每个源的“孤立首次捕获”方向：若希望源zk在r-shell上复用单峰winner机制，必须在相应cone方向中排除其它源：
\[
 \|2(z_\ell-z_k)/r-\omega\|_\infty>1
 \quad(\ell\ne k),
\tag{14}
\]
或完整捕获质量增长仍不足以使更大cube胜出。(14)的**所有**方向条件是交集，不能把各点的单独角概率相加为来源付款。若删除它，却继续用单峰(4)，构造已不是同一真实winner输入。多峰若能共享背景且保留足够多(14)方向，才可能摆脱(11)的指数补偿；本文尚未证明该packing不可能，也没有确认一般反例。B38容量搜索下一步必须保留(12)、D的费用及源方向交集，而不只优化非负质量表。

## 5. 三轮新压缩资格守卫

先登记n16/64/256、J=n²、η=1/100，仅Fraction检查source/winner参数、pointwise Θ(n)上下常数、cone exact体积比、背景与spike完整质量，以及(11)压缩比值。避免枚举J个巨大分母profile；使用(5)解析预算与rational powers，b=q^J保持symbol，只用Bernoulli的b≥3/2与b<2，不计算q^(Jn)大整数。无MC，无旧数据重跑。

终态为 **82/82 Fraction PASS**，每轮27项加注册1项，约0.15秒。压缩(11)上界 \([64+\eta(1+n/2)^2](2/3)^n\) 三轮显示值分别0.0986693、4.02331×10^-10、1.91928×10^-43；近似不参与任何断言。它们只核原中心几何公式、完整背景质量与单峰反例资格，不证明一般source-tail强合同或actual门。没有连续L1积分或多峰空间样本。

专属 [registration](source_spike_background_tail_audit_20261007_registration.json)、[script](source_spike_background_tail_audit_20261007_exact_guard.py)、[results](source_spike_background_tail_audit_20261007_results.json)已保存。执行script SHA256为dc8a75a3cb195d34f74efe2e5d18ec562f25875f3528c455bfaabea833a2a89a。无live handle，未改已执行script。§4随后纠正selected-query并集U与更强b膨胀条件的区别，不影响单峰解析或守卫，不为该文案修订重跑。
