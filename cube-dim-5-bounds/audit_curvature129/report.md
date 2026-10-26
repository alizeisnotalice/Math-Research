# audit_curvature129：锚点时间后验曲率的独立理论审计

**结论：所提新公式成立。** 对本报告固定模型的每个正截时、每个实锚点和任意实位移，有
\[
 \boxed{\log H_\epsilon(c+u)\ge \log H_\epsilon(c)
 +\nabla\log H_\epsilon(c)\cdot u-q_\epsilon(c)|u|^2,\qquad
 q_\epsilon(c)=\frac{\int_\epsilon^\infty p_t(c)t^{-2}dt}{16H_\epsilon(c)}.}
 \tag{A}
\]
这包括完整 Green 的全部奇异子空间。证明使用固定锚点的概率权重和指数 Jensen；没有沿空间路径冻结后验期望，也没有把随时间变化的指数因子从积分中取出。075 §5 的警告仍然正确，与本结论不冲突。

对 \(\epsilon_0=1/128,\epsilon_1=1/64\)，若持有同一实锚点的有效端点
\[
 0<L_{64}\le H_{1/64}(c),\qquad H_{1/128}(c)\le U_{128}<\infty,
\]
则合法的单锚点二次罚项是
\[
 \boxed{q_\epsilon(c)\le q_U(c):=
 1024\left(1-\frac34\frac{L_{64}}{U_{128}}\right),\qquad
 256\le q_U(c)<1024.}\tag{B}
\]
原 H128 的 Delta 和 log-gradient 区间可以原样使用，16 个实际闭角点仍足够。该理论结论不证明任一实际新盒通过，不是完整 g110 认证，也不改变正式 \(C_{\rm cert}\le C_5\le32\)。主代理保留最终整合和证明责任。

## 1. 审计范围和实际来源

本审计从实际来源的正产品 Cauchy 积分重新推导，不以作者的 passed 标签作依据。读取了根 AGENTS、clone 的 STATE 和主目录最新 STATE；后者含尚未提交的014状态，与独立 clone 的已提交快照不同，二者都已分别登记。独立 clone 基于 `dd413e1cf61a0e376d27ba45fdbaf7ac27fc4a0c`，`--no-hardlinks`，仅 checkout 本项目及根文件；没有提交、推送或修改其他项目目录。

| 实际来源 | SHA256 |
|---|---|
| `audit_frame075/report.md`，完整纤维、微分及警告 | `3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225` |
| `audit_mixed099/report.md`，异质截时与闭角点接口 | `96e874223e5fbc78f9688f7eec574cdcdde184344b1fa809114bc87553d48e85` |
| `upper_alpha1_029/report.md`，完整一般 L1 传递 | `298f0de7277ea6a7a516b85344c5189562bf7026f3f8de43d4d3f304ff6fed81` |
| `upper_shorttime076/report.md`，完整短时与奇异集接口 | `9468f82c4870472a573f608e9d121316c30f0f7ced857e35dbe62ed8dd0f75f8` |
| `primary/finite_source_transfer.md`，原始有限源证明 | `37bf79a4d260c853781881fd371a2ef366fd97a8dae952b0030066b73d8d846f` |

以上五个源在 clone 与主目录读取时字节一致。登记含全部绝对路径及 SHA，不复制旧核、树或大状态。本报告没有重新批量读取098/105证书，没有独立重算全量 ancestry；这些实际消费门须由对应作者/数值独审另闭合。未读取或使用另一位128作者的审计结论。

## 2. 固定模型及每条纤维的规范密度

固定五维 Lebesgue 测度、中心轴平行 cubes、所有坐标同一个正半径、一般 L1 输入。令
\[
 V=\{(1,\pm1,\pm1,\pm1,\pm1)\},\quad B=[v_1\cdots v_{16}],\quad
 A=B^T/16,\quad K=\ker B.
\]
方向不单位化，16 个方向时间均为 \(t\)、权均为1，Fourier 逆因子为 \((2\pi)^{-5}\)。符号交叉和给
\[
 BB^T=16I_5,\quad BA=I_5,\quad |Au|^2=|u|^2/16.
\]
定义
\[
 s_t(z)=\frac{t}{\pi(t^2+z^2)},\quad f_t(z)=\prod_{j=1}^{16}s_t(z_j).
\]
16个独立 Cauchy 的特征函数推前为
\(\widehat p_t(\xi)=\exp[-t\sum_v|v\cdot\xi|]\)，故与上述 Fourier 模型一致。

K 上使用11维 Euclidean Hausdorff/Lebesgue 测度。\(\sqrt{\det(BB^T)}=1024\)，所以 coarea 公式为
\[
 p_t(x)=\frac1{1024}\int_K f_t(Ax+w)dw.\tag{1}
\]
为排除仅 a.e. 的异常纤维版本，固定取全正列和四个各翻转一个尾坐标的列组成 C，其 \(|\det C|=16\)，其余列组成 D。用剩余11坐标 y 参数化纤维，其 Euclidean 体积因子为
\(\sqrt{\det(I+D^TC^{-T}C^{-1}D)}=64\)，于是处处的正积分为
\[
 p_t(x)=\frac1{16}\int_{\mathbb R^{11}}
 \prod_{j=1}^{5}s_t((C^{-1}(x-Dy))_j)\prod_{j=1}^{11}s_t(y_j)dy.
\]
后11个因子形成概率测度，前5个各不超过 \((\pi t)^{-1}\)，因此
\[
 0<p_t(x)\le\frac1{16\pi^5t^5},\qquad
 0<H_\epsilon(x):=\int_\epsilon^\infty p_t(x)dt
 \le\frac1{64\pi^5\epsilon^4}.\tag{2}
\]
相同支配给每个正时间的连续性，故(1)确实表示处处规范的有限连续密度。029 另用 det48 子基给更小的粗峰值常数，两者都是有效上界，不存在归一化矛盾。

## 3. 可微性、固定时间 Hessian 和全空间二次式

写 \(\ell_t=\log s_t\)，直接微分得
\[
 \ell'_t(z)=-\frac{2z}{t^2+z^2},\qquad
 \ell''_t(z)=\frac{2(z^2-t^2)}{(t^2+z^2)^2},\qquad
 |\ell'_t|\le t^{-1},\quad |\ell''_t|\le2t^{-2},\quad
 \ell''_t\ge-2t^{-2}.\tag{3}
\]
这里不要求任何 Cauchy 绝对矩。对方向 u 固定取 \(a=Au\)，在同一个 K 上平移，有
\[
 p_t(x+su)=1024^{-1}\int_K f_t(Ax+w+sa)dw.
\]
令
\(\sigma_a=\sum a_j\ell'_t(z_j)\)、\(\tau_a=\sum a_j^2\ell''_t(z_j)\)。它们满足
\[
 |\sigma_a|\le\|a\|_1/t,\quad |\tau_a|\le2|a|^2/t^2,\quad
 \tau_a\ge-2|a|^2/t^2.
\]
对 \(|s|\le R,t\ge\epsilon>0\)，有
\[
 f_t(z+sa)\le e^{R\|a\|_1/\epsilon}f_t(z).
\]
零、一、二阶空间导数由该右侧再乘仅依赖 \(a,R,\epsilon\) 的有限常数控制，因为二阶导数为 \(f_t(\sigma_a^2+\tau_a)\)。由(2)，这个原被积函数在 \([\epsilon,\infty)\times K\) 上总质量有限；固定 t 时的纤维积分同样有限。因此 DCT 合法交换一次、二次空间导数与纤维/完整时间积分。局部混合导数亦有同样支配，故 \(p_t,H_\epsilon\) 为 C2，导数连续。端点不随 x 移动，没有空间微分的时间端点项。

固定 \(t>0,x\)，以 \(f_t(Ax+w)/(1024p_t(x))\) 归一化 K 上概率，得到
\[
 D_u^2\log p_t(x)=\mathbb E\tau_a+\operatorname{Var}(\sigma_a)
 \ge-2|a|^2/t^2=-|u|^2/(8t^2).\tag{4}
\]
方差有限且非负；纤维条件后的坐标无需独立。由固定 t 的全空间界(4)，沿 \(c+su\) 积分二阶余项，得到
\[
 \log p_t(c+u)\ge\log p_t(c)+\nabla\log p_t(c)\cdot u
 -|u|^2/(16t^2).\tag{5}
\]
此处 \(1/8\) 是 Hessian 常数，\(1/16\) 是二次罚项，因子2无误。

## 4. 固定锚点概率与 Jensen：新公式的完整证明

对任意实 c 定义一次性固定的时间概率
\[
 d\nu_c(t)=\frac{p_t(c)}{H_\epsilon(c)}dt,\quad t\in[\epsilon,\infty).
\]
其归一化、正性和有限性已由(2)给出。令
\[
 X_t=\nabla\log p_t(c)\cdot u-|u|^2/(16t^2).
\]
由上节的 score 公式，
\[
 |\nabla\log p_t(c)\cdot u|\le\|Au\|_1/t\le |u|/t,
\]
因 \(\|Au\|_1\le4|Au|=|u|\)。所以 X 在整个时间区间上有界，\(X,e^X\) 都对 \(\nu_c\) 可积；不需要假定 \(\log p_t(c)\) 自身可积，也不用概率重尾矩。固定锚点的一阶导数交换给
\[
 \mathbb E_{\nu_c}\nabla\log p_t(c)
 =\frac{\int_\epsilon^\infty\nabla p_t(c)dt}{H_\epsilon(c)}
 =\nabla\log H_\epsilon(c).\tag{6}
\]
把(5)指数化、按 \(\nu_c\) 积分，再对指数函数使用 Jensen，有
\[
 \frac{H_\epsilon(c+u)}{H_\epsilon(c)}
 =\mathbb E_{\nu_c}\frac{p_t(c+u)}{p_t(c)}
 \ge\mathbb E_{\nu_c}e^{X_t}
 \ge e^{\mathbb E_{\nu_c}X_t}.
\]
两侧都是正有限实数，取 log 并代入(6)即得(A)。其中
\[
 0<\int_\epsilon^\infty p_t(c)t^{-2}dt
 \le\epsilon^{-2}H_\epsilon(c)<\infty,
\]
故 q 正且有限。实际上后一个不等号严格，因为 \(p_t(c)>0\) 且 \(t^{-2}<\epsilon^{-2}\) 于 \(t>\epsilon\)。

也可完全不先求 \(D^2\log p_t\)：由(3)直接得到
\[
 f_t(Ac+w+Au)\ge f_t(Ac+w)
 \exp\{\sigma_{Au}(t,Ac+w)-|u|^2/(16t^2)\}.
\]
对固定 anchor 的联合概率
\(f_t(Ac+w)dt\,dw/(1024H_\epsilon(c))\) 使用同一次 Jensen，
其 \(\sigma\) 期望等于 \(\nabla\log H(c)\cdot u\)，再得到(A)。这是相同结论的直接产品层验证。

075 §4 已有点态
\(D^2\log H(x)[u,u]\ge-|u|^2\mathbb E_x[t^{-2}]/8\)。其中 \(\mathbb E_x\) 随 x 变化，不能直接在整条路径上换成 \(\mathbb E_c\)。本证明没有作这个替换。075 §5 警告不能从 \(\int e^{|x|^2/(16t^2)}p_t(x)dt\) 抽出 t 依赖修正，本证明也没有作这一步；是在 anchor 概率下对已经有界的指数变量使用 Jensen。因此旧警告不是新式的 gap。

## 5. 有限截时序列、端点方向和严格改善

记 \(K_{\epsilon_0}(c)=\int_{\epsilon_0}^\infty p_t(c)t^{-2}dt\)。对有限序列
\(0<\epsilon_0<\cdots<\epsilon_m<\infty\)，令
\(a_j=\epsilon_j^{-2}\)，并令阶梯函数在
\([\epsilon_j,\epsilon_{j+1})\) 取 \(a_j\)，最后尾部取 \(a_m\)。它点态不小于 \(t^{-2}\)。非负积分及有限和给
\[
 K_{\epsilon_0}(c)\le a_0H_{\epsilon_0}(c)
 -\sum_{j=1}^m(a_{j-1}-a_j)H_{\epsilon_j}(c).\tag{7}
\]
每个减项系数正；它是对同一 anchor 的尾质量积分分段重组，没有微分下端，也没有截断无限尾。由末尾 \(t>\epsilon_m\) 的严格差和处处正密度，(7)事实上严格。

若持有 \(0\le L_j\le H_{\epsilon_j}(c)\) 与
\(H_{\epsilon_0}(c)\le U_0<\infty\)，则
\[
 q_{\epsilon_0}(c)\le\frac1{16}\left(a_0-
 \sum_{j=1}^m(a_{j-1}-a_j)\frac{L_j}{U_0}\right).\tag{8}
\]
方向是上界：正的真实比值 \(H_{\epsilon_j}/H_{\epsilon_0}\ge L_j/U_0\)，它在(7)除以 \(16H_{\epsilon_0}\) 后带负号。一般式若所有 L 合法，右侧至少 \(a_m/16>0\)，因为各 \(L_j\le H_{\epsilon_0}\le U_0\)。

特取128/64，\(a_0=16384,a_1=4096,a_0-a_1=12288\)，除以16正好给(B)。更明确地，设 \(H_0=H_{1/128},H_1=H_{1/64}\)，则端点替换残差为
\[
 q_U-\left(1024-768\frac{H_1}{H_0}\right)
 =768\frac{U_{128}(H_1-L_{64})+L_{64}(U_{128}-H_0)}{H_0U_{128}}\ge0.\tag{9}
\]
并且
\[
 q_U-256=768\frac{U_{128}-L_{64}}{U_{128}}\ge0,
 \qquad 1024-q_U=768\frac{L_{64}}{U_{128}}>0.\tag{10}
\]
所以新罚项严格小于旧1024，且保持正性；不能把新式解释为 H128 的真实 Hessian 处处不小于 \(-2q_U(c)I\)。它是(A)在固定 c 的全空间下支撑式。

**L64 的实际资格必须来自 H64。** 取 \(L_{64}=21/200000\) 时，c 落在某个同模型、同坐标的已认证 H64 整个闭叶中即可；在内部是更强条件，闭边界也合法。完整 g105 覆盖为 mixed 截时，故仅凭 \(G\ge21/200000\) 或某个 H128 叶的证书不能推出 H64(c) 的这个下界。必须核对真正的旧 H64 叶、确切几何包含关系、相同 anchor 和对称映射。这里不需要旧叶锚点与 c 相同，因为消费的是旧证书已覆盖整叶的函数下界。

U128 必须是同一 c 的 H128 完整时间尾上端，不可用 H64 上端替代。若 \(L_{64}>U_{128}\)，应判源/绑定不一致，不能把结果静默截为256。准确有理端点可直接用(B)作 Fraction 运算；若用数值区间，要把 \(L/U\) 向下包、最终 q 向上包。必要时与旧1024取最小值仍是合法上界，但不能替代源一致性检查。

## 6. 原 Delta、梯度与16个闭角点

取 H128 同一 anchor 的合法保存数据
\[
 \Delta_L\le\log(H_{128}(c)/g),\qquad
 \partial_j\log H_{128}(c)\in[\gamma_{j,L},\gamma_{j,U}],\quad j=2,\ldots,5.
\]
在正面位移 \(u_1=0\)，由(A)和 \(q_U\ge q\) 有
\[
 \log(H_{128}(c+u)/g)\ge
 \Phi(u):=\Delta_L+\sum_{j=2}^5\min(\gamma_{j,L}u_j,\gamma_{j,U}u_j)
 -q_U(c)\sum_{j=2}^5u_j^2.\tag{11}
\]
qU 在这个证书中只由 anchor 确定，是固定非负常数。每项 min 是凹函数，负二次亦凹，所以 \(\Phi\) 全局凹。闭四维矩形每点是最多16角点的凸组合，因而
\[
 \min_{\text{16 actual closed corners}}\Phi\ge0
 \Longrightarrow H_{128}\ge g\text{ 在整闭盒上};
\]
角点余量严格正则全盒严格大于 g。退化边和重复角点无害。穿越坐标符号面、锥节点重合面或 S2 不影响该判据。

Delta 和梯度始终是 H128 数据；H64 下界只进入 qU 的尾质量比。没有混合两个截时的梯度，没有重新计算 log，没有以 L64 代替 H128 的 log/gradient 分母。一个证书使用一个 anchor，不能对不同角点临时切换 anchor 而仍声称同一个凹式覆盖全盒。

对同一角点，新旧保存下界差准确为
\[
 \Phi_{\rm new}(u)-\Phi_{\rm old}(u)=(1024-q_U)|u|^2.
\]
非零位移严格改善；u=0没有曲率改善。旧待证盒不会因此自动通过，仍须用准确实际闭角点重算最差余量，并保持未通过盒在 pending 集中。理论上非负余量已足够得到 \(\ge g\)；若生产约定使用严格正门，应保持严格门，不以浮点正号代替外向余量。

## 7. 奇异子空间、共同尺度及一般 L1 迁移

所有新步骤只在 \(t\ge\epsilon>0\) 上使用处处正、有限、C2 的 \(p_t,H_\epsilon\)。完整
\(G(x)=\int_0^\infty p_t(x)dt\in(0,\infty]\) 在有限并的 S2 奇异平面上可为无穷；这不在微分或 Jensen 中出现。对每个点都有扩展非负身份
\[
 G(x)=H_\epsilon(x)+\int_0^\epsilon p_t(x)dt\ge H_\epsilon(x).
\]
未对无穷值取 log、求梯度或作 \(\infty-\infty\)。075、029、076关于这些集合的处理与新公式兼容。

若异质截时闭叶最终覆盖正面 \(F=\{(1,y):y\in[0,1]^4\}\)，每叶各自通过同一个 g，则逐点 \(G\ge H_{\epsilon_i}\ge g\)。闭边、对称边和奇异点可任取包含它的叶，无需连续选择 epsilon 或跨叶求导。符号与置换对称性扩展到 cube 边界。半群共同时间缩放给
\[
 G(\lambda x)=\lambda^{-4}G(x),\qquad
 H_\epsilon(\lambda x)=\lambda^{-4}H_{\epsilon/\lambda}(x).
\]
只用 G 的齐次性延至单位 cube 内，原点 G=∞。固定 epsilon 的 H 不能当作齐次函数；本式不作这种假设。

核对029及原始有限源传递，真正消费的只是同模型 essential cube 下界 g 和完整正面
\(I=\int_{[0,1]^4}G(1,y)^{5/4}dy\)。它们不含统一 epsilon 要求。为明确本新模块未破坏全尺度/一般输入接口，传递链为：

1. 同一一阶半群的 Green 弱 L5/4 和局部 L1 给 \(W=(G-h)_+\in L^1\)，\(0<h<g\)。正部半群不等式给全空间有限正测度 \(\mu=LW+\delta_0\)，质量1；不剔除 S2。
2. \(N=32(g-h),K=W/N\) 满足 \(K\ge1_{[-1,1]^5}/32\) a.e.。一阶生成元和共同半径 dilation 满足
   \(K_r*Lu=(\mu_r*u-u)/(Nr)\)。这一步 r 是所有坐标共用的尺度。
3. 对有界紧支撑非负输入，029的障碍构造给同一个 \(u,\sigma,D\)，\(Lu=\sigma-f\)、\(0\le\sigma\le\kappa\)、\(\kappa|D|\le\|f\|_1\)。取 \(\kappa=\lambda/\|K\|_1\)，正测度卷积使 \(K_r*f\le\lambda\) 在同一 D 外对每个固定 r 成立。先合并正有理尺度的可数零集，再用 K 的 L1 缩放连续性覆盖全部正 r。
4. 对一般 L1，\(f_j=\min(|f|,j)1_{[-j,j]^5}\uparrow|f|\)；正卷积 MCT 和 \(\sup_r\sup_j=\sup_j\sup_r\) 给一般输入及全部尺度的同一弱型界。没有改成独立坐标尺度、非中心 cubes 或非 Lebesgue 背景。
5. G 的 \(-4\) 齐次性及10个面、16个正象限给
   \(\|W\|_1=128Ih^{-1/4}\)，所以
   \(\|K\|_1=4Ih^{-1/4}/(g-h)\)。唯一最优 \(h=g/5\)，得到
   \[C_5\le5^{5/4}I/g^{5/4},\qquad \operatorname{cost}^4=3125I^4/g^5.\]

这条链中没有任何步骤消费 H64/H128 的梯度或 q，也没有要求 G 在 S2 有限。完整 I 必须继续使用完整 G 的上端，不能以 H128 面积分代替。故若后续实际 g110 全闭，新曲率证可合法进入029；在闭合之前它仍只是局部认证工具。

## 8. 本机实际有限核验、证明状态与后续验收边界

登记先于运行，登记 SHA256 为
`01ead6f45f6c212c55cc1976d55668863cda9efaad1e7aaa5b11ea16d9f1c41d`。
唯一新符号脚本 SHA256 为
`6ae552014ffff811b14e627034243fdcd969f360a77c21ada8048a7df98d421a`。
实际执行本机 Python3.15 / SymPy1.14，PID **35727**，exit **0**，脚本记录计算时间 **0.02243462501792237 秒**。四组检查为：两截时系数、端点方向残差分解、正性/严格改善残差、多截时阶梯的三带伸缩身份。没有旧075矩阵/基线重放，没有 H、Green、梯度、log 或超越函数求值，没有使用共享 MATLAB，没有启动后台作业。

有限符号检查只为准确代数身份提供实际计算证据，不是测度论形式化证明；(A)的严格证明是第2—4节，(7)—(11)的严格证明是第5—6节。没有发现导致所提新公式失效的反例或理论缺口。

后续生产仍须逐对象核验：实际 H128 上端与原 jet 同模型/同 epsilon/同 anchor；真实 H64 整叶包含 c 且源锁一致；q 向上端点；实际闭角点余量；完整根树与所有 pending/inflight/fault 的闭覆盖。这些不是本报告已完成的全量数值核验。尤其不能从本报告的理论通过，预称6237旧 pending 已通过、g110全闭、联合成本已取得或 C5 已改善。任务通知给出的当前完整 I 和条件成本未在此重算，正式界保持原状。

最小白名单为本报告、`registration.json`、`check_new_identities.py`、`results.json`、`manifest.json`。结果入口拒绝覆盖已有运行收据；未授权重放旧基线。复核代数可读取脚本和结果；若要新执行应另登记输出文件，避免覆盖冻结收据。
