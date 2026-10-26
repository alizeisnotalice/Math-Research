# B062：高端正事件完整时间预算独立审查 PASS

**PASS_INDEPENDENT_HIGH_FULL_TIME_EVENT_BUDGET。** 对原 body32 同时间 Cauchy 模型、rho=1/64、q=6/5，以及全部 31 个 unsigned tube 的 mean_B>3/4 高端 union 或其任意可测子集 H，独立数学审查和实际纯有理计算接受

\[
 \|G\|_{L^{6/5}(H)}
 <\frac{27779670841319}{10^{18}}
 =0.000027779670841319.
\]

支付范围包括两个真实参数事件和全部 t>0；允许 tubes 交叠、cap 交叠及随后扣除既有空间单元。这里只接受指定 high 全时间预算，不接受全 face I、global g 或新 C_n 结论。

最终作者源锁到达前仅作只读定理审查，无新 scientific 运行。最终 manifest/report 及 26 个白名单/来源条目通过字节核验后，先登记 8 个独立配置，再冻结独立程序并执行一次。没有导入或运行作者 evaluate/check_saved，没有重算旧 M、cap、055/056/060 或 full kernel，没有新方向/平面/triple 全枚举。

## 1. 实际事件及完整余下几何

设 v0=ones，v 为原 body32 的 unoriented 代表 (1,s2,…,s6)。高端 tube 的正交分解为 x_A=alpha·1_A+e，x_B=u·1_B+z，e,z 各自零和，|e|²+|z|²<rho²，u>3/4。由 x1=1，alpha=1−e1；由 A 内坐标≤1，e1≥ei。因此

\[
 x_i\ge1-\sqrt2\rho>47/64\quad(i\in A),\qquad
 x_i>3/4-\rho=47/64\quad(i\in B).
\]

17²>2 即核定前一个严格比较。singleton block 的残差为零，结论仍成立。固定新事件阈值 b=1/4 后，对 |a|<b 的残点 y=x−a v0，有全部 yi≥d=31/64>0。

**所有余下 lines 的统一 gap。** 写 y=d v0+delta，delta_i≥0。若 v 有 h 个负号，1≤h≤5，则 s=v·v0=6−2h，且

\[
 \operatorname{dist}(y,\mathbb Rv)^2
 =\tfrac23h(6-h)d^2
 +2d\sum_i(1-sv_i/6)\delta_i
 +\|P_{v^\perp}\delta\|^2
 \ge\tfrac{10}3d^2.
\]

中间系数全部非负；这证明的是整个正坐标区域，而非只在 d·ones 处算距离或只检查某些翻号类型。也有 |y|²≥6d²。

**所有余下 pair planes 的统一 gap。** 对不同 body lines v,w，把坐标分成 agreement 和 disagreement 两个非空支持块；平面是两个带符号支持向量的正交直和。若它含另一 body 向量 av+bw，则两个块给 |a+b|=|a−b|=1；相减得到 4ab=0，故 body line 必为原来两条之一。这也证明任意三个不同 body lines 独立。

若两个支持块上的 signs 都常数，则该平面包含 ones，与 v,w 都不为 ones 矛盾。所以至少一块存在相反符号的两个坐标 i,j。对平面中的任一点，它在这两个坐标取 c 与 −c（整体换号不影响），故残差至少

\[
 \inf_c\{(y_i-c)^2+(y_j+c)^2\}
 =\tfrac12(y_i+y_j)^2\ge2d^2.
\]

其余坐标只加非负残差。因此 dist(y,span(v,w))²≥2d²，覆盖所有 mixed-sign plane，不假设 pair plane 是 unsigned、不只看与原 tube 共分区的 plane、不漏 mixed 两块或交叠层。无需逐个重枚举 465 个 planes。

## 2. 两个真实参数事件及三次删参的合法顺序

令 r_t 为删去 v0 后 31 个参数之和的六维密度，f_t(a)=t/[pi(a²+t²)]。严格按实际参数事件拆分

\[
 p_t=P_t+B_t,\quad
 P_t=\int_{|a|\ge b}f_t(a)r_t(x-av_0)\,da,\quad
 B_t=\int_{|a|<b}f_t(a)r_t(x-av_0)\,da.
\]

边界参数概率为零，两个非负事件互斥且穷尽。其他 31 个参数始终无限制，因此多个大参数、负大参数和抵消事件都已包含；这不是“只有一个大跳”的假设。

已独审 A041/F041 的逐点 coarea 给直线积分=q5/sqrt6，且删掉 v0 不改变 v0 正交商。其逐点性依赖一删参后的可积直线尾及连续性，不能仅从 r_t 连续就跳步。沿用锁定的 Q5/R1 后

\[
 P_t\le \frac{D t^{-4}}{b^2+t^2},\quad
 D=\frac{R_1}{32\pi^6},\quad R_1=\frac{67320683}{3251298631680}.
\]

这里只有一个选定参数，没有再乘 32。

在 B 的 31 列余项上，在已删 k=0,1,2 列的商 W 中，余下 projected frame 的 Gram 是

\[
 32I_W-(P_Wv_0)(P_Wv_0)^T\le32I_W,
\]

不是 tight-frame 等式。下界≥26I_W 保证所用商密度连续；上界32给正删参乘积式。N=31−k 依次给时间系数 992、960、928；前节几何保证这些商中被删列非零，三次 coarea 的乘积 Jacobian 是 1/sqrt(detGram(triple))。

关键顺序为：先在 |a|<b 内用前节三个距离把分母替换成不依赖 a 的下界；把分母移出积分后，仍保留真实 f_t(a)，最后才把非负 a 积分扩大至全实轴。这恢复原 body32 的三维商密度及已审 Q3，而不是使用新的 omit-ones 商峰。若 v0 落在被删 triple span 中，投影为零，积分 f_t=1，恢复公式仍成立。扩张后的 a 无须继续满足正坐标 gap，因为 gap 分母已经先替换。

继承 A041/F041 已审完整 26,970 个余下 ordered triples 的固定和，不重枚举或重算商峰：

\[
 B_t(x)\le\frac{C_{\rm rem}}
 {(6d^2+992t^2)((10/3)d^2+960t^2)(2d^2+928t^2)},\qquad
 C_{\rm rem}=\frac{44743079}{2970\pi^6}.
\]

因此 B 是完整补事件且每个 t>0 都受控，未把未付事件藏在误差项里。

## 3. 真 face 体积、交叠及 orthant 因子

对 k=|A|、m=6−k，A 的 face deficits 到零和 e 的 Gram 为 I−11ᵀ/k，沿 ones 的特征值为1/k，其余为1，故 inverse face Jacobian 为 sqrt(k)。B 的 mean u 不是单位正交坐标，dx_B=sqrt(m)du dz。于是 dy=sqrt(km)du de dz，正规空间总维度4。

A 交换对称性将径向球体积分分成 k 个最大坐标 cone，除零测 ties 外每个占1/k；这仅是径向体积的 cone 比例，不能错误换成 orthant 因子2^{−k+1}。k=1、m=1 均由零维残差自然包括。

保留 u∈(3/4,1] 长度1/4 和 maxcone，放宽其他 box clipping，只会增大非负积分。四维球体积为 pi²rho⁴/2，所以单管上体积为(pi²/8)rho⁴sqrt(m/k)。tube 类型权重 binom(5,k−1)=(1,5,10,10,5) 合计31。正 union 上界因此是

\[
 |H|\le V=\frac{\pi^2}{8}\rho^4
 \sum_{k=1}^5\binom5{k-1}\sqrt{\frac{6-k}{k}}.
\]

这不是 clipped H 的精确体积，也不要求 tubes 互不交叠。下述全部空间范数都只用这一个全 union 的 V，不向每条 tube 重复支付整个 P 分支。

**源文“positive face slice”的措辞澄清。** 完整密度 p_t 在五个自由坐标上分别反射对称，第一坐标边缘为 Cauchy(32t)，所以准确等式是

\[
 \int_{[0,\infty)^5}p_t(1,y)\,dy
 =\frac1{32}\frac{32t}{\pi(1+1024t^2)}
 =\frac{t}{\pi(1+1024t^2)}.
\]

有界 face [0,1]^5 及 H 只需使用不等式。因 P_t≤p_t，得到 ∫_H P_t≤t/[pi(1+1024t²)]≤t/pi。P_t 自身没有被假定反射对称，不能对它直接乘1/32。作者式(7)的上界正确；此说明只精确化其 slice 用语，冻结作者文件不改。

## 4. 所有正时间一次支付

令 C=D/b²。q=6/5 的 L1–L∞ 插值给

\[
 \|P_t\|_q\le(Ct^{-4})^{1/6}(t/\pi)^{5/6}
 =C^{1/6}\pi^{-5/6}t^{1/6}.
\]

所以0到epsilon的预算为 S=(6/7)C^{1/6}pi^{−5/6}epsilon^{7/6}。对其余所有时间保留 b²+t²，而不是把 P 的长尾丢掉：

\[
 \int_\epsilon^\infty\frac{t^{-4}}{b^2+t^2}dt
 =\frac1{3b^2\epsilon^3}-\frac1{b^4\epsilon}
   +\frac{\arctan(b/\epsilon)}{b^5}=:T_b(\epsilon).
\]

一般 b>0 的导数恒等式经乘正分母后成为多项式恒等式，且该式在无穷远趋零。完整 P 尾预算是 V^{5/6}D T_b；arctan 用正分支。

只选固定 epsilon=(C pi V)^{1/5}。它来自松尾 C/(3epsilon³) 的预算：导数符号由 epsilon^{25/6}−(C pi V)^{5/6} 决定，唯一变号，且两端松预算发散。因此这是唯一松式正极小点，不是准确尾的优化宣称、扫描或新的物理参数搜索。

B 令 z=(6,10/3,2)，A=(992,960,928)，r_i=sqrt(z_i/A_i)。独立重算先换元 t=d·s，把其完整积分写成

\[
 J_B=d^{-5}\frac{\pi(r_0+r_1+r_2)}
 {2A_0A_1A_2r_0r_1r_2(r_0+r_1)(r_0+r_2)(r_1+r_2)}.
\]

三正二次因子公式由两因子积分 pi/[2ac(a+c)] 作一次 partial fraction 得到；关键分子 b(b+c)−a(a+c)=(b−a)(a+b+c)。正重复根由局部共同可积支配连续延拓，固定实现没有小根差除法。于是完整 B 预算为 V^{5/6}Crem J_B，确实包括0到无穷。

最终总预算是 S+V^{5/6}D T_b+V^{5/6}Crem J_B。两个参数事件穷尽，P 的两个时间区间穷尽，B 全时间仅付一次；空间密度无需互不重叠。先在正的有限时间截断上使用 Minkowski，再用单调收敛，得到完整 Green Lq 界。原 body32 在奇异零测集上允许 G=∞，并不要求该 union 的所有点逐点有限；有限 Lq 上界及上述极限论证已经合法处理该情形。

## 5. 独立准确算术及严格端点

本审计使用本机 Python 标准库 Fraction，所有区间运算向外舍入至2^{-192}网格。pi 由 Machin 的100/32项交错级数夹逼，固定小参数 atan 用80项及下一项余项；每个整数根上下端分别严格核对幂不等式，共20个端点根证书。使用整数 Newton root 后再由相邻整数幂夹住根，不使用“根和幂同时对称膨胀”逻辑。

与作者按原分数幂表达式计算不同，独立 P 路径先用 epsilon^5=C pi V 化为

\[
 P_{\rm budget}=\frac1{\pi V^{1/6}}
 \left(\frac{25}{21}\epsilon^2-\frac{\epsilon^4}{b^2}
       +\frac{\epsilon^5}{b^3}\arctan(b/\epsilon)\right).
\]

其中短时项为同一因子乘(6/7)epsilon²，长尾项为其余三项；atan(b/epsilon)=pi/2−atan(epsilon/b)，且独立夹逼证实0<epsilon/b<1/2。B 使用前节 t=d·s 变换及 V^{5/6}=V/V^{1/6}。这重算固定新公式，而不是读取作者小数或只相加其 saved endpoints。

实际独立总上端为

\[
 U_B=\frac{10898488752906381904555848060538100068463422708034091}
 {392318858461667547739736838950479151006397215279002157056}.
\]

其40位小数向外上界为0.0000277796708413170630462732622949266784，严格小于声明值0.000027779670841319。准确差的分子为

```
2898801708094822681483067904009501978476004099530381
```

分母为

```
1496577676626844588240573268701473812127674924007424000000000000000000
```

两者均为正整数。严格比较使用该 Fraction，不由展示小数重建。作者 V、epsilon、Pshort、Plong、Bfull 和 total 的完整 saved dyadic 区间全部包含于本次独立 enclosure 中；这比只核区间相交更强。源数据只用于一致性比较，独立上界本身不依赖作者数值。

声明值也严格小于条件性的全 face 参照4/48750，但这个比较只说明 high 占用的预算较小，不能把尚未支付的空间单元忽略。

## 6. 源锁、实际运行和交付边界

作者最终 manifest SHA-256 为 `d2b25bf37284a34054df29065ca30fe3c04b12116cac7e676304c274b10d2ecf`；report 为 `ad221142eb0e0dc6c99561ea73444ccc1d53dad4c32b4f47931de640e12298ec`；theorem 为 `1d38f4b32dafb8935a2ab3905e16498d63ad5f224f4b6dd1e5b9aa57322c13d9`。固定 A041/F041 及054依赖来源在作者manifest中逐项锁定，运行前后均复核。原 preparation 的“尚未运行”字样保留其历史阶段含义。

独立8配置依次覆盖源锁/回执、全余下几何、真实face体积、全时间代数、pi/root/唯一V与epsilon、完整P、完整B、严格total与端点一致性。一次实际子进程 PID65574，工具 chunk9c7235，exit0，无stderr；程序内耗时0.04807812499348074秒，外层guard墙钟0.07925329098361544秒。guard实际调用 communicate(timeout=180)，超时路径只处理该自有进程组并 wait/reap，无自动重试。本次正常结束且已reap，无活句柄。

results SHA-256 为 `30abeab6ac9a33326a78bc714e63f973c32b4724f3aba44ad6dee8e0ea6d96bb`，由 execution_receipt绑定；input_manifest绑定实际执行源码。review.py和run_guard.py都拒绝覆盖既有运行回执。后续 `python3 verify_seal.py` 仅核封存哈希，不重复8配置或旧计算。

**无阻断性缺口。** 主代理可将指定 high 全时间项计入已支付账本；最终空间分区拼合、middle长时、far、global g、whole I及C_n仍由主流程分别证明。没有修改任何冻结作者源或旧批次，没有stage、commit、push、递归代理或MATLAB。
