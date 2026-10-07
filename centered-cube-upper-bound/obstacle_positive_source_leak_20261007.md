# 原障碍正来源的单坐标面泄漏：域外 trace 与一次来源付款

2026-10-07。新增一般冻结原核子支；完整 ordered 弱端点及实际 cube 余项未闭合。不重复已付负部/截高、不修改总账或旧数值。依据 J03 保留整个原域外与真实跳落点；下文付款是 receiver Lebesgue 费用，不是标签行概率。

## 1. 同一正来源和未付目标

读取 [上轮域外差](ordered_frozen_exterior_difference_20261007.md)、[正预解重装](resolvent_obstacle_repacking_20261007.md) 和根的 [Gamma score 候选](resolvent_score_next_interface_20261007.md)。固定原坐标对称保质量正核 \(G_i=G_{c,i}\)，
\[
 S=\sum_i(I-G_i),\quad Q=\sum_iG_i,\quad
 K_r=\prod_i[(1-r)I+rG_i],\quad H_r=e^{-rS}.
\]
原障碍仍为同一 \(u\ge0\)、\(u\in L^1\cap L^2\)、\(\Omega=\{u>0\}\)，
\[
 Su=\nu_b-\mu_b,\qquad
 0\le\mu_b\le\kappa,\quad\mu_b=\kappa\text{ 于 }\Omega,
 \quad\operatorname{supp}\nu_b\subset\Omega,\quad
 \int\nu_b=\int\mu_b=W_b.
\tag{1}
\]
固定 \(t>0\)，
\[
 R_t=(I+tS)^{-1},\qquad
 h_t=R_tSu=(u-R_tu)/t,\qquad
 f_t=(h_t)_+,\qquad
 N_{r,t}=(K_r-H_r)(I+tS).
\tag{2}
\]
已证 \(f_t\) 支撑于原 \(\Omega\)、\(\int f_t\le W_b\)。本轮只研究
\[
 \kappa\big|\{x\in\Omega^c:\sup_{0\le r\le1}
                  [N_{r,t}^{\rm nl}f_t(x)]_+>\kappa\}\big|.
\tag{3}
\]
纯 holding 项在 \(\Omega^c\) 为0。正来源不是任意可重选的 \(L^1\) 密度；它由原同一 \(u,\nu_b,\mu_b\) 生成。

## 2. 原域外的首跳 trace，未丢来源

由 \(R_tu\ge0\)，
\[
 0\le f_t\le u/t,\qquad
 Qf_t\le Qu/t=\mu_b/t\le\kappa/t
                   \quad\text{于原 }\Omega^c.
\tag{4}
\]
这里只在原域外使用 \(u=0\)；没有把域外 cap 延伸到全部空间。式 (4) 是一个真实核输出的点态 trace，不是独立性假设或按 receiver 重新分配来源。

每个恰好一次坐标跳跃 \(G_i\) 的系数相同，准确为
\[
\begin{split}
 a_{1,t}(r)
 &=(1-r)^{n-1}[(1+t(n+1))r-t]\\
 &\quad+e^{-nr}[t-(1+tn)r]
 =d_1(r)+t\beta_n(r),\\
 d_1(r)&=r[(1-r)^{n-1}-e^{-nr}],\\
 \beta_n(r)&=[(n+1)r-1](1-r)^{n-1}+(1-nr)e^{-nr}.
\end{split}
\tag{5}
\]
这只是总 degree 为1；单坐标面上的 \(G_i^m\)、\(m\ge2\) 仍须另处理，不能把“访问一个坐标”误写成“一次跳跃”。

对所有 \(n\ge2\)、\(0\le r\le1\)，
\[
 |d_1(r)|\le1/n,\qquad |\beta_n(r)|\le18/n.
\tag{6}
\]
第一式因两个非负数 \(r(1-r)^{n-1}\) 与 \(re^{-nr}\) 均不超过 \(1/n\)。第二式的自含证明如下。

若 \(r\le1/2\)，令 \(q=nr\)、\(b=(1-r)^{n-1}\)。用
\[
 0\le-\log(1-r)-r\le \frac{r^2}{2(1-r)}\le r^2
\]
分开比较 \(b\) 与 \(e^{-(n-1)r}\)，再与 \(e^{-nr}\)，得
\[
 |b-e^{-nr}|\le(r+nr^2)e^{-(n-1)r}
                         \le(r+nr^2)e^{-q/2}.
\]
因此
\[
 |\beta_n|\le\frac{q^3+2q^2+2q}{n}e^{-q/2}
 \le\frac{216/e^3+32/e^2+4/e}{n}<\frac{18}{n}.
\tag{7}
\]
用 \(e>8/3\) 即可有理地核最后常数：右侧分子不超过 \(1113/64<18\)。若 \(r\ge1/2\)，
\[
 |\beta_n|\le(n+2)2^{-(n-1)}+(n+1)e^{-n/2}\le8/n.
\]
其中 \(n(n+2)2^{-(n-1)}\le4\) 从 \(n=2\) 开始递减；
\(n(n+1)e^{-n/2}\le(3/2)n^2e^{-n/2}\le24/e^2<4\)。
所以端点也覆盖，无除以 \(r\) 的失效公式。

由 (4)–(6)，整个一次跳跃项 \(B_{1,r,t}=a_{1,t}(r)Qf_t\) 满足
\[
 \sup_r|B_{1,r,t}|\le\gamma_{n,t}\kappa
       \quad\text{于原 }\Omega^c,\qquad
 \gamma_{n,t}=\frac{18+1/t}{n}.
\tag{8}
\]
它无需额外 Lebesgue 费用。\(t\ge1,n\ge512\) 时 \(\gamma_{n,t}<1/16\)。这是真正从原障碍结构得到的空间删项。

## 3. 单坐标两跳项严格非正

每个 \(G_i^2\) 的系数为
\[
 a_{2,t}(r)
 =tr[e^{-nr}-(1-r)^{n-1}]
            -\frac{(1+tn)e^{-nr}r^2}{2}\le0
 \quad(n\ge2,t>0).
\tag{9}
\]
证明在 \(r=0\) 是等号。若 \((n+1/t)r\ge2\)，写成
\[
 a_{2,t}=tr e^{-nr}[1-(n+1/t)r/2]-tr(1-r)^{n-1}
\]
即可。否则 \(nr<2\)。对 \(n=2\)，
\(e^{2r}(1-r)\ge1-r=1-nr/2\)。
对 \(n\ge3\)，\(r<2/n\le(n+2)/(2n+1)\)，上面的 log 下界给
\[
 nr+(n-1)\log(1-r)
 \ge r-\frac{(n-1)r^2}{2(1-r)}
 \ge-\frac{nr}{2}.
\]
故
\[
 e^{nr}(1-r)^{n-1}\ge e^{-nr/2}
       \ge1-nr/2\ge1-(n+1/t)r/2,
\]
亦得 (9)。因此 \(a_{2,t}G_i^2 f_t\) 可从正上界删除，不收来源。这里保留原 convolution power \(G_i^2\)，没有将它换成 \(G_i\)。

## 4. 所有单坐标重复尾有一份强 L1 预算

对 \(m\ge3\)，\((I+tS)K_r\) 没有 support 为单坐标、degree 为 \(m\) 的项。完整差的该项系数是
\[
 a_{m,t}(r)=e^{-nr}\frac{r^{m-1}}{m!}
                     [tm-(1+tn)r],
\]
所以
\[
 [a_{m,t}(r)]_+
 \le t e^{-nr}\frac{r^{m-1}}{(m-1)!}
 \le\frac{t}{n^{m-1}}\qquad(0\le r\le1).
\tag{10}
\]
第二个界在整个 \(r\ge0\) 上最大化；令 \(k=m-1\)，最大值为
\(t(k/e)^k/(n^k k!)\le t/n^k\)，最后由正级数
\(e^k\ge k^k/k!\) 得到，不依赖 Stirling 渐近。

定义同一固定正来源的包络
\[
 P_tf_t=\sum_{i=1}^n\sum_{m\ge3}\frac{t}{n^{m-1}}G_i^mf_t.
\]
所有项非负，Tonelli 和原保质量给
\[
 \boxed{\|P_tf_t\|_1
       =\frac{t}{n-1}\|f_t\|_1\le\frac{t}{n-1}W_b.}
\tag{11}
\]
没有按 \(r\)、坐标或重复跳次数重新领一份 \(W_b\)。这是全连续参数的包络，单参数 TV 界没有被误当 maximal 界。合并 (8)–(11)，**全部单坐标面**的正泄漏，在原 \(\Omega^c\) 上被
\(\gamma_{n,t}\kappa+P_tf_t\) 支配。

## 5. 精确新补集与参数选择

按真实坐标 support 分解
\[
 N_{r,t}^{\rm nl}=N_{r,t}^{[1]}+N_{r,t}^{[\ge2]}.
\]
后者是至少两个坐标被访问的面，包含重复跳与 full-dimensional 项。奇异面按其独立测度分解，不混为全空间连续密度。各 \(r\) 的级数在 signed-kernel TV 中绝对收敛；单面 positive envelope (11) 亦给共同可测版本。

对于 \(0<a<1-\gamma_{n,t}\)，
\[
\begin{split}
 \kappa|\{x\in\Omega^c:\sup_r[N_{r,t}^{\rm nl}f_t]_+>\kappa\}|
 &\le\frac{t}{a(n-1)}W_b\\
 &\quad+\kappa|\{x\in\Omega^c:
 \sup_r[N_{r,t}^{[\ge2]}f_t]_+>
                  (1-\gamma_{n,t}-a)\kappa\}|.
\end{split}
\tag{12}
\]
取 \(n\ge512,t\ge1,a=1/2\)，第一项为 \(2tW_b/(n-1)\)，剩余阈值至少 \(7\kappa/16\)。可保留更准确的 \((1/2-\gamma_{n,t})\kappa\)。

固定 \(t=1\) 已有常数费用；取预定 \(t=O(\log n)\) 则 (11) 费用为 \(O(\log n/n)W_b\)，且已有封顶负来源费为 \(40(1+t)^2W_b\)。本证明不显示增大 \(t\) 能改善多坐标面，不能把其 smoothing 当成新端点。same-\(u,\Omega\) 重装来源 \(\nu'_t\) 的质量仍至多 \(2W_b\)，但未被多次使用。

## 6. 多面仍未付，Gamma score 的准确作用范围

根候选的 Gamma/Exp 潜变量推导可核对：对
\(p_r=(1-r)\delta+rG\)、\(q_r=G*p_r\)，标记空间 RN 密度在 holding 为0，在 Exp 标签为 \((1-r)/r+T\)。均值为1、方差为 \(1/r+r-1\)；条件 Jensen 给实际卷积核的同向 \(\chi^2\) 上界。坐标张量 latent score 及 Poisson 总跳数分别给
\[
 \|SK_r\|_{\rm TV}\le\sqrt{n(1/r+r-1)},\qquad
 \|SH_r\|_{\rm TV}\le\sqrt{n/r}.
\]
这些是固定参数核的结论，未控制连续 sup，更未控制同一域外多面差。尤其分别取三角还丢失 \(K-H\) 的 cancellation，不能用来支付 (12) 的剩余项。

\(t\asymp\log n\) 也不使 \(R_t\) 仅剩 allvisited：它是指数混合时钟的原 \(H\)；在时钟 \(q\le(\log n)/2\) 时 allvisited 概率至多 \(e^{-\sqrt n}\)。若 \(t\le C\log n\)，proper-face 总质量至少
\[
 (1-e^{-1/(2C)})(1-e^{-\sqrt n}),
\]
仍为常数量级。它是预解核的质量诊断，不是 \(f_t\) 的独立后验，也不是 (12) 的弱反例。

目前对多面只有 (4) 的首跳 cap。作用了其他坐标前缀后 \(Qf_t\) 的原域外不等式不能沿自由再进入继承；单独某个 \(G_i^m\) 的空间密度与 \(G_i\) 的点态比也没有统一上界。若释放总质量 \(O(W_b)\) 的退出来源，再调用原 ordered continuation 最大费用，则循环使用尚未证明的端点。新 \(R_t\) killed 障碍没有把 joint multiplier 变成单生成元正预解函数。

所以本轮确切新增付款是 (12) 的全部单坐标面子支；其阈值有已知有理 margin，余项不是“任意退出源”而仍是原 \(f_t=(R_tSu)_+\)、原 \(\Omega^c\) 上的多坐标面 signed max。尚无 polylog 空间付款，也未移植到 actual FIRST/moving-hard 账。

## 7. 附加工具：短总 degree 的空间 trace（不是诊断 K）

根指出的更一般 trace 同样保持原域。由 \(Su\ge-\kappa\)，全空间
\(Qu\le nu+\kappa\)。正性及 \(Q1=n\) 迭代给
\[
 Q^ju\le n^ju+j n^{j-1}\kappa,\qquad
 Q^jf_t\le j n^{j-1}\kappa/t\quad\text{于 }\Omega^c.
\tag{13}
\]
不把右式按每条路径独立复制。在这里 \(j\) 是核展开的真实总跳次数，不是原 actual 门的辅助 Bernstein 诊断 \(K\)，也不是 early continuation 活动数。

记 \(N_{r,t}^{(j)}\) 为总 degree 恰好 \(j\) 的核，包括所有 support。令 \(e_j(G)=\sum_{|A|=j}\prod_{i\in A}G_i\)。正算子有
\[
 e_j(G)\le Q^j/j!,\qquad Qe_{j-1}(G)\le Q^j/(j-1)!.
\]
精确 degree-\(j\) 项为
\[
\begin{split}
 (I+tS)K_r\big|_j
 &=(1+tn)r^j(1-r)^{n-j}e_j(G)\\
 &\quad-tr^{j-1}(1-r)^{n-j+1}Qe_{j-1}(G),\\
 (I+tS)H_r\big|_j
 &=e^{-nr}\left[\frac{(1+tn)r^j}{j!}
               -\frac{tr^{j-1}}{(j-1)!}\right]Q^j.
\end{split}
\]
对 \(t\ge1,1\le j\le n/2\)，用 (13) 及
\(\sup_{q\ge0}e^{-q}q^k/k!\le1\)，H 项绝对值至多 \(3j\kappa\)。
K 项使用
\[
 \sup_{0\le r\le1}
 \frac{(nr)^k}{k!}(1-r)^{n-j}
 \le(1-j/n)^{-k}\le e^{2j^2/n}
 \quad(0\le k\le j)
\]
及第二项更大的 \(n-j+1\) 指数，得到
\[
 \sup_r|N_{r,t}^{(j)}f_t(x)|
       \le6j e^{2j^2/n}\kappa
                      \quad\text{于原 }\Omega^c.
\tag{14}
\]
因此预设 \(1\le J\le\sqrt n/2\)，整片短 degree 满足
\[
 \sup_r\left|\sum_{j=1}^J N_{r,t}^{(j)}f_t\right|
       \le6J(J+1)\kappa.
\tag{15}
\]
这里 \(e^{2J^2/n}\le e^{1/2}<2\)，常数覆盖求和。holding 在原域外已消失。

(15) 是真正空间核输出界，但在当前固定阈值 \(\kappa\) 下不能宣称该整片自动低于阈值。若在固定障碍入口**预先**用 \(\kappa=\lambda/[48J(J+1)]\)，短片至多 \(\lambda/8\)，且原活跃域的 \(\lambda|\Omega|\) 费用至多 \(48J(J+1)W_b\)。取 \(J\asymp\log n\) 是 polylog 费用。它给可选的短 degree 空间付款接口；剩余总 degree \(>J\) 的同一正来源 signed max 未付，不是已闭合方案。没有将多个 J 的原来源预算相加，亦未因 \(j\) 大就免费假定 allvisited。

## 8. 新证书范围

先保存 [注册](obstacle_positive_source_leak_registration_20261007.json)，随后一次执行 [新脚本](obstacle_positive_source_leak_exact_guard_20261007.py)，[结果](obstacle_positive_source_leak_results_20261007.json) 已终态 PASS。三轮 \(n=8,32,128\)，3003 项精确 Fraction/有理区间断言，无失败；各轮 \(n\) 同时核 \(t=1\) 和预设有理 \(t=\lceil\log_2(n+2)\rceil\)。脚本 SHA256 为 2b31f212e6624e73059b12a58a603e0720571a30c653a3a2cc4ee225714b5ff9。无 live handle。

Exp 用正 Taylor 下和及首遗漏项几何尾上界，\(e^{-q}\) 由 reciprocal 反向包住；没有浮点 exp、积分或拟合。覆盖 (6)/(9)/(10) 的新系数符号与包络、单面无穷几何尾、短 degree 参数及原全空间/域外 \(Q^j\) 代数组件。后者是有限正 Markov 模型，明确不是原 \(\mathbb R^n\) 核样本；系数断言则适用于原任意 \(c>0\) 的核展开，因为其结论只使用相同系数和原正保质量核。

原谱、旧负密度及 Rt 障碍守卫均不重跑。不求解新特殊输入、不将 coefficient 检查称 actual FIRST 样本。一般连续参数命题由正文证明承担，有限证书不认证多面弱端点或 actual geom 费用。
