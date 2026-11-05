# R055 额外敌意审查：三个指定连接点

2026-11-05。审查人：Astra。结论：**本次指定范围 PROOF-VALID，无必要数学修补。**

实际读取并以 SHA 绑定的源：

* scratch/r055-eighth-complement/whole_height_complement.md
* SHA256 **0684d01f3665a5a064b4786b6eb1d78ea5780ccdae64f96df3ad01647224a1a1**

本审查独立重构下列三个关键点，未先阅读另一位全稿审查人的报告。它不重审 R053 正常积分、不识别八阶系数，也不代替另行进行的全稿审查。R049/R050/R051/R052 的已审同层供货作为依赖使用；另实际读取 R036 adaptive soft-slot 源与 R055 成功的 clip 算核。未修改作者源、旧冻结稿或仓库。

## 1. Whole-height 的上截断：三近面图确已覆盖

写一张 global-min 图的 base 为
\[
 \Omega=\prod_{\ell\ne q}[0,b_\ell],\qquad
 b_\ell=w'_{i\ell}/2\asymp vR_i,\qquad
 H(y)=\min\{b_q,\min_{\ell\ne q}(b_\ell-y_\ell)\}.
\]
固定一个非正常目标，其临界轴对为 \(J=\{k,l\}\)。它的完整名义项 \(g_e(y)\ge0\) 只依赖 \(J\setminus\{q\}\) 的底变量，正支撑迫使这些变量不超过 \(C_Tv^3R_i\)。在其正支撑上，
\[
 b_q\ge c vR_i,\qquad b_j-y_j\ge c vR_i/2
                    \quad(j\in J\setminus\{q\})
\]
对固定 \(T\)、充分小 \(v\) 成立。全 packet 高度的统一上界是
\(\tau=C_Tv^4R_i\)。因此上截断事件只能由
\[
 F=\{1,\ldots,n\}\setminus(J\cup\{q\})
\]
中的某个自由变量落在 \(b_\ell-\tau<y_\ell\le b_\ell\) 造成。

由于 \(\Omega\) 是精确乘积域、\(g_e\) 独立于这些自由变量，逐目标 Fubini 直接给
\[
 \int_\Omega g_e\,\mathbf1_{\{H<\tau\}}
 \le\left(\sum_{\ell\in F}\frac{\tau}{b_\ell}\right)
                \int_\Omega g_e
 \le C_Tv^3\int_\Omega g_e.
 \tag{A1}
\]
对所有图及目标相加便是作者 (13)。这里使用的是完整非正常名义成本，不是 gross 交量，也不需要正常与非正常层的支撑互不相交。

特别核对最容易遗漏的情形：若 \(q\notin J\)，两个临界底变量仍受上述短支撑约束。**在 \(n=3\) 时 \(F=\varnothing\)，这整张三近面图完全没有上截断损失。** 在 \(n\ge4\) 时只出现剩余自由变量的中分面薄条。该论证没有删除三近面主层，更没有把其八阶成本算作误差。共同高度中其他目标如何同时出现，仅改变 \(\tau\) 的固定 \(C_T\)，不改变 (A1) 的因子化。

此外，稿中 (21) 的最浅轴在 pair 内的积分确实保留完整源有效横宽：对各横向最近面符号求和，把两个半宽恢复为一个全宽。最浅轴在 pair 外的单轴贡献与前者之比为
\[
 \frac{P Q}{(P+Q)w'_{iq}}\le C_Tv^2.
\]
零帽时两项均为零，不存在除以极小正宽造成的统一性缺口。因此 (22) 的相对比较与上面的 clip 论证相容。

## 2. 任意后选固定窗的名义成本：没有 \(T\) 循环

旧实际 \(Y_N^{\rm old}\)、正常 packet、原源及稳定常数 \(C_*\) 均先固定。新窗只添加目标，不改这些对象。由 R051 的同坐标精确包含及必要性，
\[
 Y_N^{\rm old}\subset B_N(\mathscr D)\subset B_A(\mathscr D),\qquad
 |B_A(\mathscr D)\setminus\mathscr D|_1=O_T(av^9S).
\]
所以无需假设未侵蚀的模型层在 \(\mathscr D\) 内，就有
\[
 |B_A\setminus B_N|_1
 \le|\mathscr D\setminus Y_N^{\rm old}|_1
                       +|B_A\setminus\mathscr D|_1.
 \tag{A2}
\]
乘 \(\lambda_v\) 正是作者 (14)，旧 raw 的系数为一。

新退款按固定旧正常层和其 raw 余量分开。前者逐坏边最多 \(C_Tv^7R_e^n\)；后者总量最多 \(C_Tv\mathcal R_*\)。在新窗 \(p,q\ge m>0\) 内，双向正部凸性和真实横宽下界给完整名义成本至少 \(c_m v^6R_e^n\)，即使真实根差被 clip 亦成立。因此作者 (15) 的
\[
 \mathfrak F_B\le C_{T,m}v\,\mathsf M_{T,m}+C_Tv\mathcal R_*
\]
并未预先使用新 \(\mathsf M_{T,m}=O(v^8S)\)。结合精确增量分解及 (A1)，得到
\[
 \frac{\mathsf M_{T,m}}{v^8S}
 \le
 \frac{(1+v)^{n-1}+C_Tv}{1-C_{T,m}v-C_Tv^3}
       \frac{\mathcal R_*}{v^8S}+O_T(v^2).
 \tag{A3}
\]
对任一先固定的 \(T,m\)，分母最终为正、前系数趋于一。因此
\[
 \limsup_{v\to0}\frac{\mathsf M_{T,m}}{v^8S}\le C_*.
 \tag{A4}
\]
增长的 \(C_{T,m}\) 只影响起效阈值与消失项；它不出现在 (A4) 的首常数中。这恰足以允许后续先选任意大但固定的 \(T\)，再选任意小但固定的 \(m\)。稿中没有使用 \(T(v)\) 或 \(m(v)\)；若改成同时变窗，(A3) 将不再供应所需结论，但当前陈述没有作该更改。

## 3. 一次完整 \(J\)、常数一转移和固定槽欠量

共同侵蚀及同一个真实上截断给
\(B_N^-\subset B_A^-\subset\mathscr D\)。因此实际上先有集合恒等
\[
 |\mathscr D\setminus B_N^-|_1
 =|\mathscr D\setminus B_A^-|_1+|B_A^-\setminus B_N^-|_1.
\]
第二项至多为完整坏名义高度的积分；乘 \(v/a\)，得到作者 (10) 的因子
\((1+v)^{-(n-1)}\le1\)。此处不存在两次购买同一 raw 的隐含常数。

模型 \(\mathsf P_N=\sum_N I_e-\lambda_v|Y_{2,N}|_1\) 没有减掉全体 \(J\)。作者先从整源层得到
\[
 \lambda_v|B_A(\mathscr D)|_1
 =\lambda_v|Y_{2,N}|_1+\mathsf M_{T,m}+o_{T,m}(v^8S),
\]
再只对 \(B_A^-\) 使用一次强 CDF，并同时保留其外 raw 与端移的固定正份额。这推出 (19)，再与非正常标量严格折扣合并为 (4)。没有把 normal-only 的完整源成本与另一个完整补集成本相加。

长尾 (25) 此时是依据已成立的 \(\mathcal R_*\le C_*v^8S\) 所得的数值 gross 上界；它不再从 (19) 抽取一份 \(J\)。其 \(C_L/T\) 来自避开固定正常小管的实际源条，常数独立于后选新窗。近零角的两支也直接成立：
\[
 p q\le Tm\,\omega_e,\qquad
 \omega_e=\min(1,C_{\rm soft}\min(p,q)),\quad C_{\rm soft}\ge1/m.
\]
饱和时使用 \(\min(p,q)<m\)；未饱和时使用
\(T\min(p,q)=(T/C_{\rm soft})\omega_e\le Tm\omega_e\)。

R036 的 clean 阈值是每边实际
\(\lambda_e=\min L/(1000r_i)\)，其 probe、helper 与碰撞片不依赖软权常数 \(C_{\rm soft}\)。因此改变软权只改变片上的数值权和误差常数，正常片仍取权一，旧正常占用及 \(D_N\) 均可保持不变。于是
\[
 W_B,\ W_{\mathcal D}\le\mathcal C-W_N+O_{C_{\rm soft}}(v^3S)
                      =D_N+O_{C_{\rm soft}}(v^3S)
\]
可分别用于作者的两项上界；近零 gross 已先放入任意 \(\delta v^8S\)，并未再额外消费一份槽容量。

参数顺序因此合法：先由固定 \(C_*,C_L\) 选 \(T\) 使长尾常数小；再选 \(m\) 使 \(Tm C_*\) 小；随后固定软权，最后令 \(v\to0\)。任何 \(O_{C_{\rm soft}}(v^3S)\) 容量误差乘六阶价格都为 \(o(v^8S)\)。从 (4) 到 (6) 的残余槽项为
\[
 -[q_N-(1-\gamma_\nu)A_6]D_N,
\]
其括号至少为固定正数乘 \(v^6\)。这保留了可能真实达到八阶的负槽损失，而不是把它丢入误差。正常模型和 \(q_N\) 均没有随上述辅助阈值重定义。

## 4. 证据与结论范围

本机已实际读取 r055-computation/clip_increment.wl/log；两项违反集为 False，其作用仅是确认非负高度的 clip 增量标量关系。

另在本次审查中实际用 Python Fraction 核了自由变量数 \(0,1,2,3\) 的精确联合薄条比例
\(1-\prod(1-\tau/b_\ell)\le\sum\tau/b_\ell\)，并核三个 \(M_3/M_2\) 有理边界（等宽、极短一边、零帽），均通过。输入为 \(v=1/20\)、\(\tau=7v^4\)、自由半宽 \(v/2,3v/5,4v/7\)；这些有限检查不替代 (A1) 的乘积域证明。

三项委托重点均未发现实质缺口，源无需补丁。本结论限定于同一实际八阶竞争族、固定 \(n\ge3\)、先固定辅助窗再取 \(v\) 极限的 R055 补集结算。未使用 R054 微孔或内部回收接口，也未审定 R053 的正常八阶系数及最终全局总装。
