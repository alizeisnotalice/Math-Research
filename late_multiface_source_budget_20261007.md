# 晚参数多面：同源 signed 回并与重复跳扣项

2026-10-07。读取 current_joint_budget 最后合同与单面稿；只编辑本稿和同前缀新证书。当前实际 cube 主账仍为 \(R_{\rm angle}\)，没有原 FIRST/moving-hard 的新付款。冻结问题本轮得到严格 signed 回并、全单面绝对包络及 trace 交换障碍；晚参数多面 polylog 弱费用仍未证明。

## 1. 当前目标与来源不变

固定原 \(G_i=G_{c,i}\)、\(S=nI-Q\)、\(Q=\sum_iG_i\)，
\[
 K_r=\prod_i[(1-r)I+rG_i],\quad H_r=e^{-rS},\quad D_r=K_r-H_r.
\]
原障碍 \(u\in L^1\cap L^2\)、\(u\ge0\)、\(\Omega=\{u>0\}\)，
\[
 \sigma=Su=\nu_b-\mu_b,\quad
 \operatorname{supp}\nu_b\subset\Omega,\quad
 0\le\mu_b\le\kappa,\quad\mu_b=\kappa\text{ 于 }\Omega,
 \quad\int\nu_b=\int\mu_b=W_b.
\tag{1}
\]
固定 \(R_t=(I+tS)^{-1}\)、\(h_t=R_t\sigma\)、\(f_t=(h_t)_+\)，
\(N_{r,t}=D_r(I+tS)\)。没有重新选择来源或 obstacle。

当前总账采用 \(t=1\)、\(L=\lceil\log_2(n+2)\rceil\)，已付 \(r\le L/n\) 的差核强 maximal 与晚参数正单面片。剩余是
\[
 \kappa\big|\{x\in\Omega^c:
 \sup_{r\in[L/n,1]}[N_{r,1}^{[\ge2]}f_1(x)]_+
                           >7\kappa/16\}\big|.
\tag{2}
\]
support \(\ge2\) 是真实核访问坐标的集合，既包含 squarefree，也包含重复跳；不是 actual 诊断 Bernstein \(K\) 或 early 活动数。若 (2) 有 polylog \(W_b\) 付款可接已显示冻结条件账；它尚不能接实际 \(R_{\rm angle}\)。

### 本轮首选替代：直接 D 正源，不额外领 N score 费

根提供的更简洁路线可独立核对：由 \(0\le\mu_b\le\kappa\)，
\[
 |D_r\mu_b|\le\kappa,\qquad
 \{\Omega^c:\sup_r[D_r\sigma]_+>2\kappa\}
 \subset\{\Omega^c:\sup_r[D_r\nu_b]_+>\kappa\}.
\tag{1a}
\]
因此直接以原 \(\nu_b\) 为入口可省略 N 负源/overlap 的费用；下列 N 重组是替代审计，不能与此新路线累加账。

D 的单面项为
\[
 D_r^{[1]}
 =d_1(r)\sum_iG_i
       -e^{-nr}\sum_i\sum_{m\ge2}\frac{r^m}{m!}G_i^m,
\quad d_1=r[(1-r)^{n-1}-e^{-nr}].
\]
重复项全部非正。可证明全连续 \(r\) 上
\[
 [d_1(r)]_+\le14/n^2\quad(n\ge2).
\tag{1b}
\]
在 \(r\le1/2\)，沿单面稿的比较，
\[
 |d_1|\le(q^2+q^3)e^{-q/2}/n^2
       \le(16/e^2+216/e^3)/n^2<14/n^2,
\quad q=nr.
\]
用 \(e>8/3\) 即核 \(873/64<14\)。在 \(r\ge1/2\)，
\([d_1]_+\le r(1-r)^{n-1}\le2^{-n}\le2/n^2\)。
于是固定正包络 \(P\nu_b=(14/n^2)Q\nu_b\) 满足
\[
 \|P\nu_b\|_1\le14W_b/n.
\tag{1c}
\]
域外 holding 项为0（对正源其 coefficient 亦非正），因此同一原源正 difference 可仅留下 \(D^{[\ge2]}\)；全重复多面仍保留 signed，未取每面的绝对值。

完整条件回代可直接选原 obstacle 的 \(\kappa=\lambda/4\)：原 \(\Omega\) 付 \(4W_b\)；其外
\[
 K_r\nu\le\kappa+H_r\nu_b+P\nu_b+
                               [D_r^{[\ge2]}\nu_b]_+.
\]
若真正有
\[
 \kappa|\{\Omega^c:\sup_r[D_r^{[\ge2]}\nu_b]_+>\kappa\}|
 \le C_D(n)W_b,
\tag{1d}
\]
则已有 fixed \(H\) 弱费 \(A_{\rm fix}\) 给
\[
 A_{\rm ord}\le4+4A_{\rm fix}+56/n+4C_D(n).
\tag{1e}
\]
这是首选的替代条件合同，\(C_D\) 尚未证明。晚参数可限制在当前未付窗口；初期原 difference 已另有费用。(1e) 的 all-\(r\) 式明确其完整回代常数，不能将新的多面 \(C_D\) 当原 N 的 \(C_{\rm late}\) 使用相同系数，亦不能混入实际 cube 账。

## 2. 同一正来源保留整个核取消

不先按面取正部。定义
\[
 b_t=R_t\mu_b-(h_t)_-
       =\min(R_t\nu_b,R_t\mu_b).
\]
原保常数正核给
\[
 0\le b_t\le\kappa,\qquad \int b_t\le W_b,\qquad
 f_t=R_t\nu_b-b_t.
\tag{3}
\]
因此准确在 \(L^1\cap L^2\) 上
\[
 \boxed{N_{r,t}f_t=D_r\nu_b-N_{r,t}b_t.}
\tag{4}
\]
这是原同一 \(\nu_b,\mu_b,u\) 的等式；\(b_t\) 不是另领一份 arbitrary 来源，更没有把 \(\nu_b\) 改为 \(R_t\nu_b\) 重新归一化。

由已有联合谱 maximal 工具，
\[
 \left\|\sup_r|N_{r,t}b_t|\right\|_2^2
       \le40(1+t)^2\kappa W_b.
\tag{5}
\]
任意 \(\theta>0\) 因而满足
\[
 \kappa|\{\sup_r|N_{r,t}b_t|>\theta\kappa\}|
       \le\frac{40(1+t)^2}{\theta^2}W_b.
\tag{6}
\]
原 \(\Omega^c\) 上 holding 项作用于 \(f_t\) 为0，但它不作用于 (3) 的两项各自为0；(4) 要在全核层先恢复取消。不能把 \(N^{[\ge2]}R_t\) 写成 \(D^{[\ge2]}\)：卷积 \(R_t\) 会新增坐标，改变面标签。

## 3. 全部单坐标面的 signed 包络也只有一份来源费

前一 [单面稿](obstacle_positive_source_leak_20261007.md) 只为正删项采用极小费用。本节需要绝对包络，以允许 (4) 合法回并；不是重复收此前实际交通。

记单面总项 \(N_{r,t}^{[1]}=\sum_i\sum_{m\ge1}a_{m,t}(r)G_i^m\)，\(t\ge1,n\ge2\)。已有精确度1系数给
\[
 \sup_r|a_{1,t}|\le(18t+1)/n.
\]
度2系数虽非正，其绝对值也满足
\[
\begin{split}
 |a_{2,t}|
 &\le tr[e^{-nr}+(1-r)^{n-1}]
       +\frac{(1+tn)e^{-nr}r^2}{2}\\
 &\le\frac{2t}{n}+\frac{2(1+tn)}{n^2}
       \le\frac{5t}{n}.
\end{split}
\tag{7}
\]
这里分别使用 \(\sup re^{-nr}\le1/n\)、\(\sup r(1-r)^{n-1}\le1/n\)、\(\sup r^2e^{-nr}\le4/n^2\) 与 \(1\le tn/2\)。

对 \(m\ge3\)，完整系数包含相同 Poisson 项，故
\[
 \sup_r|a_{m,t}|
 \le \frac{t}{n^{m-1}}+\frac{1+tn}{n^m}
 =\frac{2t+1/n}{n^{m-1}}.
\tag{8}
\]
以固定正原来源 \(f_t\) 及保质量卷积求和，Tonelli 给连续全参数的强 \(L^1\) 包络
\[
 \boxed{
 \left\|\sup_r|N_{r,t}^{[1]}f_t|\right\|_1
 \le C_tW_b,\quad
 C_t=23t+1+\frac{2t+1/n}{n-1}\le26t+1.}
\tag{9}
\]
这是所有单坐标面与全部重复次数，不是一个固定 \(k\) 子族。没有按面重新分配 \(f_t\)，也没有以每个面弱费用相加。\(\sum_{m\ge3}\) 的无穷几何尾显式包含于 (9)。

## 4. 晚参数多面可回并到真实原源差，但仍等价于未付端点

设 \(\mathcal I=[L/n,1]\)。在原 \(\Omega^c\)，
\[
 N_{r,t}^{[\ge2]}f_t
   =D_r\nu_b-N_{r,t}b_t-N_{r,t}^{[1]}f_t.
\tag{10}
\]
因此任何 \(a>0\) 的多面 weak 事件，若按三个阈值各 \(a\kappa/3\) 分割，满足
\[
\begin{split}
 \kappa|\{\Omega^c:\sup_{\mathcal I}
                   [N^{[\ge2]}f_t]_+>a\kappa\}|
 &\le\kappa|\{\Omega^c:\sup_{\mathcal I}
                   [D_r\nu_b]_+>a\kappa/3\}|\\
 &\quad+\left[\frac{360(1+t)^2}{a^2}
                   +\frac{3C_t}{a}\right]W_b.
\end{split}
\tag{11}
\]
它保留原 signed cancellation，且真正支付了封顶 overlap 源与全单面恢复误差。取 \(a=7/16,t=1\) 是确定常数；大 score-TV 障碍不能直接移植到 (11) 的新主项，因为完整 \(N R_t=D\) 已精确取消。

但
\[
 [D_r\nu_b]_+\le K_r\nu_b,\qquad
 K_r\nu_b\le [D_r\nu_b]_++H_r\nu_b,
\]
所以该原 \(\Omega^c\) 上的 difference 弱预算，与 ordered 初值 \(\nu_b\) 的域外弱预算只差已有 fixed \(H\) 费用和阈值常数。它不是新完整解，也不是可将 \(A_{\rm ord}\) 循环代入的退出源接口。原 \(\Omega\) 的体积已由 \(W_b/\kappa\) 支付；\(\nu_b\) 与其同一 obstacle 域相配。

固定 \(r\)，\(\|D_r\|_{\rm TV}\le2\) 给常数费用；receiver 选时仍不是固定参数平均。根的 Gamma \(r^3\) 取消与 log-\(r\) 占用只付此前窗口；其积分不能控制 (11) 的晚参数 sup。新弱预算尚缺此项，不以 (11) 的重装称阶数改善。

## 5. 保留重复跳扣项的 squarefree 递推

令 \(e_k(G)=\sum_{|A|=k}G_A\)，其中 \(G_A=\prod_{i\in A}G_i\)，并定义正算子
\[
 \mathcal R_k=\sum_{|A|=k}\sum_{i\in A}G_i^2G_{A\setminus\{i\}}.
\]
注意 \(\mathcal R_k\) 是算子标签，非新的 source。\(\mathcal R_0=0\)。精确计数为
\[
 Qe_k=(k+1)e_{k+1}+\mathcal R_k,
\]
故原同一势上
\[
 \boxed{(k+1)e_{k+1}u
   =n e_ku-\mathcal R_ku-e_k\nu_b+e_k\mu_b.}
\tag{12}
\]
完整 \((I+tS)K_r\) 的 degree 展开及该递推保留两种负扣项：
\[
\begin{split}
 (I+tS)K_r
 &=\sum_{k=0}^n r^k(1-r)^{n-k}
             [(1+tn)e_k-t\mathcal R_k]\\
 &\quad-t\sum_{k=1}^n
             k r^{k-1}(1-r)^{n-k+1}e_k.
\end{split}
\tag{13}
\]
它不是把每个面取绝对值累加，重复跳与 squarefree 仍来自同一 \(u\)。令 \(u=0\) 于外部仅能在 \(k=0\) 启动 trace；\(e_ku\) 不是新的零外边界势。

若丢掉 \(-\mathcal R_ku\)、\(-e_k\nu_b\)，并使用
\(e_k\mu_b\le\kappa\binom nk\)，只剩原已知的正上界。若用 \(e_k\le Q^k/k!\) 再套 \(Q^ku\le n^ku+k n^{k-1}\kappa\)，without-replacement 的比较价为
\[
 \frac{n^k}{(n)_k}
 =\prod_{\ell=0}^{k-1}(1-\ell/n)^{-1},\qquad
 \log\frac{n^k}{(n)_k}\ge\frac{k(k-1)}{2n}.
\tag{14}
\]
在 \(k\asymp n^{2/3}\) 时是 \(\exp(c n^{1/3})\)。真正晚窗口包含该范围；不能像短 degree 工具那样将此因子当 polylog。缺的是保留 (12) 中原 repeat 扣项和同源 forcing 后的空间 budget，而非再给一个行上限。既不能假设 \(G_i^2u\ge G_iu\)，也不能在各 coordinate prefix 后继承原 \(\Omega^c\) 的 cap。

## 6. R 的均匀 count 标签不等于共同位移来源

根提供的标签身份正确：\(R_1=\int_0^\infty e^{-\tau}H_\tau d\tau\)，原连续跳核使 visited mask 活跃概率 \(p=1-e^{-\tau}\)；\(e^{-\tau}d\tau=dp\)。故每个 mask \(A\) 的质量
\[
 \int_0^1p^{|A|}(1-p)^{n-|A|}dp
 =\frac{1}{(n+1)\binom n{|A|}},
\]
总 visited count 均匀于 \(0,\ldots,n\)。

但是 \(R\) 在给定 mask 内含 truncated-Poisson 重复跳及其 Gamma 时间，\(\int_0^1K_rdr\) 则每轴只有一次；位移核不相等。即使以 \(f_1/R\nu_b\in[0,1]\) 正分配这些来源标签，接受密度也在输出位置上截断，而非与后续卷积交换。不能据此把 \(N R=D\) 的取消放在每个已裁剪来源标签内，或将均匀 count 变成 actual 后验独立。式 (3)–(4) 是避免这个非法交换的完整同源替代。

## 7. 与新加权峰接口的关系及最小缺口

Tensor 正研究的全局 winner 参数 \(p\)、加权峰高度预算属于 \(K_p\nu_b\) 的真实固定家族。可与 (11) 的原来源主项直接接，而不是对 \(f_t\) 或 face 标签另取一个 winner。它的 overshoot 分支若可付，会留下 moderate 高度；本稿没有从 (12) 得到该剩余高度与 mask/count 的额外可积空间关系。

这里 \(p\) 是 receiver 的核参数、\(k\) 是真实核展开 degree 或 visited mask 大小（依具体式），均不等于 actual 诊断 \(K\)。旧所有 FIRST/fullfuture/出生/LCA/history 门没有被认证移植，因此本轮冻结式也不改变实际 \(R_{\rm angle}\)。

可证明的新增源一次预算止于 (9)/(11)。尚缺：

* 直接对原同一 \(\nu_b,\Omega\) 证明 (11) 的晚 difference 空间弱费；或
* 将 (12) 的 repeat 与 source 扣项转为可求和空间/首次新增坐标退出预算，同时保留同一 forcing，禁止先释放来源再引用 \(A_{\rm ord}\)。

只看正 \(Q^j\) trace、均匀 mask count 或 fixed-\(r\) TV 都不提供该预算。本文不构造特殊输入冒充 actual 反例，也不提出未经证明的漂亮强式。

首选直接 D 路线在 (1b)–(1e) 将上述预算进一步简化，不需要重复使用 (9)/(11)。剩余是原 \(\nu_b,\Omega\) 上的 signed 多面 difference。其平方自由系数
\(r^k[(1-r)^{n-k}-e^{-nr}]\) 若为正则 \(k>nr/2\)：
\(-\log(1-r)\ge r/(1-r/2)\)，由逐幂系数 \(1/\ell\ge2^{-(\ell-1)}\) 可核；端点 \(r=1\) 只有 \(k=n\) 正项，亦满足该下限。这是核 coefficient 资格，不是 actual 诊断门，更不是空间费。删除所有 repeated 的负核只得到 positive Bernstein 外放松，其全参数强包络仍有既有 \(\sqrt n\) 费用；不将新下限称 polylog 端点。

## 8. 终态后解析升级：所有低 support 一次支付

根在本稿守卫执行终态后提供的一般升级，独立解析核对通过；**没有修改已执行脚本/结果，也不声称下式包含于其 523 项**。根另负责该新合同的独立守卫。

对每个非空 mask \(A\)、\(k=|A|\)，精确有
\[
 D_{r,A}=d_k(r)G_A
 -e^{-nr}\!\!\sum_{\substack{\operatorname{supp}m=A\\|m|\ge k+1}}
       \frac{r^{|m|}}{\prod_i m_i!}G^m,\qquad
 d_k=r^k[(1-r)^{n-k}-e^{-nr}].
\tag{15}
\]
此式保留全部重复跳为负项，不取其绝对值。因而任意固定非负 \(L^1\) 来源 \(\rho\) 都有
\[
 \left\|\sup_r[D_r^{1\le|A|\le K}\rho]_+\right\|_1
 \le W_\rho\sum_{k=1}^K\binom nk\sup_r[d_k(r)]_+
 \le \frac{K(K+1)(K+2)}{3(n+1)}W_\rho.
\tag{16}
\]
最后的全参数 bound 自含如下。令
\(\eta=(n-k)\log(1-r)+nr\)。正系数要求 \(\eta>0\)，并且 \(\eta\le kr\)。故
\[
 [d_k]_+\le k r^{k+1}(1-r)^{n-k}.
\]
右侧最大点 \(r_k=(k+1)/(n+1)\)。乘 \(\binom nk\) 后不超过
\(k r_k\)，因为
\(\binom nk r_k^k(1-r_k)^{n-k}\le1\) 是 binomial PMF。
求和 \(\sum_{k\le K}k(k+1)=K(K+1)(K+2)/3\) 得 (16)。
\(r=0\) 是0；\(k<n,r=1\) 的 positive coefficient 为0；
\(k=n,r=1\) 按 polynomial 端点解释亦成立，无 log 除零。

因此 \(K=\lfloor n^{1/3}\rfloor\) 时右侧不超过 \(W_\rho\)；
费用没有 source atom count、mask 数或深度因子。这覆盖随 \(n\) 增长的所有低 support、整个连续参数以及其全部 repeated-negative 部分，明确不是固定 \(k\) 特殊族。适用于原同一 \(\rho=\nu_b\)，不是另选 mask 依赖输入。

直接 D 路线可用
\[
 D_r\nu_b\le P_K\nu_b+[D_r^{|A|>K}\nu_b]_+,
 \qquad \|P_K\nu_b\|_1\le C_KW_b,\quad
 C_K=\frac{K(K+1)(K+2)}{3(n+1)}.
\]
若剩余高 support 合同
\(\kappa|\{\Omega^c:\sup_r[D^{|A|>K}\nu_b]_+>\kappa\}|
 \le C_{\rm high}W_b\) 成立，则与 (1e) 同样回代给
\[
 A_{\rm ord}\le4+4A_{\rm fix}+4C_K+4C_{\rm high}.
\tag{17}
\]
高 support 晚参数部分仍未付；(17) 是条件预算，不是一般上界。此前 \(14/n\) 单面费用被 (16) 的统一低 support 支替换，不能另加一笔。

作为范围检查，取 \(b=K/(4n)\)，\(r\le b\) 时 \(k>K\) 的 squarefree 系数在该窗单调。因此其正包络质量不超过
\[
 W_b\,\mathbb P\{\operatorname{Bin}(n,b)>K\}
 \le\frac12(2/3)^K W_b.
\]
用 \(2^{\rm count}\) 的 MGF 与
\((1+b)^n\le e^{K/4}\le(4/3)^K\) 即得。
配 (16) 给 \(r\lesssim n^{-2/3}\) 的 constant 强费，但该早窗本身已经是 Gamma 稿 product-telescoping 导数预算删去 \(I+S\) 后的推论；不重复登记早窗为新主进展，也不为此重跑数值。

## 9. 新证书终态与范围

先登记 [注册](late_multiface_source_budget_registration_20261007.json)，再一次执行 [新脚本](late_multiface_source_budget_exact_guard_20261007.py)，[结果](late_multiface_source_budget_results_20261007.json) 已终态 PASS，523 项精确 Fraction/整数/形式符号断言，无失败。三轮有限代数维度 \(3,4,5\)，压缩维数 \(512,32768,2097152\)，各轮 \(89,153,281\) 项。脚本 SHA256 为 d14e077e5077a97cdaf939307bf3699b4e6eeac9c144f0609eefad691850efde。无 live handle。

注册在首次运行前加入根的 direct-D 14/n 系数范围检查；尚无任何早版本结果。三轮 overlap 用 \(G_i=(2I+\text{coordinate flip})/3\) 的非幂等有限正 Markov 模型，只核 (3)/(4) 的 source/cap 身份；\(H\) 保留共同 \(\exp(-rs)\) 形式符号，准确比较其 constant/exp 两个系数，没有替换为另一 semigroup。

(14) 压缩比价三轮认证下界分别为 \(2^3,2^{15},2^{63}\)；保存巨大整数的 bit-length/SHA，未输出其全部位数。它只证明 \(Q^j\) 正比较可能产生的预算损失，不是原输入、weak 或 actual FIRST 反例。既有核/单面/Gamma 数值均未重跑。

守卫核 overlap-source signed 身份、全单面恢复费用、direct-D 解析常数及 (14) 的压缩大维比较价。**§8 是运行终态后的新解析推论，不属于523項守卫。** 普遍命题由本文证明承担，有限模型不认证晚多面弱费用或 actual cube 门。
