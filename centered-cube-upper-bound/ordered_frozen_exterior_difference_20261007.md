# 原 ordered–frozen 域外差：共同 resolvent 来源、负部支付与真实坐标面

2026-10-07。接 [上轮谱差定理](ordered_frozen_spectral_difference_20261007.md)，本文在原 obstacle 来源合同内支付整个负来源分支，并把剩余正来源准确缩到非局部核。原 \(G_1\) 的严格非局部负密度证书排除正核捷径，不是原弱端点反例。尚未得到 \(A_{\rm ord}=\mathrm{polylog}(n)\)，没有加入 cube 主账。不重做 Abel/obstacle，不使用已被 spike 排除的 potential-energy cap。

## 1. 同一来源与 Markov resolvent

固定原物理尺度、\(c>0\)，沿用
\[
A_i=I-G_{c,i},\quad S=\sum_iA_i,\quad
K_r=\prod_i(I-rA_i),\quad H_r=e^{-rS},\quad 0\le r\le1.
\tag{1}
\]
原对称保质量轴向卷积的乘子是 \(g_c(\xi_i^2)=[1+cB(\xi_i^2)]^{-1}\)，其中 \(B(v)=v/\log(1+v)-1\)。使用已构造的原 fixed obstacle：\(u\in L^1\cap L^2\)、\(u\ge0\)、\(\Omega=\{u>0\}\)，以及
\[
\sigma=Su=\nu_b-\mu_b,\quad \nu_b=\nu\mathbf1_\Omega,\quad
0\le\mu_b\le\kappa,\quad
\int\nu_b=\int\mu_b=W_b,\quad |\Omega|\le W_b/\kappa.
\tag{2}
\]
完整输入仍有 good 部分 \(\nu_g\)，且 \(\mu_{\rm tot}=\nu_g+\mu_b\le\kappa\)、\(\nu=\mu_{\rm tot}+\sigma\)。不能把含 good 部分的总 \(\nu\) 当 killed 方程内部源。不添加共同帽，不让来源随 receiver 或 \(r\) 更换。

令
\[
R=(I+S)^{-1}=\int_0^\infty e^{-t}H_t\,dt,\qquad
h=R\sigma=(I-R)u.
\tag{3}
\]
\(R\) 正、对称、保质量，故
\[
\|h\|_1\le2W_b,\quad \int h=0,\quad h\ge-\kappa.
\tag{4}
\]
最后一式来自全空间 \(\sigma\ge-\kappa\)、\(R1=1\)，不是只用域外符号。在 \(\Omega^c\)，\(h=-Ru\le0\)。因此
\[
\boxed{\operatorname{supp}h_+\subset\Omega,\quad
0\le h_-\le\kappa,\quad
\int h_+=\int h_-\le W_b,\quad
\|h_-\|_2^2\le\kappa W_b.}
\tag{5}
\]
还保留 \(h_+\le R\nu_b\)、\(h_-\le R\mu_b\)。正部是同一 \(\nu_b,\mu_b,u\) 生成的来源，不能任意重选。

相比旧 \(h_{\rm old}=(I-H_1)u=\int_0^1H_t\sigma\,dt\)，新 \(h\) 使后续只需乘有限阶 \(I+S\)。旧 \(h_{\rm old}\) 其实也满足 \(h_{\rm old}\ge-\kappa\)，且域外非正；上轮“没有 cap-L2”指整个 signed source/高正部没有该界，不妨碍负部单独封顶。本稿不改已冻结旧稿及收据。

## 2. 负来源完整支付，正来源保持原合同

定义 \(N_r=(K_r-H_r)(I+S)\)，准确有 \((K_r-H_r)\sigma=N_rh\)。上轮谱界除以谱值 \(s\)，给
\[
\int_0^1\!\left(|k_r-e^{-rs}|^2+
|r\partial_r(k_r-e^{-rs})|^2\right)\frac{dr}{r}
\le40\min(1,s^{-2}).
\]
乘 \((1+s)^2\)，使用 \((1+s)^2\min(1,s^{-2})\le4\)，再用同一 log-\(r\) 基本微积分，得到
\[
\boxed{\left\|\sup_{0\le r\le1}|N_rg|\right\|_2^2
\le160\|g\|_2^2.}
\tag{6}
\]
零谱的实际差为0。all-\(r\) 版本沿用 bounded-\(S\) exp 幂级数的正 majorant，给几乎处处紧参数上的绝对一致收敛，不只依赖 strong \(L^2\) 连续。

由(5)、Chebyshev，
\[
\kappa|\{\sup_r|N_rh_-|>\kappa\}|\le160W_b.
\tag{7}
\]
原域外事件
\[
E_{\rm diff}=\{x\in\Omega^c:
\sup_r[N_r(h_+-h_-)(x)]_+>2\kappa\}
\]
严格阈值对应的事件包含于
\[
\{x\in\Omega^c:\sup_r[N_rh_+(x)]_+>\kappa\}
\cup\{\sup_r|N_rh_-|>\kappa\}.
\tag{8}
\]
不宣称集合包含严格，两个比较量等于阈值时仍合法；无需 \(N_r\) 正。真正剩余的充分合同是对(2)–(5)同一 obstacle 正部证明
\[
\kappa|\{x\in\Omega^c:\sup_r[N_rh_+(x)]_+>\kappa\}|
\le C_nW_b.
\tag{9}
\]
若它成立，上轮域外差系数可取
\[
\boxed{B_n\le2C_n+320.}
\tag{10}
\]
这里 \(320\) 来自总事件的 \(2\kappa\) 乘(7)体积，不是 \(160\)；没有扩大正部合同为任意非负来源。旧 \(\kappa=\lambda/6\) 的条件回代为
\[
A_{\rm ord}\le2A_{\rm fix}(n)+6+3B_n
\le2A_{\rm fix}(n)+966+6C_n.
\tag{11}
\]
这是条件结论。已付的新分支是(7)，不是(11)。

## 3. 域外 holding 项消失及真实多跳系数

记 \(Q=\sum_iG_i\)，则 \(I+S=(n+1)I-Q\)、\(H_r=e^{-nr}e^{rQ}\)。原 \(G_c\) 有连续一维跳密度，故 \(N_r\) 的纯原子系数为
\[
\alpha_n(r)=(n+1)[(1-r)^n-e^{-nr}]\le0,\qquad
|\alpha_n(r)|\le2.
\tag{12}
\]
当 \(r\le1/2\)，\(-n\log(1-r)-nr\le nr^2\) 给
\[
|\alpha_n(r)|\le(n+1)nr^2e^{-nr}
\le2(nr)^2e^{-nr}<2.
\]
当 \(r\ge1/2\)，界为 \((n+1)e^{-n/2}\le2\)，由 \(e^{n/2}\ge1+n/2\) 即得。端点包含。因此 \(N_r=\alpha_n(r)I+N_r^{\rm nl}\)，且
\[
N_rh_+=N_r^{\rm nl}h_+\quad\text{于 }\Omega^c.
\tag{13}
\]
这是域外真正删项，不是由 \(u=0\) 推断所有跳项消失。整个负部已由(7)支付，无需其原子项再收来源质量。

对多重指标 \(m=(m_i)\)、\(j=|m|\)，令 \(G^m=\prod_iG_i^{m_i}\)，
\[
b_r(m)=
\begin{cases}r^j(1-r)^{n-j},&m_i\in\{0,1\}\ \forall i,\\0,&\text{否则},\end{cases}
\qquad p_r(m)=e^{-nr}\frac{r^j}{\prod_i m_i!}.
\]
写 \(d_r=b_r-p_r\)，有
\[
N_r=\sum_m a_r(m)G^m,\qquad
a_r(m)=(n+1)d_r(m)-\sum_{i:m_i\ge1}d_r(m-e_i).
\tag{14}
\]
固定 \(r\) 时级数在有限 signed-measure/算子 \(L^1\) 范数中绝对收敛，并有 \(\sum_m|a_r(m)|\le4n+2\)、\(\sum_m a_r(m)=0\)。前者只是固定参数的线性粗界，后者是质量 cancellation，均不能支付连续 maximal 的 polylog 费用。

独立核对的闭式为
\[
[(I+S)H_r]_m=e^{-nr}\frac{r^j}{\prod_i m_i!}
\left(n+1-\frac jr\right),\quad r>0.
\tag{15}
\]
对 \((I+S)K_r\)，squarefree support \(k\ge1\) 的系数是
\[
(n+1)r^k(1-r)^{n-k}-kr^{k-1}(1-r)^{n-k+1}
=r^{k-1}(1-r)^{n-k}[(n+k+1)r-k].
\tag{16}
\]
仅一个坐标为2、support 为 \(k\)、总 degree 为 \(k+1\) 的系数是 \(-r^k(1-r)^{n-k}\)；其余非 squarefree 形式均为0。\(r=0\) 用 polynomial recurrence，不用(15)除法。尤其(16)括号为 \(n+k+1\)，不是 \(n+1\)。

按 support 合并时，每个非空 support 是一个真实坐标面，\(G_i^{m_i}\) 的原 convolution-power 密度必须保留。重复跳不能当成同一个 \(G_i\)，不同维坐标面的奇异测度不能混为一个全维密度。坐标标签来自生成元内部跳跃，不假设源坐标独立，也不是辅助 killed operator 或下界模型的源标签。

## 4. 原 \(c=1\) 非局部负密度：正核捷径失败

此节只排除 \(N_r^{\rm nl}\ge0\)，不排除(9)。原 \(c=1\)、单位物理尺度时
\[
g_1(\xi^2)=\frac{\log(1+\xi^2)}{\xi^2}
=\int_0^1\frac{dt}{1+t\xi^2}.
\]
故原密度为 Laplace 混合
\[
w(x)=\int_0^1\frac{e^{-|x|/\sqrt t}}{2\sqrt t}\,dt,\quad
w(0)=1,\quad \|w^{*m}\|_\infty\le1\ (m\ge1).
\tag{17}
\]
直接积分两 Laplace 密度 overlap，得
\[
w*w(0)=\int_0^1\!\int_0^1\frac{dt\,ds}{2(\sqrt t+\sqrt s)}
=\frac43(1-\log2).
\tag{18}
\]
两个不同坐标 \(i,j\) 的真实坐标面上，\(r^2\) Taylor 首项为
\[
-r^2\left[2G_iG_j-\frac12G_i^2G_j-\frac12G_iG_j^2\right].
\tag{19}
\]
这是从完整 \(K_r-H_r=-r^2\sum_iA_i^2/2+O(r^3)\)，乘 \(I+S=(n+1)I-Q\) 后提取两面 support 所得，不是逐跳猜符号。该面原点的首项密度为
\[
-r^2\left[2-\frac43(1-\log2)\right]<-r^2.
\tag{20}
\]
连续性给邻近负密度开集，剔除坐标轴后仍有正面测度；不把原点取值本身当正体积事件。

完整 \(N_r\) 的 \(r^3\) 及以上系数绝对和被
\[
\begin{split}
\mathcal R_n(r)
&=(2n+1)\sum_{\ell\ge3}
\left[\binom n\ell 2^\ell+\frac{(2n)^\ell}{\ell!}\right]r^\ell\\
&\le\frac{2n+1}{3}(2nr)^3e^{2nr}
\end{split}
\tag{21}
\]
上界，二项式项在 \(\ell>n\) 取0。每个两面 convolution-power 密度峰值至多1，故(21)也上界该面剩余密度的绝对值。这里用明确系数和及原密度峰值，不将一般 TV 界误当 pointwise density 界。

选 \(r_n=1/(100n^4)\)、\(n\ge2\)。由 \(e^{2nr_n}\le(1-2nr_n)^{-1}\)，
\[
\frac{\mathcal R_n(r_n)}{r_n^2}
\le\frac{(2n+1)8n^3r_n}{3(1-2nr_n)}<\frac12.
\tag{22}
\]
因此完整两面密度仍严格负。其他 support 的测度不能消除此面负开集：一面在坐标轴上，更高维连续面测度在这两面上质量0。原 \(N_r^{\rm nl}\) 确为 signed kernel，不能假设正。合法 \(c=1\) 子族已足以否定全 \(c,r\) 非局部正性；来源分解及(6)–(13)仍 uniform \(c>0\)。

## 5. 可调 resolvent 与精确未付项

可选 \(t>0\)：
\[
R_t=(I+tS)^{-1},\quad h_t=R_t\sigma=(u-R_tu)/t,\quad
N_{r,t}=(K_r-H_r)(I+tS).
\tag{23}
\]
同样有 \(\|h_t\|_1\le2W_b\)、\(\int h_t=0\)、\(h_t\ge-\kappa\)，域外 \(h_t=-R_tu/t\le0\)，两 signed masses 至多 \(W_b\)。谱界是
\[
\|\sup_r|N_{r,t}g|\|_2^2\le40(1+t)^2\|g\|_2^2,
\tag{24}
\]
因 \((1+ts)^2\min(1,s^{-2})\le(1+t)^2\)。负部在 \(\kappa\) 阈值下支付
\[
\kappa|\{\sup_r|N_{r,t}(h_t)_-|>\kappa\}|
\le40(1+t)^2W_b.
\tag{25}
\]
回代 \(2\kappa\) 总差事件体积须乘2。取 \(t=O(\log n)\) 仅使此已付分支增加 polylog 费用，允许更长 smoothing；不能由此称 \(R_t\) 是 allvisited。有限 \(t\) 的 holding 与 proper faces 仍存在。本轮未为此扩展实验。

真正未付的是
\[
\kappa|\{x\in\Omega^c:\sup_r[N_r^{\rm nl}h_+(x)]_+>\kappa\}|
\stackrel{?}{\le}\mathrm{polylog}(n)W_b.
\tag{26}
\]
它保留 \(h_+\subset\Omega\)、\(h_+\le R\nu_b\)、原 cap 来源与 joint symbol。局部项删掉、负部 L2 支付和总质量0均不足以证明它。按各面/多跳取绝对值只留下线性维数的粗 coefficient bound；现有连续 square function 只付封顶负部，高正部没有 \(\kappa W_b\) 的 \(L^2\) cap。

旧 fixed free–killed Abel exit measure 是 \(S_\Omega+\kappa/u\) 的单谱 resolvent 展开，\(K_r\) 则是 joint spectrum，不是 \(S\) 的单变量函数，压缩坐标也不再交换。新 exit/intertwining 预算必须控制(26)，不能先将高正源未来 maximal 交回未证 \(A_{\rm ord}\)。真实有序 killed 出口有质量 \(W_b\)，但其 receiver 选时仍返回此循环。本文不平均条件 weak 常数、不逐窗重收 \(W_b\)，也未将原时变 \(G_{c(1-r),i}\) 替为固定 \(G_c\) 后宣称 actual geom 已付。

## 6. 新三轮精确证书

先写 [注册](ordered_frozen_exterior_difference_registration_20261007.json)，再执行新 [guard](ordered_frozen_exterior_difference_guard_20261007.py)，终态 [结果](ordered_frozen_exterior_difference_results_20261007.json) 保存全部断言与区间。三轮 \(n=8,32,128\)，exp/log positive-series 项数32/128/512；参数 \(r=0,1/(100n^4),1/(2n),1/n,1/2,3/4,1\)。两面证书只用原 \(c=1\) 的(17)–(22)，不是 reset proxy，也不重跑旧原核数值。

6030条 exact Fraction/有理区间断言全部通过，各轮330/1170/4530。覆盖160/320封顶常数、原子符号/界、全部 support 大小的 squarefree/repeated recurrence、完整两面负密度区间。三轮密度除 \(r_n^2\) 的外包区间约为
\[
[-1.647532,-1.534194],\quad
[-1.645030,-1.536696],\quad
[-1.644405,-1.537321].
\]
十进制仅供阅读，认证字段保存精确有理端点。Exp 用正 Taylor 和作下界、首遗漏项除 \(1-x/(M+2)\) 作尾上界；\(\log2=2\operatorname{artanh}(1/3)\) 用正奇幂和及几何尾；负指数用精确 reciprocal 反向。没有浮点 quadrature/Gamma 或事后拟合。

首执行 exit1 仅因 Python 默认4300位整数转录限制：写很长 Fraction 端点时触发 ValueError，没有结果文件或冻结 hash。新脚本加入 sys.set_int_max_str_digits(0) 后，同一注册策略完整执行 exit0；明确保留此 serialization 修订记录。不改旧脚本/数据，不运行新 obstacle solver 或正部 maximal 实验。

注册 SHA256：dcf3b3a2e2e9cc30240651ebb5679088f54505e9182ca87cf7cc43640d55e06b。
终态脚本 SHA256：c2836fedb3795ea2738b0ec18be664ca6609c55080b15726855460f550ac323d。
这是负部支付与原 signed face 结构的证书；不认证(26)、完整 ordered 弱端点或 actual cube 余项。
