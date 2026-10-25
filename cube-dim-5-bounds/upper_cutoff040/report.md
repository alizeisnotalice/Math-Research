# upper_cutoff040：穿过零投影的一个截时势闭盒

本轮只生产一个四变量闭盒，固定 \(\varepsilon=1/512\)、半宽 \(h=1/2048\)、Taylor 总阶数 3、384-bit Arb。得到的对象是

\[
H_\varepsilon(x)=\int_\varepsilon^\infty p_t(x)dt,
\]

不是完整 \(G\)，也不是正部 \((G-c)_+\)。设

\[
x_0=\left(1,\frac{78547}{104729},\frac{78545}{104729},
              \frac{52273}{104729},\frac{13113}{104729}\right),\qquad
B=\{x_0+(0,u):|u_j|\le1/2048\}.
\]

这个闭盒严格位于正单位面内，准确体积为 \(2^{-40}\)，并穿过真顶点投影超平面 \(x_2=x_3\)。全部 3504 锥、完整截时区间 \([1/512,\infty)\) 的严格结果为

\[
\boxed{0.000218951499<H_{1/512}(x)<0.000235519403\quad(x\in B),}
\]
\[
\boxed{\int_BH_{1/512}^{5/4}
<\frac{53}{2\cdot10^{18}}=2.65\times10^{-17}.}
\]

更紧的积分上包表达式为 \(2.6491550704986359416\ldots\times10^{-17}\)。全部可用端点是 JSON 中的准确 dyadic；显示小数不是证书端点。没有全域 mesh、完整 \(I\)、全域 \(g\) 或新 \(C_5\) 结论。

## 1. 预登记、中心拒绝与只读输入

`registration.json` 在实际计算前锁定固定参数与最多两个备选。原建议中心的第四坐标分子为 52271；准确中心资格检查发现其 280 个投影只有 278 个互异，不能用本轮普通中心 DD 递推。因此程序在构造任何 stencil 之前记录拒绝，采用预登记的第一个备选：只将该分子加 2 成为 52273。该中心全部 280 个投影互异，最小绝对整数投影为 18。第二个备选（减 2）没有执行，没有调整盒宽、阶数、位数或截时参数。

原中心并非使 \(H_\varepsilon\) 无定义；拒绝只是本实现使用互异中心节点的准确 DD 递推的资格条件。截时势本身在重复/零投影处也光滑。

只读使用 017 已锁的表：280 个整数顶点 \(z=144r\)、3504 个锥、各锥正行列式分子 \(d_C=|\det(z_1,\ldots,z_5)|\)。表 SHA256 为
`8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9`。
035 的报告及代码作为只读数学/实现参考；本轮没有执行 017、018、035 生产器、复制输入快照、重算旧点或旧几何。所有这些源的字节 hash 都在登记与结果中。

原定中心与备选的准确资格在 `qualification.json`。此前向主代理的即时消息曾过早称“原中心直接合格”，读回该 JSON 后已明确更正；数值始终对应上述第一个预登记备选，代码与结果没有变动。

## 2. 截断 DD 公式与全部尺度

固定五维 16 条未单位化 body 方向、权 1、\(\alpha=1\)、Lebesgue 测度、Fourier 逆因子 \((2\pi)^{-5}\)。017 的锥积分给

\[
p_t(x)=(2\pi)^{-5}\operatorname{Re}\sum_Cd_C
          \prod_{j=1}^5(144t-i\,x\cdot z_j)^{-1}.
\]

令 \(D=104729\)、\(k_v=D x_0\cdot z_v\in\mathbb Z\)、
\(l_v=D(z_{v,2},\ldots,z_{v,5})\)，并作变量代换

\[
s=144Dt,\qquad a=144D\varepsilon=942561/32>0,
\quad F=\frac{D^4}{144(2\pi)^5}.
\]

于是整个盒上的截时势为

\[
H_\varepsilon(x_0+(0,u))=
F\operatorname{Re}\sum_Cd_C\int_a^\infty
\prod_{j=1}^5(s-i k_j-i l_j\cdot u)^{-1}ds. \tag{1}
\]

五个分母提供 \(D^5\)，时间代换提供 \((144D)^{-1}\)，所以是 \(D^4/144\)，没有额外 \(144^4\) 或 \(5!\)。下节的正绝对积分界同时证明逐锥在这个截时区间的绝对可积性。

对互异实节点 \(b_j\)，部分分式系数为
\(1/\prod_{m\ne j}(b_j-b_m)\)：五因子产生的 \(i^4=1\)。积分至无穷时共同的 \(\log s\) 因这些系数之和为零而消失。因此

\[
\operatorname{Re}\int_a^\infty\prod_j(s-i b_j)^{-1}ds
=-\operatorname{DD}_4 h_a(b_1,\ldots,b_5),\qquad
h_a(b)=\tfrac12\log(a^2+b^2). \tag{2}
\]

在原未缩放变量中，同式使用
\(h_\varepsilon(b)=\tfrac12\log(\varepsilon^2+b^2)\)。由
\(h_\varepsilon(b/(144D))=h_a(b)-\log(144D)\)，共同常数被四阶 DD 消去，DD 节点缩放给 \((144D)^4\)，与 (1) 的因子一致。

因为 \(a>0\)，\(h_a\) 在全部实轴解析。重复节点对应 Hermite 延拓，零节点同样有限；不应把 \(\log|b|\) 的全 Green 原函数代入这里。式 (2) 的中心普通 DD 只用于准确提取 jet；盒内有效性由 (1) 的 resolvent 证明，不要求节点在盒中继续互异。

## 3. 四变量区间 jet：有理系数和 log 常数

记 \(d=a^2+k^2\)。中心所需的三个导数准确为

\[
h_a'(k)=k/d,\quad h_a''(k)=(a^2-k^2)/d^2,\quad
h_a'''(k)=2k(k^2-3a^2)/d^3.
\]

对 \(1\le|\alpha|\le3\)，

\[
[u^\alpha]h_a(k+l\cdot u)
=\frac{h_a^{(|\alpha|)}(k)}{\prod_j\alpha_j!}
 \prod_jl_j^{\alpha_j}.
\]

这些是准确有理数。常数项写成
\(\tfrac12\log(A^2+(32k)^2)-\log32\)，其中 \(A=942561\)；最后的共同常数在 DD 中准确消失。每个 DD 除数是 \(d_0+d_1\cdot u\)，其中中心 \(d_0\ne0\)。若 \(Q=P/(d_0+d_1\cdot u)\)，按总次数递推

\[
Q_\alpha=\frac{P_\alpha-\sum_{j:\alpha_j>0}d_{1,j}Q_{\alpha-e_j}}{d_0}.
\]

`cutoff_box.py` 使用 python-flint `fmpq`，对全部 3504 锥先准确合并 35 个总次数不超过 3 的系数 stencil，再以 384-bit Arb 计算 log 与 \(F\)。最终形式是

\[
A_\alpha=F\left(c_{\alpha,0}+\sum_{n>0}c_{\alpha,n}\log n\right),
\qquad c_{\alpha,n}\in\mathbb Q.
\]

“准确有理”指 log stencil 的权重，不是把含 log 和 \(\pi^{-5}\) 的 Taylor 系数称为有理数。本轮 4934 个非零 stencil 项，最大整数 29399 bits。全部 35 个最终系数球已保存；完整准确 stencil 留在外置 checkpoint。

中心实际得到

\[
H_{1/512}(x_0)=0.00022723545107106519949366601909964609\ldots,
\]

其球仅认证这个中心截时值及算术，不是完整 \(G(x_0)\)。

## 4. 允许零投影及重复节点的完整截时余项

对 \(s\ge a\) 和整个闭盒，

\[
\left|\frac{i l_v\cdot u}{s-i k_v}\right|
\le\frac{hD\|z_{v,\rm tail}\|_1}{\sqrt{a^2+k_v^2}}.
\]

代码用 Arb 平方根的准确有向下端 \(b_v\in\mathbb Q\) 满足

\[
a\le b_v\le\sqrt{a^2+k_v^2},\qquad
\bar q_v:=hD\|z_{v,\rm tail}\|_1/b_v.
\]

所有关系在有限 280 个顶点上实际验证。结果准确有理 \(\bar q_{\max}<1/32\)，显示约
\(0.03124999416493176648\)。这里既没有除以投影 gap，也没有要求 \(|k_v+l_v\cdot u|>0\)。

对一个锥，取 \(q_C=\max_{v\in C}\bar q_v\)。五个 resolvent 几何级数的总次数大于 3 的尾满足

\[
|\text{尾}|
\le\prod_j|s-i k_j|^{-1}T_3(q_C),
\quad
T_3(q)=(1-q)^{-5}-(1+5q+15q^2+35q^3). \tag{3}
\]

总次数 \(n\) 的项数为 \(\binom{n+4}4\)，所以 (3) 是四变量空间 Taylor 的完整高阶尾，不是再截掉某段时间。\(T_3(q_C)\) 用准确 `fmpq` 计算。

将该锥的五个 \(b_j\) 升序排列。由于 \(s\ge a\)，

\[
|s-i k_j|^{-1}\le\min(b_j^{-1},s^{-1}).
\]

在 \(a,b_1,\ldots,b_5\) 分段积分，得到

\[
\begin{aligned}
\int_a^\infty\prod_j|s-i k_j|^{-1}ds
\le J_C={}&\frac{b_1-a}{\prod_jb_j}
 +\frac{\log(b_2/b_1)}{\prod_{j=2}^5b_j}\\
&+\sum_{m=2}^4
\frac{b_m^{1-m}-b_{m+1}^{1-m}}
{(m-1)\prod_{j=m+1}^5b_j}
 +\frac1{4b_5^4}. \tag{4}
\end{aligned}
\]

相等断点时该段为零。式 (4) 的有理部分准确计算，log 用 384-bit Arb；取整个式子的 dyadic 上端 \(\bar J_C\)。因此这个正项上包不依赖数值时间积分。最终采用

\[
R_{\rm cert}=F\,\operatorname{up}_{\rm dyadic}
\left(\sum_Cd_C\bar J_CT_3(q_C)\right).
\]

于是对所有 \(x=x_0+(0,u)\in B\)，

\[
|H_\varepsilon(x)-P_3(u)|\le R_{\rm cert}.
\]

\((1-q_C)^{-5}\bar J_C\) 同时提供一致可积支配，允许时间积分与空间展开交换。中心形式 DD 系数由 Taylor 唯一性等于 resolvent 系数，因此即使盒内 DD 的个别差分母为零，也无需在该点执行普通 DD。\(s\ge a>0\) 保证真实 resolvent 无零分母，所有零投影和重复节点都被覆盖。

## 5. 准确跨零见证

017 真顶点按表中从 0 起编号，125 为 \(z=(0,-9,9,0,0)\)，154 为其相反数。中心整数投影为 \(k=-18\)，整个盒的该整数投影区间为

\[
[-960993/1024,\;924129/1024],
\]

严格包含 0。取盒内点

\[
x_*=(1,78546/104729,78546/104729,52273/104729,13113/104729).
\]

相对中心的两个非零扰动为 \(\pm1/104729\)，严格小于 \(h\)。两个相反真顶点在 \(x_*\) 的投影都准确为零，因此也给了重复节点见证。实际共有 4 个真顶点投影在盒内严格跨零。这是锥表达的零投影平面，不据此声称 \(G\) 或 \(H_\varepsilon\) 在那里奇异。

## 6. 闭盒范围和盒积分上包

取

\[
V_3=\sum_{1\le|\alpha|\le3}|A_\alpha|h^{|\alpha|}.
\]

用准确有向端点得到

\[
L=\max(0,\operatorname{down}(A_0-V_3-R_{\rm cert})),
\quad U=\operatorname{up}(A_0+V_3+R_{\rm cert}).
\]

真实热核为非负密度，故 \(H_\varepsilon\ge0\) 可用于下端裁切；本例本已严格为正。实际：

|量|上包表达式或范围端点近似|
|---|---:|
|\(V_3\)|\(3.17621539477132\times10^{-7}\)|
|\(R_{\rm cert}\)|\(7.96633010936616\times10^{-6}\)|
|uniform lower|0.000218951499422222|
|uniform upper|0.000235519402719908|

下/上界表达式球各自很窄，只表示一个统一端点表达式的算术包围；不表示整盒函数值落在某个窄球里。余项明显主导本例范围宽度，不把它当成真实函数变动。

闭盒的多项式均值准确利用奇次矩为零、偶次矩

\[
\frac1{|B|}\int_Bu^\alpha du
=\prod_j\frac{h^{\alpha_j}}{\alpha_j+1}
\quad(\text{全部 }\alpha_j\text{ 为偶数}),
\]

得到均值球
\(\bar P=0.00022723545440329464470873928396488070\ldots\)。因此 \(H_\varepsilon\) 的均值上包取
\(m=\min(U,\operatorname{up}(\bar P+R_{\rm cert}))\)。对 \(q=5/4\) 和所有 \(z\in[L,U]\)，凸弦给

\[
z^q\le L^q+\frac{U^q-L^q}{U-L}(z-L).
\]

所以

\[
\int_BH_\varepsilon^q
\le |B|\left[L^q+\frac{U^q-L^q}{U-L}(m-L)\right].
\]

这一步用的是凸弦上界，不是将 Jensen 的均值下界反向使用。所有幂、斜率与求和向外包围，最后取 dyadic 上端。结果为

\[
[2.64915507049863594161549142340530357748\ldots\times10^{-17}
\ \pm9.00\times10^{-131}].
\]

这个球包围可用上界表达式；真实积分只保证不超过其严格上端。仅用常数值域的 \(|B|U^q\) 约为 \(2.65358826264700\times10^{-17}\)，本轮均值弦略紧。没有推断全域积分或以本例成本外推 mesh 总成本。

## 7. 实际运行、锁定与复现

实际一次生产：PID 44037，全部 3504 锥，约 9.6196 秒，峰值 RSS 51494912 bytes，python-flint 0.9.0，384 bits。每 128 个完整锥保存外置 checkpoint；未触发 120 秒停止，也没有为了测试续算重跑同批。

外置准确缓存为
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cutoff040/checkpoint.json`，2474868 bytes，SHA256
`c02912118c9f40bba1d3a2b91ad5eec3505f450a5da91dc2075535d31c69954e`。
只在本目录归档 35 个系数球、单盒端点、积分上包、资格与源 hash，不归档该可重建缓存。

`readback_checks.json` 用 Fraction 读回可用 dyadic，准确验证 \(q_{\max}<1/32\)、本报告三个短有理端点、体积 \(2^{-40}\)、已登记备选、跨零见证及源码/登记/外部源 hash。读回不是第二次生产。复现命令是

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_cutoff040/cutoff_box.py \
  --work-dir /path/outside-git/new-cutoff040-work
```

若实际到达 120 秒，程序完成当前锥后退出 75；以相同源码、登记、输入和 work 路径加 `--resume`，只继续下一未完成锥。源码/输入指纹改变则拒绝缓存。不要将已用缓存目录在无 `--resume` 时重用。源路径是登记中的只读项目路径，移机需映射同 hash 字节；旧源未运行。

本轮全部新文件只在 `upper_cutoff040/`；034、035、017/018、主 clone 与他人文件未修改。无 MATLAB、编译二进制、空间 mesh、额外中心搜索或提交。被省略的小时间 \((0,1/512)\) **不属于当前对象**；未来若要认证完整 \(I_{\rm face}\)，必须另接 034 已明确的完整短时范数和奇异管/覆盖接口，不能把本盒 \(H_\varepsilon\) 上包称为全 Green 上包。
