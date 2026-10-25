# 第005轮：删一参数的 whole-face 短时界

本节替换第002轮在 face slab 上按32个大跳事件作并集的估计。新论证先在**全空间密度**层面求和，再限制到正 face；因而不会在奇异子空间或重叠处重复计数。

设 \(V=\{v_i=(1,v_i'):1\le i\le32\}\)，令 \(Z_i\) 独立且密度 \(f_t(z)=t/[\pi(t^2+z^2)]\)，并令 \(X=\sum_i v_iZ_i\) 的六维密度为 \(p_t\)。对实 \(z_1,\ldots,z_N\) 和 \(s=\sum_i z_i\)，有
\[
 \prod_{i=1}^N f_t(z_i)\le {Nt\over\pi(s^2+N^2t^2)}
 \sum_{i=1}^N\prod_{j\ne i}f_t(z_j). \tag{1}
\]
除以左端后，所需不等式等价于
\[
 {N\sum_i(z_i^2+t^2)\over s^2+N^2t^2}\ge1,
\quad N\sum_i z_i^2-s^2=N\sum_i(z_i-s/N)^2\ge0. \tag{2}
\]
所以这是点态恒等式余项，非 union bound。

固定 \(i\)，由 \(z_i=s-\sum_{j\ne i}z_j\) 消元，有
\[
 y-sv_i'=\sum_{j\ne i}(v_j'-v_i')z_j. \tag{3}
\]
令 \(q_{i,t}\) 是这31个变量在五维商方向的 push-forward density。delta/coarea 表示中第一行 \(\sum z_i=s\) 对被消变量导数是1，故没有额外 Jacobian；对紧支撑测试函数先用 Fubini、再用非负 Tonelli，得到
\[
 p_t(s,y)\le {32t\over\pi(s^2+1024t^2)}
 \sum_{i=1}^{32}q_{i,t}(y-sv_i'). \tag{4}
\]
对任何 \(i\)，符号反射重标将 quotient 方向变成 \(\{0,-2\}^5\setminus\{0\}\)。第002轮的 Fourier estimate \(\psi_Q(\zeta)\ge12\|\zeta\|_1\) 因而统一给出
\[
 \|q_{i,t}\|_\infty\le(12^5\pi^5t^5)^{-1}. \tag{5}
\]
在 \(s=1\) 代入 (4)：
\[
 \|p_t(1,\cdot)\|_\infty
 \le {32^2\over12^5\pi^6}t^{-4}
 ={1\over243\pi^6}t^{-4}. \tag{6}
\]
两个32分别来自 (1) 的系数及其求和，故没有第三个32。这严格优于旧的 \(32/(243\pi^6)\)。

边际的精确等式仍是
\[
 \int_{\mathbb R^5}p_t(1,y)\,dy={32t\over\pi(1+1024t^2)}\le {32t\over\pi}. \tag{7}
\]
它可先积分五个自由物理坐标，或直接由非负联合密度 Tonelli 得到，所以不依赖 face slice 的点态定义。free-coordinate 反射不变性把全空间 \(L^{6/5}\) norm 转到正正交体时带来 \(32^{-5/6}\)，恰与 (7) 的 \(32^{5/6}\) 相消。令 \(T_\epsilon(1,y)=\int_0^\epsilon p_t(1,y)\,dt\)，Minkowski 及 \(\|h\|_{6/5}\le\|h\|_\infty^{1/6}\|h\|_1^{5/6}\) 给
\[
 \|T_\epsilon(1,\cdot)\|_{L^{6/5}([0,1]^5)}
 \le K_{\rm new}\epsilon^{7/6},\qquad
 K_{\rm new}={6\over7}(1/243)^{1/6}\pi^{-11/6}. \tag{8}
\]
逐时被积范数为 \(O(t^{1/6})\)，所以 \(t=0\) 可积。这是在 \(t>0\) 全空间平滑密度上完成的证明，未逐点排除或忽略 \(S_2\) 及其交叠。

本改进是 universal whole-face tail budget，可替换第002轮 tail 常数。由于 \(K\) 含 \(C_f^{1/6}\)，相对旧 \(K\) 的因子是 \(32^{-1/6}=2^{-5/6}\)，不是 \(1/2\)。它尚未产生 \(H_\epsilon\) 的全 face 积分、cube 的 \(g\) 下界或新的 \(C_6\) 上界。short_time_product_audit.py 只复核准确常数与代数（零 kernel evaluation）。
