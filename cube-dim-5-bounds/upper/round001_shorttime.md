# 第 001 轮：五维 body frame 的全方面短时预算

令
\[
 V=\{(1,\epsilon_2,\epsilon_3,\epsilon_4,\epsilon_5):\epsilon_j\in\{-1,1\}\},
 \qquad\psi(\xi)=\sum_{v\in V}|v\cdot\xi|.
\]
每条方向权重为一，全部一维 Cauchy 变量共享时间 \(t\)。本文件只建立后续截时积分的接口；它没有给出五维 Green 势闭式、方体下界或新的弱型上界。

所有十六个方向的第一坐标均为一，故真实热核的第一坐标边际精确为
\[
 \int_{\mathbb R^4}p_t(1,y)dy=\frac{16t}{\pi(1+256t^2)}\le\frac{16t}{\pi}. \tag{1}
\]
令 \(k=16\)。若随机表示命中 \(x_1=1\)，至少一个 Cauchy 参数的绝对值至少 \(1/k\)。`round001_shorttime.py` 用整数 Bareiss 消元逐项穷举证明：每条方向均可补四条其他 body 方向组成一个行列式绝对值为 \(48\) 的基。该命中参数的密度至多 \(k^2t/\pi\)，四个补坐标密度各至多 \((\pi t)^{-1}\)，余变量仅产生概率卷积。因此
\[
 p_t(1,y)\le C_Ft^{-3},\qquad
 C_F=\frac{k^2}{\pi^5}\sum_{v\in V}\frac1{48}
 =\frac{256}{3\pi^5}. \tag{2}
\]
这里是对空间密度的事件并集估计，不把零测截面当作正概率事件。

令 \(q=5/4\)。由 (1)--(2)，
\[
 \|p_t(1,\cdot)\|_{L^q(\mathbb R^4)}
 \le C_F^{1/5}(16/\pi)^{4/5}t^{1/5}.
\]
故对 \(T_\varepsilon(1,y)=\int_0^\varepsilon p_t(1,y)dt\)，先在有界区域使用 Minkowski，再单调扩张，
\[
 \|T_\varepsilon(1,\cdot)\|_{L^q(\mathbb R^4)}
 \le\frac56 C_F^{1/5}(16/\pi)^{4/5}\varepsilon^{6/5}. \tag{3}
\]
最后四坐标独立反射对称，正正交象限占 \(L^q\) 次积分的 \(1/16\)，于是
\[
 \|T_\varepsilon(1,\cdot)\|_{L^q([0,\infty)^4)}
 \le16^{-4/5}\frac56\left(\frac{256}{3\pi^5}\right)^{1/5}
 \left(\frac{16}{\pi}\right)^{4/5}\varepsilon^{6/5}. \tag{4}
\]

这处理完整面及其无穷尾，包含一切面内奇性。因而可对平滑的
\(H_\varepsilon=\int_\varepsilon^\infty p_tdt\) 做连续计算认证，再以 (4) 支付短时余项。

几何枚举同时显示：16 条 rank-1 方向线、120 个由方向对张成的 rank-2 平面，且任意三条不同方向 rank 为 3。尺度律在 rank \(r\) 子空间处为 \(t^{2r-5}\)：方向线为 \(t^{-3}\)，方向对平面为 \(t^{-1}\)。所以五维 Green 势在 rank-2 平面可出现真实对数发散，不能把四维“只在线上发散”的连续性论证照搬。式 (3)--(4) 不要求此连续性。

## 复现

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper/round001_shorttime.py
```

输入是脚本内十六个整数方向，输出为锁定 JSON 及 SHA-256。运行不使用浮点算术。
