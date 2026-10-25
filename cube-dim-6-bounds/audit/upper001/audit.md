# 独立审计：六维 body frame、二维奇异层与 whole-face 短时预算

审计日期：2026-10-25。审计对象是理论组的
`upper/theory/001_n6_body_frame_and_singular_strata.md`，读取时 SHA-256 为
`6e23fed601cd491a7651ad607e0932760aa954f44abdb7f6a80a1e29c342dcd6`；其有限
诊断脚本 SHA-256 为
`00e6affbb3df4d3b05394ea9bcf39bab7520471afca7e6f73f94118da03591f9`。
本审计不产生 \(g_0\)、\(I\) 或新的 \(C_6\) 上界。

## 结论

**条件性上界接口、奇异层判定和 short-time whole-face 预算接受；数值上界仍未认证。**

`check.py` 是独立有限复核，实际用精确 SymPy 算术检查了：

* 指定六方向基的 \(|\det|=160\)；
* 全部 32 个一元、496 个二元、4960 个三元 body 子集的秩分别为 1、2、3；
* 每一个可能的 active body 方向都可由坐标翻转把上述基送到一个仍有
  \(|\det|=160\) 且含该方向的基；
* \(C_f=1024/(5\pi^6)\)、\(\epsilon^{7/6}\) 系数及成本代数消去。

该核验只证明存在 determinant-160 基，完全没有把浮点扫描当作全局最大
行列式证明。

## 奇异结构

令 \(S_2=\bigcup_{v,w\in V}\operatorname{span}\{v,w\}\)。对每个紧
\(K\Subset\mathbb R^6\setminus S_2\)，令其到有限并 \(S_2\) 的距离为
\(d>0\)，取 \(30\sqrt6\delta<d\)。如果所有 32 个一维 Cauchy 坐标中至多
两个有绝对值超过 \(\delta\)，余项范数最多 \(30\sqrt6\delta\)，总和不可能
落在 \(K\)。三个超过阈值的坐标对应三个不同 body 方向；独立秩审计表明它们
线性无关，并可补成六方向基。对该基的 fiber convolution 表示中，三个 active
密度至多 \(t/(\pi\delta^2)\)，其余三个至多 \((\pi t)^{-1}\)，余下 26 个变量
作为概率卷积。有限三元组求和给 \(p_t\) 在 \(K\times(0,1]\) 一致有界。
这证明短时可积；determinant-160 的 \(t^{-6}\) 长时界处理 \([1,\infty)\)，
从而 \(G\) 在 \(S_2\) 外有限且连续。

反向包含也成立：若 \(x=av+bw\) 落在一个独立二平面，选含 \(v,w\) 的
body 基并把另外 26 项限制在一个缩放为 \(t\) 的、在基坐标中有正概率的盒内。
当 \(a,b\ne0\) 时指定两项各提供 \(\asymp t\)、其余四个基坐标各提供
\(\asymp t^{-1}\)，所以 \(p_t(x)\gtrsim t^{-2}\)。零系数只增加发散阶；原点
由缩放和正性处理。因此 \(G=\infty\) 正好发生在 \(S_2\)。

审计要求把以上“fiber convolution 的事件分割”写成显式 Fubini/非负积分引理，
以免把不同补基的密度上界误读为同一坐标表示中的逐点乘法；这是写作完整性要求，
不是目前发现的反例。

## whole-face 预算与成本

在 \(x_1=1\) 上，32 个方向的第一坐标都为一，故 fiber 上至少一项
\(|z_v|\ge1/32\)。逐一选择含该 active 方向的 determinant-160 基，求和得到

\[
 \|p_t(1,\cdot)\|_\infty\le
 32\,{1024t/\pi\over160(\pi t)^5}
 ={1024\over5\pi^6}t^{-4}.
\]

第一坐标边际是 scale \(32t\) 的 Cauchy 密度，所以 whole face 的 \(L^1\)
范数为 \(32t/[\pi(1+32^2t^2)]\le32t/\pi\)。插值到 \(q=6/5\)，再利用后五个
坐标反射对称、Minkowski 和 \(\int_0^\epsilon t^{1/6}dt=(6/7)\epsilon^{7/6}\)，
严格得到

\[
 \left\|\int_0^\epsilon p_t(1,\cdot)dt\right\|_{L^{6/5}([0,1]^5)}
 \le {6\,2^{5/3}\over7\,5^{1/6}\pi^{11/6}}\epsilon^{7/6}.
\]

独立代码以符号算术核对了该常数。原稿中将未来 long-time 项记为
\(H_\epsilon\) 而未在本段重新定义；未来连续证书应明确它是何种
\(\int_\epsilon^\infty p_tdt\) 或同等余项，才可把此预算与其合并。

若未来严格证明 \(G\ge g_0>0\) a.e. 于 cube，并证明同一真实 Green 函数的
face 积分 \(I\) 上界，则 \(c=g_0/6\)、
\(K=(G-c)_+/[64(g_0-c)]\) 的质量严格化简为
\[
 \int K=6^{6/5}I/g_0^{6/5}.
\]
该式和一般 \(L^1\)、全实尺度接口代数正确：\(K\) 支配归一 cube 核，
Fourier 齐次性给 \(K_r*Lu=(\mu_{c,r}*u-u)/(64(g_0-c)r)\)，在同一接触集
补集非负；有理尺度合并零集、\(L^1\) dilation 连续性和单调截断依次给所有尺度
和一般 \(L^1\)。本审计接受其为**条件性**机制，不把模型成本或 (14) 记为
\(C_6\) 改进。

## 复现

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/audit/upper001/check.py
```

实际结果为 PASS；脚本 SHA-256 由其 JSON 输出为
`7839faf07804d788f8994447847e553e64f9c7ee60d31c4c901ca28ec3e3b664`。
