# Astra015：唯一 Cauchy Young 候选的严格必要上限排除

固定 `alpha=1`、`beta=1/2`、同一 16 条未单位化 body 方向。唯一候选
\[
 Z=a u,\qquad u=(W+2)/12,\qquad a\ge0
\]
在必要角点的全部 `a`、全部 `sqrt(d)` 下逼近的最优值，严格不超过
\[
 0.000104704438710474493840\ldots <\boxed{21/200000}=0.000105.
\]
这低于 Astra014 当前完整积分预算要求的 `0.00015198362980537...`。因此按照预登记第一阶段即停止：**未运行**高次 `d` 矩、Hermite、81项空间证书。此为该固定 Young 家族的障碍，不是真实 Green 的上界，不排除 Cauchy 模型，不产生新 `g` 或 `C5` 数值。

## 1. 固定模型及 Young 方向

沿用014的 `R_(1/2)`，定义
\[
 M=\sum_vP_vvv^T,\quad d=\det M,\quad
 B(x)=x^T\operatorname{adj}(M)x,\quad
 W(x)=\sum_{i<j}M_{ij}x_ix_j.
\]
对角元 `Mii=1`，`0<d<=1` 几乎处处。若 `x` 是 cube 顶点，则
\[
 x^TMx=\sum_vP_v(v\cdot x)^2\ge1
\]
因为五个符号之和为奇整数。因此 `W=(x'Mx-5)/2>=-2`。`W(x)` 对各 `xi` 分别仿射，整个 cube 的最小值可取在顶点，故全域 `W>=-2`；同理 `W<=10`，所以 `0<=u<=1`。

对 `y=sqrt(d)>0`、`B>0`，单变量极值给
\[
 \frac{y^3}{B^2}=\max_{Z\ge0}(3yZ^2-2BZ^3).
\]
因此对任意 `a>=0`，
\[
 H=d^{3/2}B^{-2}\ge3a^2\sqrt d\,u^2-2a^3Bu^3. \tag{1}
\]
式 (1) 保留形状与 `B` 的相关性。原点按 `G(0)=infinity` 处理，其余点 `B>0`。

在唯一角点 `x=(1,1,1,1,1)`，令
\[
 U=\mathbb E_R\sqrt d\,u^2,\quad V=\mathbb E_RBu^3,
 \quad A=\mathbb E_Rdu^2,\quad E=\mathbb E_Ru^2.
\]
由 Cauchy–Schwarz，`U<=sqrt(AE)`，而 `V>0`。准确优化 `a` 给
\[
 \max_{a\ge0}(3a^2U-2a^3V)=U^3/V^2\le(AE)^{3/2}/V^2. \tag{2}
\]
若用任何 `ell(d)<=sqrt(d)` 替换正项，只会降低式 (2)。故这个必要上限同时排除所有此种 Hermite/多项式改善。

## 2. 新 beta 的准确矩后端

014 的倾斜形状变换给，在原点邻域
\[
 \mathbb E_R(1-\sum_e h_er_e)^{-1/2}
 =\frac{16}{\sum_{v\in V}(1-\sum_eh_ev_iv_j)^{1/2}},
 \qquad r_{ij}=M_{ij}. \tag{3}
\]
本批只重用已验证后端的图编码、GMP 和 checkpoint 实现，**不读取任何 beta=3/4 数值缓存**。新的常数为
\[
 D_n=\prod_{j=0}^{n-1}(2j+1),\qquad
 C_n=\prod_{j=1}^{n-1}(2j-1),\quad D_0=C_1=1.
\]
设 `k` 是十边多重指数，`n=|k|`，`B_k=D_n E r^k`。对 (3) 的分母乘积求导，若图不是 Euler 图则 `B_k=0`；否则
\[
 B_0=1,\qquad
 B_k=\sum_{0<j\le k\atop j\text{ Euler}}
 \binom{k}{j}C_{|j|}B_{k-j}. \tag{4}
\]
这是对 `q^n` 缩放导数的整数递推，`q=2`；分母乘积的非零导数为 `-C_n/2^n`。因此全部图矩整数算术精确，符号不从浮点猜测。

新头文件 `graph_moments_beta_half.hpp` SHA256 为
`bfdf502cac067d2a1d9465c96ac15076c2e897cee89f8121bfed40977fc3fbf0`。
它的 journal 协议明确含 `astra015-K5-beta1/2` 与新源哈希，旧 beta 缓存不能被接受。十条边各5位，保留总阶31和高位合法性守卫；本批实际最高阶仅7，所有 `B,C,D`、二项式乘积和空间累加系数均为 GMP 任意精度，`uint64` 只编码图。

`corner_screen.cpp` 以 Leibniz 式生成 `det(M)` 及全部伴随项，在角点准确相加成 `B`，再乘完整 `W^n`；未把 `B` 矩替换成边际均值乘积。带符号多项式在准确整数中累加，再除各总阶自己的 `D_n`。底层缓存新建于 `/tmp/astra_upper015_beta_half_cache.tsv`，首批加载数为0，实际新算28个非零状态，加上常数共29个状态；运行约0.0174秒，未触及120秒 checkpoint。若触及，缓存安全保存，角点累加从空重启。

## 3. 独立交叉和实际矩

每条边把16个方向分成8正、8负，故其 Stieltjes 变换为
\[
 \frac2{\sqrt{1+z}+\sqrt{1-z}}
 =\frac12\int_{-1}^{1}(1+zr)^{-1/2}\,dr.
\]
紧支撑矩唯一性给 `r` 均匀分布于 `[-1,1]`，从而 `E r^(2j)=1/(2j+1)`，奇矩为0。程序准确检查了次数0至6。

另由16个实际符号字符，直接构造多变量普通 Taylor 系数，再以有限几何级数倒数计算 (3)，不调用 GMP 图递推。它逐一交叉三角形、四环、邻双边和不相邻双边：

| 图 | 准确矩 |
|---|---:|
| `r01² r02² r12²` | `53/495` |
| `r01 r03 r12 r23` | `1/7` |
| `r01² r02²` | `17/105` |
| `r01² r23²` | `17/105` |

共11个独立交叉目标全部通过。这个交叉共享数学变换 (3)，但普通级数路径独立于 Euler 整数后端。角点最终准确结果为
\[
 (\mathbb EW^n)_{n=0}^2=(1,0,10/3),
\]
\[
 (\mathbb EdW^n)_{n=0}^2=(8/63,0,1520/9009),
\]
\[
 (\mathbb EBW^n)_{n=0}^3=(8/7,-32/63,80/63,-32/91).
\]
这些矩全部属于新的 `beta=1/2`。代入 `u=(W+2)/12` 得
\[
 E=11/216,\qquad A=127/27027,\qquad
 V=44/7371,\qquad AE=127/530712. \tag{5}
\]

## 4. 严格排除及解释范围

`c1=1/(32pi^3)`。由 (2)、(5) 得所有 `a` 的必要角点上限
\[
 F_*\le\frac1{32\pi^3}
 \frac{(127/530712)^{3/2}}{(44/7371)^2}
 =0.000104704438710474493840\ldots.
\]
独立的纯有理比较避免从此小数重建证书。Machin 恒等式
`pi=16 atan(1/5)-4 atan(1/239)` 配合第一项上界及六项交错下界，准确证明 `pi>333/106`。因此
\[
 F_*^2<\frac{4131500149367901037}
 {376797974936010039515676672}
 <(21/200000)^2.
\]
`21/200000` 又严格小于014所需 `g` 的 Arb 有理下端 `37995907451/250000000000000`。由于所有量正，平方比较方向准确。

角点本身是 Cauchy Green 的无限值点，这不削弱这个**候选证书**的必要障碍：式 (1) 的期望函数 `U(x),V(x)` 是有界形状量的连续函数，且这里 `V>0`。优化后的候选上限也在角点邻域连续。角点上限与目标严格分离，所以其邻近正测面片也不能由这一族单独认证。此论证没有把有限候选值冒充真实 Green 的值；若后续将其它独立 ridge 下界与本族覆盖组合，需另作新的分区证明，这不在本排除范围内。

本批上限略高于 `1e-4`，故**没有**排除该族给出某个约 `1e-4` 的辅助 `g`，更没有排除改变 `u`、改变 `Z` 的行列式依赖、不同 alpha 或更紧 `I` 后的方案。只是当前固定 `u` 与014积分预算的突破目标不可能通过继续提高 `sqrt d` 下逼近实现。

## 5. 停止记录、复现与白名单

Stage1 已满足预登记停止条件，Stage2/3 未启动。未产生高阶角点、Hermite节点、全空间矩或 Bernstein 证书。因没有新的 `g`，没有把 alpha=3/2 的全空间传递文字机械替换参数；未来任何 Cauchy 成本主张仍须核五维 `S1/S2`、一般 `L1` 与全部尺度的正部有限源接口。

```sh
clang++ -std=c++17 -O2 -Wall -Wextra -Werror -isystem /opt/homebrew/Cellar/gmp/6.3.0/include -L/opt/homebrew/Cellar/gmp/6.3.0/lib cube-dim-5-bounds/astra_upper015/corner_screen.cpp -lgmpxx -lgmp -o /tmp/astra_upper015_corner_screen
/tmp/astra_upper015_corner_screen cube-dim-5-bounds/astra_upper015/corner_screen.json /tmp/astra_upper015_beta_half_cache.tsv
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper015/verify_screen.py
```

白名单为登记、新后端、角点源码与准确结果、独立验证源码与结果、本报告和 manifest。29状态缓存约数KB、可由源码在不足一秒内再生，不归档；manifest 保留其实际哈希。未提交，未修改014或任何他人目录。主 C5 上界32不变。
