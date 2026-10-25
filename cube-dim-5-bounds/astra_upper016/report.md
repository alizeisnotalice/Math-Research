# Astra016：两种 Cauchy ridge 幂筛选与四次行列式下界的封顶

本批只考察新候选 `Z=a*u^m`，`u=(W+2)/12`，`m=2,3`，固定 `alpha=1`、`beta=1/2`。第一阶段的必要角点上限分别为 `3.402075122...e-4`、`6.806103486...e-4`，均不能排除014的所需门槛 `1.519836298...e-4`。按预登记只继续较大者 `m=3`。

第二阶段新算到26阶矩，构成唯一双节点 Hermite 下界。其实际角点值仅 `3.367169648496...e-5`，未达门槛。进一步以同一准确矩证明：**所有**形如 `ell(d)=d Q3(d)<=sqrt(d)` 的三次 `Q3`，配合此 `m=3` 和任意 `a>=0`，角点必要上限均满足
\[
 F_{\rm cubic}\le3.3671696499665860062\ldots\,10^{-5}
 <\boxed{337/10000000}.
\]
差距不是节点有理化或浮点选择造成。依登记不运行空间矩、固定 `a` Bernstein 或自适应 `a(x)` 证书。没有新全域 `g`，没有新 `C5` 上界。此结论没有排除真实 `sqrt d` 的 `m=3` 家族、其它行列式下界，也没有对 `m=2` 做第二阶段结论。

## 1. 模型、唯一形状与第一阶段必要筛选

`R_(1/2), M, d=det M, B=x'adj(M)x, W=sum_(i<j)Mij xi xj` 均与014/015相同。015证明在整个 cube 上 `0<=u<=1`。Young 的准确形式是
\[
 d^{3/2}B^{-2}\ge3\sqrt d\,Z^2-2BZ^3.
\]
在唯一必要角点 `x=(1,1,1,1,1)` 定义
\[
 E_m=\mathbb Eu^{2m},\quad A_m=\mathbb Edu^{2m},\quad
 V_m=\mathbb EBu^{3m},\quad U_m=\mathbb E\sqrt d\,u^{2m}.
\]
Cauchy–Schwarz 给 `U_m²<=A_m E_m`；准确优化 `a` 得
\[
 \max_a c_1(3a^2U_m-2a^3V_m)
 =c_1U_m^3/V_m^2
 \le c_1(A_mE_m)^{3/2}/V_m^2,\quad c_1=1/(32\pi^3). \tag{1}
\]
第一批新算 `EW^n,EdW^n` 至 `n=6`、`EBW^n` 至 `n=9`。只压缩一个不带空间变量的、在 S5 下不变的行列式或伴随多项式因子；`W^n` 仍完整按匹配 vertex parity 枚举，未独立压缩两个因子而丢失相对位置。

准确结果：

| m | `E_m` | `A_m` | `V_m` | 式 (1) 上限 |
|---:|---:|---:|---:|---:|
| 2 | `253/18144` | `553/1181466` | `2514503/11313718416` | `0.00034020751220070669...` |
| 3 | `112733/14370048` | `3231145/33941155248` | `965344619/31225862828160` | `0.00068061034863345913...` |

384-bit Arb 及有理外端分别严格大于014门槛的有理上端。这个“未被排除”不是可达性；这里的 Cauchy–Schwarz 很松。第一批最高图阶13，缓存由015的29个状态复制到新文件，新增至880状态，约0.087秒；低阶输出与015全数准确匹配。

## 2. 仅 m=3 的高阶准确矩和唯一 Hermite 选择

记正测度 `dnu=d*u^6 dR`，
\[
 \mu_j=\int t^j\,d\nu(t)=\mathbb Ed^{j+1}u^6,\qquad j=0,1,2,3.
\]
`hermite_moments.cpp` 新生产 `E[d^k W^n]`，`k=1..4,n=0..6`。行列式各幂的原始项数为73、2021、25050、171705，S5轨道项数为7、56、384、2001；与旧 alpha=3/2 共享的仅是整数代数表达式，所有图矩均来自独立 beta=1/2 后端及缓存。

最高图阶26；后端仍保留十边各5位、总阶31守卫，全部累计值为 GMP 任意精度。新 beta 图状态累计55942个，实际运行70.8352秒，在120秒 checkpoint 内完成。每一完整 `W^n` 保存结果；若中断，图缓存持久化，角点和从空重建避免重复加。`k=1` 的全部7个输出与第一阶段准确匹配。没有加载 beta=3/4 的矩值。

四个 `mu_j` 的完整有理数见 `hermite_verification.json`。正方差
\[
 \Delta_\mu=\mu_0\mu_2-\mu_1^2>0
\]
给两节点 Gauss 多项式
\[
 p_2(t)=t^2-A t+B_0,\quad
 A=\frac{\mu_0\mu_3-\mu_1\mu_2}{\Delta_\mu},\quad
 B_0=\frac{\mu_1\mu_3-\mu_2^2}{\Delta_\mu}.
\]
其节点准确区间约为 `.07030184144986368`、`.39970021160746826`。按计算前登记的最近 `1/10000` 规则，唯一有理节点取
\[
 r_1=703/10000,\qquad r_2=3997/10000.
\]
令 `h(t)=t^-1/2`。其四阶导数 `105 t^-9/2 /16` 为正，所以在两节点取值与导数的三次 Hermite 插值 `Q_H` 满足
\[
 h(t)-Q_H(t)=\frac{h^{(4)}(\xi)}{4!}(t-r_1)^2(t-r_2)^2\ge0
 \quad(0<t\le1).
\]
这即使 `t` 不在两节点之间也成立，应用区间取包含节点与 `t` 的闭区间，远离0。乘 `t` 给 `t Q_H(t)<=sqrt t`。直接插值 `sqrt t` 会因四阶导数为负而给上界，本批没有那样做。

由于 `h'(r)=-h(r)/(2r)`，每个 Hermite 系数是 `h(r1),h(r2)` 的准确有理线性组合。程序先用分母 `10^14` 的正有理端点夹住各 `1/sqrt r`，并以平方比较独立核验；再依线性组合系数的符号选下端。由 `t^k>=0`，所有系数分别向下保证
\[
 \ell(t)=t\sum_{j=0}^3 c_jt^j\le tQ_H(t)\le\sqrt t.
\]
`c_j` 可为负，均准确保留；不要求这个 minorant 点态非负。

于是正项下界 `U_L=sum c_j mu_j` 为正。按登记取
\[
 a=\lfloor10^5 U_L/V_3\rfloor/10^5=5131/500.
\]
用014已证明的 `c1` 有理下端，角点认证表达式为
\[
 c_L(3a^2U_L-2a^3V_3)
 =3.3671696484964373975\ldots\,10^{-5}.
\]
这是对候选下界在这个点的准确计算；角点真实 Green 本身为无穷，不能把候选值误当成真实点值，也不能将它称为全域 `g`。

## 3. 全部四次行列式 minorant 的严格最优上限

该批四个准确矩足以证明有理化不是瓶颈。设准确 Gauss 节点为 `z1,z2`，权重为
\[
 \lambda_1=\frac{\mu_0z_2-\mu_1}{z_2-z_1},\qquad
 \lambda_2=\frac{\mu_1-\mu_0z_1}{z_2-z_1}.
\]
程序严格核 `0<z1<z2<1` 与 `lambda_i>0`。这两个权重匹配次数0、1的矩；节点满足二次方程，且 `p2` 与1、t正交，所以二次、三次矩也准确匹配。因此对任意次数至多3的 `Q`，有精确恒等式
\[
 \int Q(t)d\nu(t)=\lambda_1Q(z_1)+\lambda_2Q(z_2).
\]
如果 `tQ(t)<=sqrt t` 对 `(0,1]` 全局成立，则 `Q(zi)<=zi^-1/2`，故
\[
 \mathbb E[dQ(d)u^6]\le
 U_G:=\lambda_1z_1^{-1/2}+\lambda_2z_2^{-1/2}. \tag{2}
\]
准确 Gauss 节点的 Hermite 插值同时是全局下界，并在两节点取等，故达到式 (2)。这证明 `U_G` 正是整个该阶下界族的最优值，包含任意 Hermite 节点或其它三次 minorant，并非只对本批选定节点成立。

再准确优化 `a`，384-bit Arb 严格核
\[
 c_1 U_G^3/V_3^2
 =3.3671696499665860062\ldots\,10^{-5}
 <337/10000000.
\]
这个上限与实际有理候选之差约 `1.47e-14`；调整节点舍入或 `a` 精度不可能弥补与014门槛的差距。

自适应 `a(x)=U(x)/V(x)` 仍受这个必要角点上限限制。候选多项式期望在角点附近连续，`V_3>0`；严格差距扩展到正测邻域。因此采用自适应 `a(x)` 本身不能挽救本批 **m3 加 dQ3 minorant** 的全域目标。与独立 ridge 下界作分区覆盖属于另一个机制，不在这个排除结论内。

## 4. 状态边界与复现

第一阶段对真正 `sqrt d` 的两个家族均未排除；第二阶段只排除 `m=3` 的 `d*三次多项式` 下界能力。`m=2` 没有做 Hermite 生产或判断；没有更高 `m` 扫描、没有空间网格、没有高于31阶的图，没有 `C5` 新结论。完整一般 `L1`、全尺度与五维 `S1/S2` 传递仍须独立闭合后才能使用任何新的全域 Cauchy `g`。

```sh
clang++ -std=c++17 -O2 -Wall -Wextra -Werror -isystem /opt/homebrew/Cellar/gmp/6.3.0/include -L/opt/homebrew/Cellar/gmp/6.3.0/lib cube-dim-5-bounds/astra_upper016/corner_screen.cpp -lgmpxx -lgmp -o /tmp/astra_upper016_corner_screen
/tmp/astra_upper016_corner_screen cube-dim-5-bounds/astra_upper016/corner_screen.json /tmp/astra_upper016_beta_half_cache.tsv
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper016/verify_screen.py
clang++ -std=c++17 -O2 -Wall -Wextra -Werror -isystem /opt/homebrew/Cellar/gmp/6.3.0/include -L/opt/homebrew/Cellar/gmp/6.3.0/lib cube-dim-5-bounds/astra_upper016/hermite_moments.cpp -lgmpxx -lgmp -o /tmp/astra_upper016_hermite_moments
/tmp/astra_upper016_hermite_moments cube-dim-5-bounds/astra_upper016/hermite_moments.json /tmp/astra_upper016_beta_half_cache.tsv
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper016/verify_hermite.py
```

两份 C++ 均以 Apple clang21 的 `-Wall -Wextra -Werror` 通过。浮点只用于建议有理节点；节点舍入区间、函数根夹和最终全部常数判断另用有理/Arb核验。全部准确矩仍共享015的已锁 GMP 后端；独立多变量 Stieltjes及单边均匀交叉见015，016不声称高阶矩全数由第二后端重产。

必要归档是两个登记、两批 C++ 与小型准确矩结果、两份验证代码及验证结果、本报告和 manifest。大图缓存留 `/tmp`，记录源哈希和实际文件哈希，不归档。未提交，未改变旧锁文件。
