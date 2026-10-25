# I027：三个新的 R 下负行列式矩 cap

本批核对父任务提出的指数及Harris路径，三点均成立，得到同一真实形状概率 (R) 下的新矩上界：

| (gamma) | (Q_\gamma) 的显式公式值 | 可直接使用的严格有理cap |
|---|---:|---:|
| (1/5) | (1.8526921133850069613860141591\ldots) | (926346056693/500000000000) |
| (3/14) | (2.0606032208164426175209144845\ldots) | (2060603220817/10^{12}) |
| (1/4) | (3.0306921209739844685944279946\ldots) | (1515346060487/500000000000) |

严格含义为 (mathbb E_Rd^{-\gamma}\le Q_\gamma<)表中有理cap。新 (Q_{3/14}) 上端严格低于旧锁定 (2.4576305439372\ldots) 球的下端。**本批没有新 (I) 上界，没有majorant优化或 (C_5) 改善；(C_5\le32) 不变。**

只写独立clone的新 `upper_I027/`。主目录、其他代理及已锁020/021/022/024文件只读，未提交推送。输入hash见登记；没有五基/flag几何重算、旧Gamma cap公式重跑或随机/空间采样。

## 1. 固定概率与原始 detA cap

仍取16条未单位化body方向、独立标准 (alpha=3/2) stable时间、Lebesgue。正subordinator变量 (S_i) 独立且

\[
 \mathbb Ee^{-zS_i}=e^{-z^{\beta}},\quad\beta=3/4,
 \quad A=\sum_iS_iv_iv_i^T,\quad T=\sum_iS_i,\quad M=A/T.
\]

写 (d=\det M)。R用 (T^{-\beta}/Z_{\rm tilt}) 倾斜并推送至形状，

\[
 Z_{\rm tilt}=\mathbb ET^{-\beta}=[16\Gamma(7/4)]^{-1}.
\]

这里的 (Z_{\rm tilt}) 与024的flag代数常数 (Z_f=3.32689476018\ldots) 是两个不同对象，不能混用。

022/024的独立五因子分解，对任意 (gamma>0) 都给

\[
 \mathbb E(\det A)^{-\gamma}
 \le16^{-5\gamma}m(\gamma)^5Z_f^{-\gamma},\qquad
 m(\gamma)=\frac{\Gamma(1+\gamma/\beta)}{\Gamma(1+\gamma)}.
\]

原因是 (mathbb EF_j^{-\gamma}=m(\gamma)^5R_j^{-\gamma/\beta})，熵AM–GM/普通Hölder的log预算为

\[
 -5\gamma\log16+5\log m(\gamma)-\gamma\log Z_f
 +\gamma\mathrm{KL}(q\|q^*),
\]

其中同一个 (q_j^*=\theta_jR_j^{1/\beta}/Z_f) **不依赖 (gamma)**。所以不需要为三点重新做权重优化，更不需要重算几何。所有正stable负矩存在，故这个原始 (S) 概率下的cap适用于所有 (gamma>0)。

## 2. 转回 R 与 Harris 方向

精确变换为

\[
 \mathbb E_Rd^{-\gamma}
 =Z_{\rm tilt}^{-1}\mathbb E[T^k(\det A)^{-\gamma}],
 \qquad k=5\gamma-\beta.
\]

若 (0<k<\beta)，(T^k) 对每个 (S_i) 递增；((\det A)^{-\gamma}) 递减。后者可直接由Cauchy–Binet非负系数，或正定矩阵rank-one增量的det单调性证明。原独立产品概率的反向关联因此给**上界方向**

\[
 \mathbb E[T^k(\det A)^{-\gamma}]
 \le\mathbb ET^k\,\mathbb E(\det A)^{-\gamma}.
\]

不假定 (T/P) 独立，也不在 (R) 下错误使用产品概率。对两非负函数分别取 (min(N,\cdot))，先用有界两副本/逐坐标归纳，再单调取极限即可；右端有限将同时证明左端有限，不循环假设联合矩存在。

正矩存在条件必须保留：

\[
 \mathbb ET^k=16^{k/\beta}
 \frac{\Gamma(1-k/\beta)}{\Gamma(1-k)}<\infty
 \quad(0<k<\beta).
\]

因此本R-cap路径的有效指数区间是 (eta/5\le\gamma<2\beta/5)，即 (3/20\le\gamma<3/10)；左端 (k=0) 为精确常数项，不必使用Harris。

若 (gamma<3/20)，两个函数同向递减，不能沿用这里的Harris上界；若 (gamma\ge3/10)，所用正矩发散或Gamma因子触及极点。原始detA负矩cap虽仍成立，不代表这个R-cap公式可以无条件外推。

## 3. 尺度指数及最终公式

合并 (Z_{\rm tilt}^{-1}=16\Gamma(7/4))、正矩和原始detA cap，16的净幂为

\[
 1+k/\beta-5\gamma
 =5\gamma(1/\beta-1)=5\gamma/3.
\]

最终

\[
 \boxed{\mathbb E_Rd^{-\gamma}\le Q_\gamma=
 \Gamma(7/4)16^{5\gamma/3}
 \frac{\Gamma(1-k/\beta)}{\Gamma(1-k)}
 m(\gamma)^5Z_f^{-\gamma}.}
\]

三点的准确参数为

| (gamma) | (k=5\gamma-\beta) | (k/\beta) | 16净幂 |
|---|---|---|---|
| (1/5) | (1/4) | (1/3) | (1/3) |
| (3/14) | (9/28) | (3/7) | (5/14) |
| (1/4) | (1/2) | (2/3) | (5/12) |

本批仅按预登记这三点核验，未扩展Gamma网格。在 (\gamma=3/20) 的符号极限，它还原024的J公式；这只是解析一致性检查，没有数值重跑024。

## 4. 实际一次有限计算

读取024已经认证的 (Z_f) 精确有理区间，不重新构造代数根或几何。一次384-bit Arb批使用未移位Mellin式

\[
 m(\gamma)=\frac{\Gamma(\gamma/\beta)}{\beta\Gamma(\gamma)},
\]

与题设移位Gamma式由递推恒等。对正 (gamma)，用 (Z_f) 上端取得Q下包围、用其下端取得Q上包围；全部Gamma与幂运算由Arb向外包围。保存球/端点的50位显示区间为

\[
 Q_{1/5}\in[
 1.85269211338500696138601415910284796776225617324733,
 1.85269211338500696138601415910284796776225617324734],
\]
\[
 Q_{3/14}\in[
 2.06060322081644261752091448451310557804816335235458,
 2.06060322081644261752091448451310557804816335235459],
\]
\[
 Q_{1/4}\in[
 3.03069212097398446859442799463492845657412770031681,
 3.03069212097398446859442799463492845657412770031682].
\]

有理caps由这些新球的严格向上端点取得。仅对同指数 (3/14) 比较旧完整Q球；旧表达式没有重新计算，且没有把不同指数cap互相当成改进比较。

复现：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_I027/negative_moment_caps.py
```

实际退出码0，约0.00087秒；python-flint0.9.0、Arb384位、Fraction准确参数与端点。完整有理区间、旧Q锁定比较、源hash在 `results.json`。没有MATLAB、后台进程、图窗、采样或参数优化。

## 5. 与007排除的关系及边界

007准确排除针对的是旧输入 (E_Rd=m,E_Rd^2=s,E_Rd^{-3/14}\le Q_{\rm old}) 形成的上包络预算族。其双点分布论证需要该分布满足相应负矩约束；换成严格更小的新cap后，原条件必须重新检验，不能自动把旧排除迁移成新cap下的排除。

反过来，负矩cap下降本身也不证明某个标量majorant现在能达到目标，更不证明实际 (I) 已下降。本批遵照范围不启动标量优化，不计算或重放旧抽象分布。父任务可根据三项新真实矩输入另行登记；本报告只认证上述cap公式及有限算术。

必要白名单仅报告、登记、脚本、结果、manifest五项，无输入快照；不修改主clone、不提交推送。新I/C5结论均为空，正式 (C_5\le32) 保持。
