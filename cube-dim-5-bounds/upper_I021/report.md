# I021：三子集二维商的两层线性stable分组

本批推导新的完整辅助界

\[
 I\le1.335145794635596923586034684\ldots\times10^{-5}
 <\boxed{\frac{6675729}{500000000000}}.
\]

它严格改善已锁I020的 (1.4214700731529485\ldots\times10^{-5})，仍未达到 (11/10^6)。新机制覆盖完整Green、全部正时间、整个正unit face；不使用随机诊断、时间截断或空间采样。独立审计前不替换主辅助账本；没有新 (g) 或 (C_5<32)。

仅写本独立clone的 `upper_I021/`；I020已锁且未修改，主和其他代理文件只读。不重复主type权重优化。有限登记在 `registration.json`。

## 1. 沿用接口及新分组对象

固定16-body、权1、独立 (alpha=3/2) stable、Lebesgue。令 (eta=3/4)、(b=3/20)、(S_i) 独立正stable，(A=\sum_iS_iv_iv_i^T)，(W=16^5)。已审019及020给

\[
 J=\mathbb E_R(\det M)^{-b}=Z^{-1}\mathbb E(\det A)^{-b},\quad
 Z=\frac1{16\Gamma(7/4)},
\]
\[
 I\le c^{10/7}\frac{\pi^2}{60}J^{10/7},\quad
 c=\frac3{2^{15/2}\pi^{5/2}}.
\]

这里不假定 (T=\sum S_i) 与 (P_i=S_i/T) 独立。总体尺度抵消仍由 (5b=\beta) 精确取得。

对一个三子集 (J\subset[16])，其方向总独立（程序对560个Gram行列式均精确核正）。正交投影至 (\operatorname{span}V_J) 的二维正交补，记投影 (q_i)。以下系数不依赖所选二维坐标：

\[
 w_{ik}=\det(V_{J\cup\{i,k\}})^2
 =D_J\det(q_i,q_k)^2,\quad
 D_J=\det(V_J^TV_J)>0.
\]

写 (B_J=\sum_{i<k}w_{ik}S_iS_k)。Cauchy–Binet给

\[
 \det A=\frac1{10}\sum_{|J|=3}\left(\prod_{j\in J}S_j\right)B_J.
\]

因为每个五方向基的三子集有 (inom53=10) 个。与020单线性完成项不同，这里保留完整二维商的正二次结构。

## 2. 共线类合并与两个独立stable线性和

只对 (q_i\ne0) 分类，将平行/反平行向量放入同一无向线类 (ell)。选代表 (r\in\ell)，记

\[
 a_i=\frac{\|q_i\|^2}{\|q_r\|^2}\quad(i\in\ell),\qquad
 X_\ell=\sum_{i\in\ell}a_iS_i,\quad
 L_\ell=\sum_{k\notin\ell}w_{rk}S_k.
\]

若 (q_i\parallel q_r)，则 (w_{ik}=a_iw_{rk})。因此

\[
 B_J=\frac12\sum_\ell X_\ell L_\ell.
\]

每个非零完成对的两条商线各计一次。(X_\ell,L_\ell) 的变量支持不交，并且二者都不含选中的三变量 (S_j)，所以这五个因子在原产品概率下相互独立。

设 (x_\ell=\sum_{i\in\ell}a_i)、(l_\ell=\sum_kw_{rk})，

\[
 C_{J\ell}=x_\ell l_\ell
 =\sum_{i\in\ell,k\notin\ell}w_{ik}>0,
 \quad\theta_{J\ell}=\frac{C_{J\ell}}{20W}.
\]

于是精确身份为

\[
 \det A=\frac1{20}\sum_{J,\ell}\left(\prod_{j\in J}S_j\right)X_\ell L_\ell,
 \quad\sum_{J,\ell}C_{J\ell}=20W.
\]

令

\[
 F_{J\ell}=\left(\prod_{j\in J}S_j\right)
 \frac{X_\ell}{x_\ell}\frac{L_\ell}{l_\ell}.
\]

则 (det A/W=\sum\theta_{J\ell}F_{J\ell})，权和1；AM–GM给 ((\det A)^{-b}\le W^{-b}\prod F_{J\ell}^{-b\theta_{J\ell}})。当 (S_i=1) 时所有 (F_{J\ell}=1)，故规范化在均匀源处取等。

## 3. 真实scope、分数Hölder与固定Gamma式

把两个线性和的正归一系数分别记作概率向量 (p^X,p^L)，

\[
 R_X=\sum_i(p_i^X)^{3/4},\quad R_L=\sum_i(p_i^L)^{3/4}.
\]

正stable的精确线性和尺度及五因子独立性给，对 (s>0)，

\[
 \mathbb EF_{J\ell}^{-s}=m(s)^5(R_XR_L)^{-s/\beta},
 \quad m(s)=\frac{\Gamma(1+s/\beta)}{\Gamma(1+s)}.
\]

真实scope为 (A_{J\ell}=J\cup\{i:q_i\ne0\})，不依赖 (ell)；某些额外落入三维span的方向完全不进入该组。精确有限几何及方向传递性给每个坐标的占用

\[
 \sum_{(J,\ell):,i\in A_{J\ell}}\theta_{J\ell}
 =\rho=\frac{247}{256}.
\]

因此 (h_{J\ell}=\theta_{J\ell}/\rho) 是产品概率超图的fractional matching。与020同一逐坐标Hölder证明应用于 (F_{J\ell}^{-b\rho})，得

\[
 \mathbb E(\det A)^{-b}
 \le W^{-b}m(b\rho)^{5/\rho}
 \prod_{J,\ell}(R_XR_L)^{-\theta_{J\ell}/5}.
\]

所有组负矩有限，足以直接对非负截断取极限，不对退化形状施加虚假的统一条件数。令

\[
 D_{21}=\prod_{J,\ell}(R_XR_L)^{-\theta_{J\ell}/5},
\]

得到闭合矩界

\[
 \boxed{J\le J_{21}=2\Gamma(7/4)
 \left[\frac{\Gamma(1+247/1280)}{\Gamma(1+741/5120)}\right]^{1280/247}D_{21}}.
\]

注意 (b\rho=741/5120)、(b\rho/\beta=247/1280)。该式仍是完整尺度抵消后的原独立 (S) 期望；没有在倾斜 (R) 下错误地使用独立性。

## 4. 精确实现与实际运行

没有扫五方向基。对每个三子集，计算 (G=V_J^TV_J)、(D=\det G)、(z_i=V_J^Tv_i)，并用整数

\[
 A_{ik}=D(v_i\cdot v_k)-z_i^T\operatorname{adj}(G)z_k
 =D(q_i\cdot q_k).
\]

完成系数直接为

\[
 w_{ik}=\frac{A_{ii}A_{kk}-A_{ik}^2}{D}\in\mathbb Z_{\ge0}.
\]

程序断言整除、对称性和非负性；(A_{ii}=0) 判定零商向量，(w_{ik}=0) 判共线。还逐组核 (w_{ik}A_{rr}=A_{ii}w_{rk})，确保合并系数正确。一个选定三子集的全部120个完成对另用独立5阶Fraction行消元复核。

有限560个三子集给2240个两层组，规范化概率仅5个类型。逐三子集核 (sum_{i<k}w_{ik}=D\cdot16^2)，逐线类核 (C_{J\ell}=x_\ell l_\ell)，总组质量 (20W=20,971,520)，16个scope占用都严格为247/256。完整类型概率、代表与权重在 `results.json`；没有保存重复全量几何副本。

本机一次固定384-bit Arb实际运行约0.094秒，退出码0，得

\[
 J_{21}=1.445587993610619085725373267\ldots,
 \qquad I_{21}=1.335145794635596923586034684\ldots\times10^{-5}.
\]

Arb严格核出 (I_{21}) 小于020旧球的下端；有理上包围取 (6675729/500000000000)。Gamma只算本批固定公式，没有重跑019/020基线，也没有做type权重优化、空间值、Monte Carlo或QMC。

复现：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_I021/two_level_negative_moment.py
```

必要文件仅报告、登记、脚本、结果、manifest五项，无输入快照；源hash在登记与结果中固定。未调用MATLAB，未提交、未推送。数学证据是全域分组与Hölder证明，有限工具证据仅核精确几何、规范化和明确常数。本批锁定后，二子集/三线性因子另入022，不覆盖021。
