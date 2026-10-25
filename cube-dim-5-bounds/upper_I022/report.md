# I022：二子集的三线性stable flag 完整上界

本批得到

\[
 I\le1.302214692587240555056040203\ldots\times10^{-5}
 <\boxed{\frac{13022147}{10^{12}}},
\]

严格改善I021的 (1.33514579463559\ldots\times10^{-5})，仍高于 (1.18\times10^{-5}) 及目标 (1.1\times10^{-5})。这是完整时间/完整正unit face的新辅助 (I) 机制，待主独审；没有新 (g) 或 (C_5<32)。仅写本独立clone的 `upper_I022/`，不修改已锁020/021或主/他人文件，不进行主负责的type权重优化。

## 1. 固定模型及成对完成多项式

模型与019/020/021相同：16条未单位化body方向，权1，独立 (alpha=3/2) stable、Lebesgue，(S_i) 为正 (eta=3/4) stable。记 (A=\sum S_iv_iv_i^T)、(W=16^5)、(b=3/20)，

\[
 J=\mathbb E_R(\det M)^{-b}=Z^{-1}\mathbb E(\det A)^{-b},
 \quad Z=\frac1{16\Gamma(7/4)},
\]
\[
 I\le c^{10/7}\frac{\pi^2}{60}J^{10/7},\qquad
 c=\frac3{2^{15/2}\pi^{5/2}}.
\]

对无序两子集 (J=\{i,j\})，令

\[
 B_J=\sum_{|K|=3,\ K\cap J=\varnothing}
 \det(V_{J\cup K})^2\prod_{k\in K}S_k.
\]

每五方向基含10个两子集，故 (det A=(1/10)\sum_{|J|=2}(\prod_JS)B_J)。其正交商维数为3。

## 2. 两次共线类合并

在三维商中记投影为 (q_i)。选非零商线类 (X) 与代表 (r)，令

\[
 a_i=\|q_i\|^2/\|q_r\|^2,
 \qquad X_S=\sum_{i\in X}a_iS_i.
\]

固定 (r) 后，再投影至 (q_r^\perp) 的二维商。残余商线类记作 (Y)，代表 (s)，令

\[
 b_j=\|q'_j\|^2/\|q'_s\|^2,\qquad
 Y_S=\sum_{j\in Y}b_jS_j,
\]
\[
 L_S=\sum_{k\notin\operatorname{span}(J,r,s)}
 \det(V_{J,r,s,k})^2S_k.
\]

第一商的平行系数使一次完成多项式按3个完成方向计数；二维商再按2个方向计数。因此精确身份为

\[
 B_J=\frac16\sum_{X,Y}X_SY_SL_S,
 \qquad\boxed{\det A=\frac1{60}\sum_{J,X,Y}
  \left(\prod_{j\in J}S_j\right)X_SY_SL_S.}
\]

所有平方行列式都不受方向代表的整体符号影响。三线性因子的支持依次为 (X)、(Y)、剩余最后完成方向，与选中的两个 (S_j) 不交。任何两条原body方向的线性span恰含这两条方向（本批也逐对精确核验），所以这五个因子的支持**覆盖全部16个变量**。

设 (x=\sum a_i)、(y=\sum b_j)、(l=\sum\det(V_{J,r,s,k})^2)，

\[
 C_f=xyl,\quad\theta_f=\frac{C_f}{60W},\quad
 F_f=\left(\prod_{j\in J}S_j\right)(X_S/x)(Y_S/y)(L_S/l).
\]

则 (sum_fC_f=60W)、(det A/W=\sum_f\theta_fF_f)。均匀 (S_i=1) 时所有 (F_f=1)，AM–GM取等。

## 3. 完整负矩公式与真实scope成本

三个正归一概率向量分别为 (p^X,p^Y,p^L)。记

\[
 R_f=\prod_{H=X,Y,L}\sum_i(p_i^H)^{3/4},\qquad
 m(s)=\frac{\Gamma(1+s/\beta)}{\Gamma(1+s)}.
\]

五因子在原 (S) 产品概率下独立，因此 (mathbb EF_f^{-s}=m(s)^5R_f^{-s/\beta})。

此时每组真实scope已是全部16变量，分数Hölder占用 (ho=1)，不能沿用020/021的 (ho<1) 损失折扣。但完整AM–GM及普通加权Hölder仍给

\[
 \mathbb E(\det A)^{-b}
 \le W^{-b}m(b)^5\prod_fR_f^{-\theta_f/5},
\]
\[
 \boxed{J\le J_{22}=2\Gamma(7/4)
 \left[\frac{\Gamma(6/5)}{\Gamma(23/20)}\right]^5D_{22}},
 \qquad D_{22}=\prod_fR_f^{-\theta_f/5}.
\]

负矩有限，积分Hölder由非负截断与单调收敛合法推广。没有在形状概率 (R) 下错误使用独立性。尽管 (ho=1)，本批实际全式仍改善021，不能只根据scope预判成败。

## 4. 无五基扫描的精确投影证书

对每一对 (J)，写其2阶Gram矩阵 (G)、(D=\det G\in\{16,24\})，(z_i=V_J^Tv_i)。用整数

\[
 A_{ik}=D(v_i\cdot v_k)-z_i^T\operatorname{adj}(G)z_k
 =D(q_i\cdot q_k).
\]

第一线类由 (A_{ii}A_{kk}-A_{ik}^2=0) 判定。选代表 (r) 后，用整数Schur residual

\[
 B_{ik}=A_{rr}A_{ik}-A_{ir}A_{rk}.
\]

第二层完成系数精确为

\[
 w_{ik}=\frac{B_{ii}B_{kk}-B_{ik}^2}{D^2A_{rr}}
 =\det(V_{J,r,i,k})^2\in\mathbb Z_{\ge0}.
\]

程序断言全部整除及非负，共线类的传递性与完全分割、五因子的支持分割，以及 (C_f=xyl\in\mathbb Z_{>0})。首个flag的64个实际系数另用独立5阶Fraction行消元核

\[
 \det(V_{J,i,j,k})^2=(A_{ii}/A_{rr})(B_{jj}/B_{ss})w_{sk}.
\]

有限120对的结果为：40对 (D=16)，三维商有7个双元素线类；80对 (D=24)，有4个双元素及6个单元素线类。共有4560个flag、5个归一概率类型；总 (C_f=62,914,560=60W)，每flag的scope严格为全16。

五个类型的 (R_f) 及权重可写成四个不同常数（两个类型的 (R) 相同）：

| (R_f) | 合计 (	heta) |
|---|---:|
| ((7+8\sqrt2+3\sqrt3)/2^{15/4}) | 1/8 |
| (3^{1/4}(1+1/\sqrt2)) | 3/16 |
| (1+\sqrt2) | 5/16 |
| (2\sqrt2) | 3/8 |

全部类型概率、代表及精确质量保留于 `results.json`，没有重复保存全量枚举或重跑021的560循环。

## 5. 实际常数与后续结构瓶颈

预登记的一次固定384-bit Arb批，实际约0.129秒、退出码0，得

\[
 D_{22}=0.83746789731743065419\ldots,
\]
\[
 J_{22}=[1.42053605699085302692149495477727755832\ldots\pm2.39\times10^{-85}],
\]
\[
 I_{22}=[1.30221469258724055505604020332797222323\ldots\times10^{-5}
 \pm1.30\times10^{-90}].
\]

严格比较使用021锁定球的下端，新旧显式 (I) 比约0.9753352015。向上有理上界 (13022147/10^{12}) 由新球上端严格比较取得。未达 (1.18\times10^{-5}) 或 (1.1\times10^{-5})，不从上包围推断真实 (I) 大小。

一子集/零子集的继续flag拆分已从理论判定为重复，**不执行空枚举**：固定一方向后，其rank4商的15个投影无平行类，否则原三方向相关；021已精确核所有560原三方向独立。所以

\[
 \det A=\frac15\sum_iS_iB_i,\qquad
 B_i=\frac14\sum_{j\ne i}S_jB_{ij}.
\]

再接上本批的 (1/3,1/2) 两层，恰为分母120的有序 ((i,j)) flag和。每个无序对只重复两次，完成系数、概率类型、scope相同；sumC由60W变120W，各copy归一权减半，合并后 (J/I) 完全不变。零子集首次拆分原16方向也全singleton，正好回到上述一子集式。它们没有提供新联合信息。

因此下一改进不能只加空的flag层。全权重熵最优属于另登记024；更实质的进一步下降需要保留至少两个非共线商类的联合负矩，而不是全部线性化后丢弃组间相关性。本批没有宣称排除更广泛的矩阵方法或真实目标。

## 6. 复现、边界与白名单

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_I022/three_factor_flag.py
```

Python标准库与python-flint 0.9.0，未调用MATLAB、采样、空间求积或参数优化。严格连续域依据是正分组、AM–GM、独立stable尺度及Hölder；实际工具核有限几何和固定常数。仅报告、登记、脚本、结果、manifest五项归档，无输入快照，未提交推送。022现锁定，后续024不改此目录。
