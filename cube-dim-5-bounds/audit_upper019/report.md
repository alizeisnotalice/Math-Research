# upper019：完整形状 I 上界及轨道 Minkowski 改进独立审计

**通过。** 第011轮的完整面消元及负矩 Gamma 上界成立；主新增的轨道平均后 Minkowski 路线也成立，并给出更强的完整证书

\[
 I=\int_{[0,1]^4}G(1,y)^{10/7}dy
 \le I_{\rm orb}=1.57928752881910986698\ldots\times10^{-5}
 <\frac{157929}{10^{10}}<\frac1{60000}.
\]

这里认证的是完整全时间、全正面 (I) 的上界，不是中段验收门槛或采样估计；不需要再加短/长时间或角点管误差。这个上包围仍高于目标 (11/10^6)，不能据此判断真实 (I) 超过目标。已知全域 (g\ge1/65000) 仍不足以使传递成本降到32，正式 (C_5\le32) 不变。

本批只写独立 clone `Math-Research-cube5-upper-audit-cont/cube-dim-5-bounds/audit_upper019/`；主目录及作者目录只读。未提交、未推送，未占用 MATLAB。没有重复4368个五方向子集或任何采样、空间求积、参数搜索。

## 1. 锁定对象与适用范围

固定16条未单位化体对角方向 (V=\{(1,\pm1,\pm1,\pm1,\pm1)\})，独立标准对称 (\alpha=3/2) stable 时间，权均为1；Lebesgue 测度与逆 Fourier 因子 ((2\pi)^{-d})。传递对象仍为中心共同半径轴平行方体、一般 (L^1)、全部正尺度。

关键源版本：

- `upper/round011_shape_integral.md`：`eda3dd844c83c8a333505f0023923fe5fccd1fc24edec9f81b98016dc6b6b96b`。
- Astra `astra_upper001/report.md`：`99d1ad30fa9b1df3f59bc88ae0dfa782c4928280efb3550ecfd11a0cee58a462`。
- 主 `round006/orbit_minkowski_proof.md`：`ba298fabe815adceac6cb9a0001268d04177862134722042cb32326bde41c3f0`。
- 主轨道源码：`6023fa367e0f4da8c535fa059efbd822d8910c1bbfd0e1bc16363de9eee2495e`。

全部13项输入、原作者根目录、项目相对路径及SHA见 `registration.json`。Astra本批只消费其第1、3节的正Gaussian/时间积分、(T^{-\beta}) 概率倾斜与形状表示；**没有审核其 trace 排除、全部多项式矩、minorant 完整性或低阶证书**。这限制了本报告对上游的验收范围。

## 2. 正表示与 R 的归一化

令 (\beta=3/4)，(S_v) 独立正stable，(mathbb E e^{-\lambda S_v}=e^{-\lambda^\beta})。写

\[
 T=\sum_vS_v,\quad P_v=S_v/T,\quad A=\sum_vS_vvv^T=TM,
 \quad M=\sum_vP_vvv^T.
\]

因为 (S_v>0) a.s.、方向张成五维，(M) a.s. 正定；(operatorname{diag}M=\mathbf1)。条件Gaussian covariance是 (2t^{4/3}A)，密度因子为 ((4\pi)^{-5/2})。先对全部正时间积分给

\[
 G(x)=\kappa\mathbb E[(\det A)^{-1/2}(x^TA^{-1}x)^{-7/4}],
 \qquad\kappa=\frac{3\Gamma(7/4)}{2^{7/2}\pi^{5/2}}.
\]

本批 SymPy 实际独立积分了 (t^{-10/3}e^{-z/(4t^{4/3})})，再乘Gaussian因子，确认此 (kappa)；(4) 没有遗漏。非负Tonelli允许上述量取无穷，方向线上Green的无穷值不破坏该表示。

矩阵函数在 (A=TM) 下恰带 (T^{-5/2+7/4}=T^{-3/4})。Gamma负矩公式给

\[
 Z:=\mathbb ET^{-\beta}=\frac1{16\Gamma(7/4)}>0.
\]

用 (T^{-\beta}/Z) 倾斜原概率，再推送至 (P) 得真正概率 (R)。于是

\[
 G(x)=c\mathbb E_RH(M,x),\quad
 c=\kappa Z=\frac3{2^{15/2}\pi^{5/2}},\quad
 H=(\det M)^{-1/2}(x^TM^{-1}x)^{-7/4}.
\]

没有假定 (T) 与 (P) 独立。原 (S_v) 的置换不改变联合律与倾斜量 (T)，因此 (R) 继承全部五坐标带符号置换对称。改变方向整体符号不改变 (vv^T)，所以使用半套16方向并不丢失交换第一坐标的对称。

## 3. 球面、cube-face 与关键 1/5 因子

置 (q=10/7)、(d=\det M)。固定正定 (M) 时

\[
 H(M,x)^q=d^{-5/7}(x^TM^{-1}x)^{-5/2}
\]

是空间的负5次齐次函数。五维Gaussian积分的径向部分是 (Gamma(5/2)/2)，故

\[
 \int_{S^4}(\theta^TM^{-1}\theta)^{-5/2}d\omega
 =\frac{2\pi^{5/2}}{\Gamma(5/2)}d^{1/2}
 =\frac{8\pi^2}{3}d^{1/2}.
\]

每个cube面上 (x\cdot n=1)，径向投影的Jacobain为 (d\omega=|x|^{-5}d\sigma)。它与函数的 (|x|^5) 缩放恰抵消，所以

\[
 \sum_{F\text{为10完整面}}\int_F H(M,x)^q,d\sigma
 =\frac{8\pi^2}{3}d^{-3/14}.
\]

完整五坐标置换及符号对称使取 (R) 期望后的160块面正交小片积分相同。因此正unit face的精确因子是

\[
 \frac{8\pi^2/3}{10\cdot16}=\frac{\pi^2}{60}.
\]

单独扩大至整平面的 (pi^2/12) 比它大5倍；这个因子不能只从面内反射获得。边界交叠为面测度零，所有非负积分均由Tonelli处理。确定性校准 (M=I) 时 (H^q=|x|^{-5})，每个正unit face积分正是 (pi^2/60)，角因子与等号一致。

第011轮只在形状期望上作Jensen，得到正确的完整上界

\[
 I\le c^q\frac{\pi^2}{60}\mathbb E_Rd^{-3/14}.
\]

后述负矩有限同时补足其可积性；不必预先排除奇异角点或假定面迹连续。

## 4. 均衡 Cauchy–Binet 权重与反向关联

写 (V) 为 (5\times16) 方向矩阵。Cauchy–Binet给

\[
 d=\sum_{|B|=5}w_B\prod_{i\in B}P_i,\qquad w_B=\det(V_B)^2\ge0.
\]

(VV^T=16I)，故权重总和 (W=16^5)。除方向传递性证明外，本批给出独立的导数核对：对 (det(\sum_i z_iv_iv_i^T)) 在 (z_i=1) 求导，

\[
 \frac{\sum_{B\ni i}w_B}{W}
 =v_i^T(VV^T)^{-1}v_i=\frac5{16}.
\]

因此加权AM–GM严格给出方向正确的全simplex内部下界

\[
 d\ge16^5\prod_iP_i^{5/16}.
\]

零权重项可删去，退化边界以连续性延拓。这不声称 (d) 有统一正常数下界。代码只检验tight frame、16个rank-one行列式身份和4个预登记确定性形状，未重扫全部基。

四个形状为均匀 (P_i=1/16)、一个 (241/256) 加15个 (1/256)、权重按整数1至16归一化，以及单方向原子边界。用精确比较 (d^{16}\ge16^{80}\prod_iP_i^5) 全部通过；均匀形状取等，原子边界 (d=0)。单重方向的精确式也核为 (d=(16\eta)^4(5-64\eta))。这些是有限退化检错，严格全域性来自上述AM–GM证明。

旧指数 (b=3/14)、(a=5b/16=15/224)、(gamma=16a-\beta=9/28\in(0,\beta)) 给

\[
 \mathbb E_Rd^{-b}\le\frac{16^{-5b}}Z
 \mathbb E[T^\gamma\prod_iS_i^{-a}].
\]

原独立概率下 (T^\gamma) 对每个坐标递增，而 (prod_iS_i^{-a}) 递减，故反向关联的方向确为上界。单坐标两副本公式的积差非正；逐坐标条件化后，两条件期望仍保持相应单调性，归纳给产品概率的不等式。对两非负函数分别取 (min(N,\cdot))，先用有界结果，再单调收敛，得到

\[
 \mathbb E[T^\gamma\prod_iS_i^{-a}]
 \le\mathbb ET^\gamma(\mathbb ES^{-a})^{16}.
\]

有限性无需循环假设：所有负矩存在；正矩的Laplace积分在零端须 (gamma<\beta)，在无穷端须 (gamma>0)，均满足。此处仍未假定 (T/P) 独立。由Mellin公式得到旧上包围

\[
 Q_*=16^{5/14}\frac{\Gamma(7/4)\Gamma(4/7)}{\Gamma(19/28)}
 \left[\frac{\Gamma(61/56)}{\Gamma(239/224)}\right]^{16},
\]

从而 (I\le1.9382520037464933788\ldots\times10^{-5}<1/50000)。这条证书有效，随后被下节的更强完整证书承接。

## 5. 轨道平均、稳定子24及 Minkowski

令 (Gamma) 是全部 (2^5 5!=3840) 个带符号置换，

\[
 \bar H_M(x)=\frac1{3840}\sum_{U\in\Gamma}H(UMU^T,x).
\]

正面小片 (F_+=\{x_1=1,0\le x_j\le1\}) 的轨道为160块。保持该小片的变换只能置换其余四坐标，稳定子恰为 (S_4)，大小24；独立改号或移动固定坐标不能保持小片内部。这是**面小片的稳定子**，不是某一特殊矩阵 (M) 的共轭稳定子；后者的大小不能代入此计数。全群中有重复矩阵共轭项也不改变平均身份。

有限Jensen、表面测度保持的正交换元以及 (3840=24\cdot160) 给

\[
 \int_{F_+}\bar H_M^q
 \le\frac1{160}\sum_F\int_F H(M,x)^q
 =\frac{\pi^2}{60}d^{-3/14},
\]
\[
 \|\bar H_M\|_q\le(\pi^2/60)^{7/10}d^{-3/20}.
\]

R的共轭不变性给 (G/c=\mathbb E_R\bar H_M)。(M\mapsto\bar H_M) 是可测正函数族，且下述矩界使其 (L^q) 范数在 (R) 下可积；因此积分Minkowski合法，给

\[
 \|G\|_q\le c(\pi^2/60)^{7/10}\mathbb E_Rd^{-3/20}.
\]

也可先对有限形状子域使用该不等式，再非负取极限；无需假定所有退化形状有统一条件数。(M=I) 的确定性校准同时使固定形状Jensen与轨道范数估计取等。

新的 (b=3/20)、(a=5b/16=3/64) 恰使 (16a-\beta=0)。故

\[
 \mathbb E_Rd^{-3/20}\le\frac{16^{-3/4}}Z
 \mathbb E\prod_iS_i^{-3/64}
 =2\Gamma(7/4)\left[\frac{\Gamma(17/16)}{\Gamma(67/64)}\right]^{16}
 =:J_*.
\]

这一步只用原 (S_i) 的独立性，**不需要旧反向关联放宽**。所有负矩有限，因而Minkowski的范数可积条件闭合。于是

\[
 I\le c^q\frac{\pi^2}{60}J_*^q=I_{\rm orb}.
\]

## 6. 不同代数路径的实际独立核验

仅执行预登记的一组384-bit Arb公式及有限精确测试。Gamma负矩使用未移位式

\[
 \mathbb ES^{-a}=\frac{\Gamma(a/\beta)}{\beta\Gamma(a)},\quad
 Z=\frac{\Gamma(1)}{16\beta\Gamma(\beta)},
\]

与生产脚本的 (Gamma(1+a/\beta)/\Gamma(1+a)) 路径不同。最终 (I) 用七次幂闭式计算，再由Fraction整数幂比较包围七次根：令

\[
 C_7=\frac{3^{10}}{2^{65}60^7\pi^{11}},\quad
 m_a=\frac{\Gamma(a/\beta)}{\beta\Gamma(a)},\quad
 p_\gamma=\frac{\Gamma(1-\gamma/\beta)}{\Gamma(1-\gamma)}.
\]
\[
 I_{\rm old}^7=C_7\Gamma(7/4)^7p_\gamma^7m_{15/224}^{112},\qquad
 I_{\rm orb}^7=C_7\Gamma(7/4)^{10}m_{3/64}^{160}.
\]

90位十进制有理根端点均由整数幂方向检验。独立的 (Z,c,Q_*,I_{\rm old},J_*,I_{\rm orb})、新旧比例及所需 (g) 共8个区间全部包含于作者/主原保存打印球；不是只检查两个浮点中心接近。

\[
 I_{\rm orb}\in[
 0.00001579287528819109866987013756481830367768,
 0.00001579287528819109866987013756481830367769].
\]

这实际证明父任务要求的 (I_{\rm orb}<157929/10^{10})，还给出略强的可用有理界 (I<1579288/10^{11})。新旧显式上包围比例约0.8147998948；两者是同一Green核、同一面测度，比较合法。

作为旁核，精确Fraction计算 (eta^5 5!/(\beta)_5=648/1463)，与作者的第一行列式矩系数一致；这不等于对Astra全部矩接口的审计。

实际命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/audit_upper019/check_interfaces.py
```

退出码0；python-flint 0.9.0、SymPy 1.14.0、Arb384位，具体运行记录及完整有理区间见 `check_results.json`。未运行作者程序、未重扫全部基、未调用MATLAB。纯解析审核覆盖连续域；有限测试只是归一化、代数和退化检错。

## 7. 一般传递与剩余缺口

本批完整 (I) 与已审 `alternative_alpha/transfer.md` 使用完全相同的 (G)、模型、Lebesgue测度、unit face和 (q)，所以可以直接输入一般 (L^1)/全尺度传递

\[
 C_5\le(10/3)^{10/7}I/g^{10/7}.
\]

不把新 (I) 再按截断门槛处理，也不额外加入第010轮的短时或角点管误差。当前已证 (g\ge1/65000) 可给一个有效但弱于32的传递成本；即使暂以旧条件 (g=1/10000) 配这个显式 (I_{\rm orb})，成本仍超过32。独立核算该公式的成本32门槛为

\[
 g_*=(10/3)(I_{\rm orb}/32)^{7/10}
 \in[0.00012829112235487426982580271893663047437939,
       0.00012829112235487426982580271893663047437940].
\]

这是这一个显式上包围的充分 (g) 门槛，不是对真实 (inf G) 的结论。原目标 (I\le11/10^6) 仍未证；上包围大于目标也不排除通过更紧负矩或减少形状损失达到目标。正式上界仍32。

白名单仅为本报告、`registration.json`、`check_interfaces.py`、`check_results.json`、`manifest.json` 五项；不存重复输入快照。检查脚本按项目相对路径或锁定作者路径读原件，强制SHA一致。导入后可用 `check_interfaces.py --hash-only --source-root <主仓库>` 只验来源而不重跑数学公式。没有独立不可替代的快照字节；本批不扩到Astra002等后续接口。
