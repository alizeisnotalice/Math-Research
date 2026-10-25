# Astra vertex017：280真顶点的实际边界 pulling 剖分

本批实际生成了 `P={xi:sum_(v in V)|v.xi|<=1}` 的完整边界剖分：

- **280个原顶点、1882个 facets、3504个边界4单纯形，零人工 ray。**
- 完整面数向量为 `(280,2320,6080,5920,1882)`，连同空面与 `P` 共16484个面。
- 59040项余维一面限制相容性检查全部通过；8760个边界3单纯形各出现两次，有向系数和均为0。
- 准确 `sum |det R|=35509/23887872`，`Vol(P)=35509/2866544640`，与旧26880锥准确一致。
- 两个指定正时间 anchor 的新旧256-bit Arb 区间均相交，差值均包含0。

完整实际顶点、facet contact 与 simplex 表在 `boundary_triangulation.json.gz`，45114字节，SHA256：
`8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9`。
这不是仅证明 pulling 存在；新表已经可直接用于3504项热核计算。本批不运行空间网格、时间积分或全域 `g`，没有新的 `C5` 上界。

## 1. 输入范围与未重复的已锁事实

模型固定 `alpha=1`，`V={(1,+/-1,+/-1,+/-1,+/-1)}`，16个方向等权且不单位化，Lebesgue 测度，Fourier 逆归一化 `(2pi)^-5`。

已锁 `audit/005_vertex_rays_exact.py` 与其结果给出三个真顶点轨道：
\[
 (1,1,0,0,0)/16,\qquad(1,1,1,1,0)/24,\qquad(2,1,1,1,1)/36.
\]
本批只展开符号/排列轨道，数目分别40、80、160；没有重新做旧11类ray的顶点分类。旧 `numerics/body_fan_geometry.py` 的7基锥完整 chamber 覆盖及其3840个符号排列像是另一已锁输入。因为旧锥覆盖整个频率空间，`P` 是所有旧归一 ray 与原点形成的有限单纯形并；所以 `P` 的真正极点必在旧ray表中。005分类因此给出全体280个顶点，而不只是候选子集。原点为内部点，不是极点。

三轨道分母的最小公倍数为144。保存整数 `z=144r`，以其五坐标字典序作为**全局唯一顶点顺序**。每个 `z` 准确满足 `sum_v |v.z|=144`。坐标绝对值最大9。

原始来源及哈希：

- audit clone 的 `005_vertex_rays_exact.py`：`87ffeb3363fe6e25babd67fa01851c0d339bf0399761f45989e8e5f41e95bdbb`。
- 其锁定结果：`ba0f857ea0d590cb1363a17a2281cd4859c7fed2a8321419ae287d22f16662f2`。
- 主 clone 的 `body_fan_geometry.py`：`84dad4a7c74297fcb4e63377276d4a56e7dfa75c038063c2e3f319954aa1ce43`。

本批读取主 `primary/finite_source_transfer.md`，其中 alpha1 的一般 `L1`、全尺度、正部有限源及允许 `S2` 上无穷值的接口已经存在。本批不改写它，也不把几何资格代替后续 `g/I` 证明。

## 2. 65536符号半空间与全部1882个facet

有限恒等式
\[
 \psi(\xi)=\sum_v|v\cdot\xi|
 =\max_{\sigma\in\{-1,1\}^{16}}a_\sigma\cdot\xi,
 \qquad a_\sigma=\sum_v\sigma_vv
\]
给出准确半空间表示 `P=intersection{a_sigma.xi<=1}`。旧源仅保存枚举方法和5817个不同法向量的计数，未找到可复用的完整接触表，所以本批为 **facet任务** 完整枚举65536符号一次，合并为5817个不同整数斜率并保留其多重数，总多重数65536。这里没有重复旧vertex分类。

对每个斜率 `a`，接触顶点条件是整数等式
\[
 a\cdot z=144.
\]
全部280个顶点同时满足所有半空间。任何非空接触集位于不经过原点的超平面 `a.z=144`，所以其仿射维数是接触向量的线性秩减1。故接触矩阵线性秩5当且仅当它是4维facet。

这也证明facet表完备：一个有限全维多面体的任一facet，在其相对内部有至少一个有效定义不等式取等；该不等式的取等面必须含此facet，因两者均是proper face而facet极大，故恰等于该facet。所有定义不等式都已枚举，没有漏掉其它方向的支持面。

实际接触秩分布为：秩0有2055个斜率，秩1有1160个，秩2有720个，秩5有1882个。facet按顶点数分布为：

| facet顶点数 | facet数 |
|---:|---:|
| 5 | 1120 |
| 6 | 640 |
| 7 | 80 |
| 8 | 10 |
| 10 | 32 |

每个facet保存准确支持斜率和全部接触顶点ID。不同facet接触集互不相同；秩5也保证其支持斜率唯一。

NumPy仅进行整数符号与矩阵乘法：斜率每坐标绝对值不超过16，`z`每坐标绝对值不超过9，所以接触点积的先验绝对值至多720，远低于 `int64` 界。没有浮点rank、容差、SciPy凸包候选或MATLAB几何调用。所有rank由 Python任意精度整数交叉消元和gcd规范化得到。

## 3. 完整面格的有限构造

面用280位顶点掩码表示，从完整顶点集出发，与全部facet接触掩码反复求交，直到队列完全处理、没有新面。空面同时保存。

掩码相交等价于几何面相交：若凸组合在支持不等式取等，则其所有正系数顶点也必须取等。因此两个面相交正是共有顶点的凸包。每个多面体面又是包含它的facets之交，故这个有限闭包恰好得到完整面格。

对所有非空proper face，维数由整数顶点向量的rank减1计算；空面维数-1，整个 `P` 维数5。每个面的余维一子面通过与全部facets相交并筛选维数 `d-1` 得出。这不是从待验证的simplex邻接反推面格。

结果为
\[
 f(P)=(280,2320,6080,5920,1882),\qquad
 280-2320+6080-5920+1882=2.
\]
各父面维数0到5的余维一incidence数依次为
`280,4640,18240,24320,11840,1882`。Euler等式只是附加一致性检查；完备性来自完整半空间及交集闭包。

## 4. 只用原顶点的确定性 pulling 与覆盖证明

对每个非空proper face `F` 按维数从低到高计算 `T(F)`。零维面只有该顶点。若 `dim F=d>=1`，令 `v` 是其全局顺序最小顶点，定义
\[
 T(F)=\{\operatorname{conv}(v,\tau):
 G\text{ 是 }F\text{ 的facet},\ v\notin G,\ \tau\in T(G)\}.
\]
程序没有增加重心、原点或任何其它边界顶点。每个生成simplex检查顶点个数 `d+1`、无重复、全部属于该面以及线性rank `d+1`；在proper face上这正是仿射独立。

覆盖与不重叠可逐维证明。对 `relint F` 中任一点 `x!=v`，从 `v` 沿 `x-v` 向前的射线在紧凸集 `F` 有唯一末端 `y`。任何含 `v` 的facet在该射线上离开 `v` 后严格满足其支持不等式，故末端落在一个不含 `v` 的facet中；于是 `x` 属于上述某个锥化simplex。相对内部点落在两个不同锥化块的内部会使同一射线有两个不同末端，矛盾；同一facet内的simplex不重叠由归纳得到。有限闭并的极限覆盖其边界。这给每个 `F` 的完整三角剖分。

为检查不同facet上的实际剖分相容，程序对 **全部59040项proper-face余维一incidence** 检查：

1. `T(F)` 各simplex与子面 `G` 的顶点交集，其所有最高维simplex恰为 `T(G)`。
2. 每个较低维交集也包含于某个 `T(G)` 的最高维simplex中。

因为 `G` 是支持面，simplex与 `G` 的几何交集恰是顶点交集的凸包。以上两个准确集合检查证明相应simplicial complex限制完全相同。沿面格的饱和链迭代，所有面之间的限制均相容。因此相邻原facet的剖分拼接成整个边界的单纯复形。

最终1882个facet产生3504个4单纯形。每facet的simplex数分布为
`1:1120, 2:448, 3:128, 4:112, 6:32, 8:10, 12:32`。全部280个真顶点都被使用。另直接逐项检查每个边界3单纯形出现两次，且由准确行列式定向后其系数和为0；共有8760个这样的3单纯形。

因为原点在 `P` 内部，从原点对这些边界4单纯形锥化，其锥覆盖 `R5`、相对内部互不交，交界为零测有限低维集合。无需对所有锥作两两数值相交检测。

## 5. 准确体积、实际表格式与新热核公式

对每个表中4单纯形，五个整数顶点为 `z1,...,z5`。Bareiss整数消元计算 `D=det(z1,...,z5)`，每次整除都断言准确，最终 `D!=0`。归一ray是 `rj=zj/144`，故
\[
 |\det R|=|D|/144^5.
\]
结果准确满足
\[
 \sum_{3504}|\det R|=35509/23887872
 =3840\sum_{7\text{ old base cones}}
 \frac{|\det N|}{\prod_j\psi(n_j)}.
\]
除以 `5!=120` 得 `Vol(P)=35509/2866544640`。覆盖已经由上一节证明，体积等式是跨表示的独立核对，不是唯一覆盖依据。

压缩JSON表的schema为：

- `denominator=144`。
- `vertices[280][5]`：字典序整数 `z=144r`。
- `facets[1882]`：整数 `normal`、完整接触 `vertices` ID、所属 `simplices` 的全局索引。
- `simplices[3504]`：五个递增顶点ID和 `det_numerator=abs(D)`。
- 源哈希、完整面格输入哈希、全局pulling顺序说明。

每个simplex属于同一facet `a.r=1`，所以其整个正锥上 `psi(sum lambda_j rj)=sum lambda_j`。由准确覆盖和绝对可积的锥参数积分，对任意实 `x`、正 `t`，
\[
 \boxed{p_t(x)=(2\pi)^{-5}\sum_{3504}
 \frac{|\det R|}{\prod_{j=1}^{5}(t-i x\cdot r_j)}}.
\]
以表中整数计算的完全等价形式是
\[
 p_t(x)=(2\pi)^{-5}\sum_{3504}
 \frac{|D|}{\prod_{j=1}^{5}(144t-i x\cdot z_j)}.
\]
没有遗漏144的幂或 `5!`。`t>0` 保证各因子实部正，投影重复或为0都不会产生除零。此公式无需generic点、partial fraction留数或人工pole相消。

## 6. 两个且仅两个资格anchor

完整几何生成后，以256-bit Arb对每个新表达式和锁定旧26880项表达式各算一次：

| x | t | 新旧共同值的显示近似 |
|---|---:|---:|
| `(1,7/8,3/4,5/8,1/2)` | `1/3` | `2.42607069231531572355e-5` |
| `(1,1,1,1,0)` | `1/8` | `0.000384577364219553057791` |

两组准确球都严格为正；新区间的虚部包含0，新旧实部区间相交，差值球分别为 `[+/-5.79e-78]`、`[+/-9.18e-77]`。后一个空间点虽属于Cauchy `S2`、全时间Green为无穷，但这里 `t=1/8>0`，热核值有限，资格交叉没有混淆两者。

新表达式各耗时约 `.00548`、`.00478` 秒，旧表达式约 `.667`、`.616` 秒；新项绝对值和除以总值分别约1.3866、7.0349，仅描述这两个固定点的相消。未在其它点、时间或区间箱进行测试；未重放旧Fejer或失败全箱。单点吻合不承担几何证明，更不承担全域 `g` 证明。

## 7. 进程、检查点、复现和剩余工作

实际四个计算进程均退出0，没有挂起进程：

| 阶段 | PID | 运行秒数 |
|---|---:|---:|
| contacts | 39657 | .1244 |
| full faces | 39805 | 2.9360 |
| pulling/table | 39912 | .6190 |
| two anchors | 40009 | 1.4193 |

各阶段均在120秒checkpoint内完成。恢复状态分别位于
`/tmp/astra_vertex017/contacts_checkpoint.json.gz`、`faces_checkpoint.json.gz`、`pulling_checkpoint.json.gz`；输入/源码指纹不匹配时拒绝续算。完整接触及面格缓存留该目录，未放进repo。anchors完成记录为同目录 `anchors_checkpoint.json`；资格交叉不自动重复。

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_vertex017/contacts.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_vertex017/faces.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_vertex017/pulling.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_vertex017/anchors.py
```

采用有限维凸包与极点技能的假设/证书核验方式；本证明只涉及显式五维多面体，不调用无限维Krein–Milman结论。没有外部浮点凸包软件或采样作为几何证据。归档白名单为登记、整数辅助与三个生成器、三个结果摘要、实际压缩表、anchor代码和结果、本报告及manifest。

未提交或推送；没有修改旧fan、原分类或其它代理目录。本轮完成的是原来缺失的真实顶点边界表和可审计热核接口。如何利用它取得完整时间/空间下界属于下一批；当前 `C5<=32` 不变。
