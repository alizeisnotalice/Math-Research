# upper_cover046：直接 Arb 一阶 API 与首个有限覆盖批

## 本批结果及未完成范围

固定 \(\alpha=1\)、\(\epsilon=1/64\)、\(g_0=1/10000\)。本轮新增直接 Arb 一阶 API，
使用017的全部3504锥，计算 \(H_\epsilon=\int_\epsilon^\infty p_tdt\) 及四个物理尾坐标导数。
只在043旧中心作一次新实现资格交叉，再执行一次预登记覆盖批。

实际128个新锚点全部完成，均为256-bit；128个完整闭盒通过043下 log 抛物面的16角点认证。
没有数值升位、细分、缓存命中、H 目标反例或数值诊断停止。按128新锚点上限停止，
状态为 `anchor_checkpoint`，没有执行第二批。

| 项目 | 实际结果 |
|---|---:|
| 初始 level5 根盒 | 52360 |
| 新锚点开始 / 完成 | 128 / 128 |
| 完整密度一阶 kernel 调用 | 129（其中1次为043旧点资格交叉） |
| 认证通过的闭叶盒 | 128 |
| 失败而细分的节点 | 0 |
| 尚未处理的完整 queue 叶盒 | 52232 |
| 几何资格拒绝的中点 | 128 |
| 首个预登记 jitter 成功 | 128 |
| 新 kernel 位数 | 全部256 bits |
| 全批墙钟 / 峰值 RSS | 10.57849775 秒 / 70500352 bytes |

所有通过盒的最差准确有理 log 余量约为 \(0.856770039443402>0\)。
128个中心的 \(H\) 包围下端最小约为 \(0.0006404722869873585\)，上端最大约为
\(0.0008289672729312678\)。这些展示小数只用于阅读；判断全部使用精确端点及分数。
新点完整 kernel 总时间约7.41401秒，单点平均约0.057922秒，范围约0.05679--0.06248秒；
总墙钟还包括资格交叉、完整 queue 的校验/持久化和树检查。

本批按登记的嵌套升序处理前128个根，首根 \((0,0,0,0)\)，末根 \((6,1,0,0)\)，
下一未处理根为 \((6,1,1,0)\)。这是尾坐标较小的区域，不能由本批通过率或耗时推断
接近018中心等困难区域的表现、全局完成时间或所需细分总数。
本轮只称 **partial certified cover + 完整剩余 queue**。没有全域 \(g\) 或新的 \(C_5\) 声明。

## 1. 完整根覆盖与传递条件

令闭 Weyl 区域
\[
 W=\{y\in[0,1]^4:y_1\ge y_2\ge y_3\ge y_4\}.
\]
对每个 \(31\ge i_1\ge i_2\ge i_3\ge i_4\ge0\)，登记完整闭盒
\[
 B_i=\prod_{j=1}^4[i_j/32,(i_j+1)/32]. \tag{1}
\]
此类四元组数为 \(\binom{35}{4}=52360\)。若 \(y\in W\)，取
\(i_j=\min(31,\lfloor32y_j\rfloor)\)，则索引仍非增，且 \(y\in B_i\)。
因此这些闭盒覆盖整个闭 \(W\)，包括相等坐标、网格边界、0和1；没有忽略低维面。
部分完整盒伸出 \(W\)，本方法对整个盒认证，因此这种外放松安全。

16条 body Cauchy 方向在四个尾坐标的置换下仅被重排，所以 \(H_\epsilon(1,y)\) 具有该置换对称性。
如果将来所有根及其替代子叶均认证，则对任意 \(y\in[0,1]^4\)，排序进入 \(W\) 可传递门槛到
整个正单位面。当前只有128个叶盒认证，不能使用这个尚未满足的全覆盖前提。

失败盒沿最长边中点二分，同长时取最小坐标序号；两个子盒闭合，交在中分面，
其并准确等于父盒。根节点用整数编号，路径每一位0/1表示该确定二分的低/高子盒。
树账本保留每个 split 的父、轴、两个子节点；通过叶及完整 FIFO queue 共同构成当前叶集合。
每次持久化均从所有根重建树，检查：父节点确为当前叶、孩子无重名、通过叶和待处理叶不交，
以及二者并等于全部当前叶。特别地
\[
 \#\mathrm{passed}+\#\mathrm{queue}=52360+\#\mathrm{splits}. \tag{2}
\]
本批还实际独立检查：通过根准确为编号0--127，待处理 queue 准确为128--52359，均无子路径。
所以当前整域保持无缺口；未认证部分仍明确待处理。

## 2. 直接 Arb 一阶 kernel

`arb_cover.py` 可安全 import（有 main guard）。API 为
`FirstJetKernel(table, epsilon)`、`eligible(point)` 和
`evaluate(point, bits=256, partial=None, deadline=None, stride=128)`。
调用者应使用已经核对 SHA 的017表。五维输入点用准确有理数；API仅对五维输入的后四坐标求导。
`evaluate` 返回完成标志、值及四个导数球、位数、耗时、相消诊断等；若触发时间检查点，
返回完整前缀的5个和球及下一个锥编号，不把部分锥和当完整密度。

沿用017/018/043的整数缩放：\(z_v=144r_v\)、\(d_C=|\det z_C|\)，
\(D\) 为中心共同分母，\(k_v=Dx\cdot z_v\)、\(l_v=Dz_{v,2:5}\)、
\(e=144D\epsilon=A/B\)，并取
\[
 F=\frac{D^4}{144(2\pi)^5},\qquad
 f_e(k)=\tfrac12\log(e^2+k^2).
\]
完整时间积分为
\[
 H_\epsilon(x)=-F\sum_{3504\ C}d_C\operatorname{DD}_4f_e(k_{C,1},\ldots,k_{C,5}). \tag{3}
\]
本轮新实现不构造有理 log-stencil。初始5元 jet 直接用 Arb 包围
\[
 \left(\tfrac12\log(B^2k^2+A^2),\
       \frac{l_1B^2k}{B^2k^2+A^2},\ldots,
       \frac{l_4B^2k}{B^2k^2+A^2}\right). \tag{4}
\]
共同常数 \(-\log B\) 的四阶差商及空间一阶导数恒为零，故 (4) 可准确省略它。
每锥中心节点准确排序。差商步若分母为 \(d+\delta\cdot u\)，则
\[
 R_0=(P_0-Q_0)/d,\qquad
 R_j=(P_j-Q_j-\delta_jR_0)/d. \tag{5}
\]
整数 \(d,\delta_j\) 都来自准确投影及斜率；所有实数函数、除法、有符号求和及最终 \(F\)
用 Arb 向外运算。整批没有用浮点有限差分或浮点投影判重。
与043不同，此处每锥先带球递推再求和，不能消除数值相消；因此仍保留位数诊断及有限升位规则。

资格交叉只在043原中心作一次256-bit新 kernel，5个球全部包含043已存的值/四梯度球，
且全都通过登记的宽度阈值，耗时约0.0567643秒。因此没有启用登记的384-bit同点重试。
这是一份新实现资格证据，不是新空间点、额外覆盖进展或独立数学审计。

新锚点默认256 bits；只在非有限、区间宽度不足，或中心 \(H\) 相对于门槛尚未分离时，
对同一坐标依次384、512 bits重试。单纯整盒几何下界失败只触发细分，不通过升位伪装成进展。
512 bits仍不可分则保存诊断并保留原盒。若任一完整中心计算已经严格给出 \(H_U<g_0\)，
则保存此固定 \(\epsilon\) 的 H 目标反例并停止，不能说完整 Green 或 \(C_5\) 失败。
本批所有这些重试/停止分支均未触发，不能声称已压力测试。

## 3. 中点退化与有限 jitter

API主动拒绝投影中出现0或重复的中心，以免用普通 DD 处理退化节点。
这只是中心实现的资格限制，\(H_\epsilon\) 本身没有对应的实空间奇点，
完整盒内的零/重复投影不会影响043的正混合定理。

给定盒边界 \(a_j,b_j\)，令中点 \(m_j=(a_j+b_j)/2\)、
\(\eta=\min_j(b_j-a_j)/4096\)。按固定次序只考虑
\[
 m,\quad m+\eta(1,\beta,\beta^2,\beta^3),
 \qquad \beta=1/8,1/16,1/32. \tag{6}
\]
每个扰动坐标小于最短边的 \(1/4096\)，故全部候选严格留在原盒内。
逐个用准确有理投影测试，取第一个具有280个不同非零投影的候选。
若四个均不合格，则诊断停止、保留此盒，不删除盒、不用普通零分母 DD。
不声称三个扰动对所有将来盒都有代数上必成功的保证；账本记录各次有限测试。

本批128个中点全部不合格，全部由 \(\beta=1/8\) 的首个扰动成功。
另外两个候选没有被评估密度或拿来搜索较好函数值；几何拒绝不计作新密度计算。
实际重新检查了128个采用锚点确实是 (6) 首个 jitter，处于原盒，280投影互异且非零。
所有新锚点均位于正单位面。

## 4. 整闭盒的可核验角点证书

只使用043已由主全文核对的全空间下 log 抛物面；不要求完整 \(G\) 在该盒光滑。
对 \(u=(0,v)\)，固定的 det48 基给
\(\|N^{-1}u\|_2^2=v^TQv\)，\(Q=I_4/4-J_4/36\)。
于是对任意实 \(v\)，
\[
 \log H_\epsilon(x+u)\ge\log H_\epsilon(x)
       +\nabla_{2:5}\log H_\epsilon(x)\cdot v-v^TQv/\epsilon^2. \tag{7}
\]
本轮没有把尚待另行归档的043独立审计视为已经完成，也未据此声称全域 \(g\)。

kernel 返回准确有向包围 \([H_L,H_U]\)、各导数包围；Arb 除以正 \(H\) 得
\(g_{j,L}\le\partial_{j+1}\log H\le g_{j,U}\)。再取准确 dyadic
\(\Delta_L\le\log(H_L/g_0)\)。对原盒的16个角点 \(c\)，以实际 jitter 锚点 \(x\) 为原点，
准确计算
\[
 M_c=\Delta_L+\sum_j\min(g_{j,L}(c_j-x_j),g_{j,U}(c_j-x_j))
 -\frac{\sum_j(c_j-x_j)^2/4-(\sum_j(c_j-x_j))^2/36}{\epsilon^2}. \tag{8}
\]
只在全部16个准确有理 \(M_c\) 的最小值严格正时认证。
真实下 log 抛物面对 \(v\) 是凹函数，盒的任一点是角点凸组合；因此角点下界推出完整闭盒
\(H_\epsilon>g_0\)。锚点不必等于盒中点，(8) 已用真实偏移，不使用对称半宽代替它。

持久 checkpoint 对每个通过叶存锚点 cache 键、准确 \(\Delta_L\)、四个梯度有向区间、
最差角点及准确 log 余量；cache 存全部5个球的准确 dyadic 端点及各精度尝试。
因此角点证书可以在不重新计算密度的情况下重放。本批额外实际进行了以下只读审查：
检查完整52360根与叶分割，全部128采用锚点的有理几何资格，384-bit标量 log 验证
\(\Delta_L\) 的方向，准确分数除法验证四个 log 梯度区间包含端点商的范围，
并独立重算全部128×16个 (8)，均与所存最差余量准确相等且严格正。
这次检查没有新密度求值，也没有重跑资格中心或128个新 kernel。

## 5. 持久化、hash 与恢复句柄

源文件、登记、三份上游来源和043四份原件全部在启动时核 hash，组成 fingerprint
`5ac95d23a36f4bfea5fabd7a2d834e9045e774a78d2cc318d7124e46ad82ca9e`。
全部043/017/018来源均只读，未导入生产器。具体冻结来源 hash 见 `registration.json` 和
`results.json`；043原件没有改变。

| 完整对象 | SHA-256 |
|---|---|
| 初始52360根索引 | `69e90a368b9114ac984cdbb15179fabdc45126bf0893675989e12efb1cc46681` |
| 初始完整 FIFO queue | `6643dd72010ef8dedb2c09adfa9b1b5e6b116bcb5eb041cfc48ccbad9fce22d9` |
| 本批完整树（根/splits/passed/queue） | `0b026f4471c143a6c3bb654ec381173b5daaa2aa4c61c1cc027ce9ea21d3ccd5` |
| 剩余完整 queue | `bf926e309f732373af5f6b1232ff6b9be7d1ac56c1ed0ba2a7fe9b5ca60bdf03` |
| 持久 state JSON 字节 | `b23f4ee65bc7ce1d8bdee7f0050b4e492b453241a6594cce8e79f96e7bcfbb2f` |

前四项使用源码中的 UTF-8 canonical JSON（排序key、紧凑分隔、无尾换行）；最后一项是实际
1957168-byte 文件的 hash，含尾换行。持久工作目录为
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover046/`。
可续算句柄是该目录的 `state.json`；另保存同字节只读历史快照 `batch001_state.json`，
以免后续继续合法更新 `state.json` 后失去本批账本。

state 保存完整根、完整 queue、split/pass树、所有已算锚点 cache、资格状态和活动求值前缀，
每个完整叶处理后做树检查并原子写入/fsync。每128锥检查120秒期限；若在 kernel 内到时，
只保存已经完成的整锥前缀（其和用准确 dyadic 区间），恢复时从下一锥接续。
位数重试也保留已失败精度的诊断，不重新执行已有完整精度结果。
本批128锚点上限先触发，没有活动前缀；时间中断、部分 kernel 恢复、细分及数值升位路径
尚未实际触发，必须与已经完成的整队列保存/哈希检查区别。

首批实际命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_cover046/arb_cover.py \
  --upstream-root /Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds \
  --work-dir /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover046
```

本轮停止后没有执行下面的继续命令；它是供主后续有限批调度使用的恢复接口。
下个有限批仍最多128新锚点或120秒，并优先处理现有 queue，不重算已通过叶或旧锚点：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_cover046/arb_cover.py \
  --upstream-root /Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds \
  --work-dir /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover046 \
  --resume --output-dir /tmp/upper_cover046-next-summary
```

新输出目录用于保留本轮冻结摘要；源/登记/上游/043字节变化时 fingerprint 拒绝续用。
迁移机器可提供 `--reference-root` 指向同 hash 的043目录。
若仅复现实验应选择新的仓库外 `--work-dir` 与 `--output-dir`，不得覆写本批历史快照。
本轮只新增046四个必要文件，外置队列/缓存留仓库外；未使用 MATLAB、未提交、未修改他人或冻结来源。
