# upper_cover061：首块零核多基转换

本轮是对不可变旧052树的有限转换，尚未完成全域覆盖，也未改变正式 $C_{cert}/32$。新核链未启动；完成后等待主审。

## 输入与范围

- 父状态：`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover052/batch0043/state.json`，743364738 bytes，SHA256 `065382625da1645e267025ef7c26d292a29a62b65cb7ab97f49dfa9c06995dbf`。
- 输入为52,360个root，passed84,802、split93,920、queue61,478，inflight/resume均0。
- $H_\varepsilon$，$\varepsilon=1/64$，$g_0=1/10000$，$x_1=1$。043/047及056/060的已冻结结论是数学前提，详见登记锁。
- 按 `(path深度, root索引, path)` 祖先优先；最多测试2048个尚未被新父覆盖的split父盒，或资格计算120秒。每16个测试盒保存恢复检查点。
- Python标准库 Fraction 进行全部有理数计算；单进程、单线程，无MATLAB、H、梯度、Green核、log、exp、Gamma新调用。

## 证书与正确性

对16个符号向量 $s\in\{-1,1\}^4$，正混合基给出曲率 $Q_s=I/4-ss^T/36$。逐点取最大合法下界，曲率变为

$$\psi(v)=\|v\|_2^2/4-\|v\|_1^2/36.$$

复用旧失败证书保存的严格下端 $\Delta$ 与有向梯度区间 $[g_j^-,g_j^+]$，计算

$$m(v)=\Delta+\sum_j\min(g_j^-v_j,g_j^+v_j)-4096\psi(v).$$

在每个闭象限内，线性项固定，$Q_s$正定，$m$是凹函数，所以可在每个象限子盒顶点取下界。独立从root及path重建物理盒，并核实anchor位于盒内；对各坐标的 `{lo-anchor,0,hi-anchor}` 去重直积逐一计算，不能仅取16个外角点。每个被测盒均实查81点。还核验缓存H区间及梯度比值包络，原16点最小值和最坏角点与旧证书严格相同。

通过父盒成为新叶，旧passed与queue的严格后代退休；旧split中该父及其所有后代退休。退休清单存不可变源索引及唯一新父；全部历史cache原样保留。未退休节点按原顺序、逐字段保留。新前沿审计包括无重复、antichain、每个root的精确二进制Kraft和1、每个内部节点的子节点/祖先闭合及确定性切分轴。

## 实际运行

PID `57649`，exit0，无子worker。实测资格计算5.611393秒，总28.730657秒；峰值RSS 3084484608 bytes。测试父盒2048，全部通过，总点数165888；最小新log裕量有理下舍入为 `97910373/1000000000`。

首父 `[8816, '']`，末父 `[18868, '']`，本块均root父。退休passed4,306、split2,258、queue0；现passed82,544、split91,662、queue61,478。该块原root子树已经完成，因而缩短证明树但尚未削减待核queue。

源文件运行后SHA仍完全相同；cache全部178,722项保持规范哈希 `f3b174e2bab3fcaed2c291090b378ddac7414393a6fb7a7bceffedca0ffaccac`；新state SHA `d0d62d99d57b3f15d3acb1f4a635b188e196b0de6e155b461460ccbd64a72d31`。完整state和退休账本留在仓库外，精简归档仅含代码、登记、摘要、核验、清单与报告。

## 恢复、审计和剩余工作

原子检查点 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover061/batch0001/checkpoint.json` 已为complete，本块无需恢复。未完成时可对原登记原代码执行 `--resume`；不启动任何新核链。后续有限块需主审后另行登记，不能把此状态交给旧052控制器自动继续。原顶层旧status/fingerprint保留作源元数据，新状态识别 `task061/status061/fingerprint061/conversion061`。

独立062审计进行中；作者运行验证通过尚不替代独立审计。下界正裕量及树转换是有理证书结论，依赖已冻结的正混合抛物面定理；不构成完整全域 $g$ 或 $C_5$ 结论。

复现本块需空输出目录：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B /Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_cover061/convert.py --output-dir /absolute/path/to/new-empty-output
```
