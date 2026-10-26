# upper_cover066：面向待证前沿的零核多基转换

本轮完成所有仍含queue后代、且尚未测试的源split父盒的有限转换。queue由61,478降至13,764，减少47,714（77.611503%）。没有重算H或梯度，没有启动新核链；完整全域覆盖尚未闭合，正式 $C_{cert}/32$ 未变。

## 输入、登记和选择

唯一完整输入为 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover061/batch0001/state.json`，SHA256 `d0d62d99d57b3f15d3acb1f4a635b188e196b0de6e155b461460ccbd64a72d31`，735,202,818 bytes。该061状态含passed82,544、split91,662、queue61,478、cache178,722，52360个roots；inflight/resume均0。061及其锁不变。056/060的正混合数学结论与062已经通过的树接入审计是冻结前提。本轮不将已有Delta下端另算新log。

保留模型 $\varepsilon=1/64$、$g_0=1/10000$、$x_1=1$。每块事先写登记：父检查点SHA、父逻辑frontier SHA、确切候选前缀、代码与完整源SHA。候选仅含仍活且未测、具有queue严格后代的源split，按 `(path深度, -queue后代数, root索引, path)` 确定排序；取前2048候选，测试至候选耗尽或资格计算120秒。通过祖先后跳过其全部后代，失败仅记局部仍未证，并继续其子树，不能视为H反例。

`initial_queue_descendants`是源计数。因为所有严格祖先先测试，未测的活候选不可能已有接受的严格后代；因此选取时其活queue后代数仍等于源计数。每个源父最多测试一次。16次测试保存一次可恢复前缀；每块完整校验成功后才登记下一块。

## 数学证书和树不变量

复用旧 `failed_certificate` 的Delta、gradient与cache，独立按root/path重建物理box，并验证anchor在盒内。16符号正混合基的逐点最大下界对应

$$\psi(v)=\|v\|_2^2/4-\|v\|_1^2/36,$$

$$m(v)=\Delta+\sum_j\min(g_j^-v_j,g_j^+v_j)-4096\psi(v).$$

每个闭象限内 $m$ 凹，故分别取子盒顶点；逐坐标集合为 `{lo-anchor,0,hi-anchor}`，本轮每个被测父盒均为81点。仅当全部点严格正才转为新叶。严格有理数Fraction运算同时核验旧16点最小裕量及最坏角点完全重现、H/梯度区间比值包络、曲率边界。继承Delta的合法性仍以冻结源证书为前提；主独立Delta复核不由本作者验证替代。

新父通过后，删除其旧passed与queue的严格后代，删除其自身及后代旧split，加入该新叶。每个退休记录给出不可变061源数组索引、node与唯一新父。每块均核验：活前沿antichain、每root精确二进制Kraft和1、源split子节点和最长边确定性分裂轴、全树祖先闭合。每块前沿由源数组索引与新证书构造，未退休记录逐字段保持原对象。全部cache逐块规范SHA不变。

## 实际运行和结果

原运行PID61162，17块，exit0，OS确认已不存在；总178.380484秒，资格计算累计98.480351秒，峰值RSS3106684928 bytes。单进程、单线程，MATLAB/H/梯度/Green核/log/exp/Gamma新调用均0。

| 项目 | 数量 |
|---|---:|
| 测试父盒 | 34038 |
| 新通过父盒 | 15143 |
| 测试仍未证父盒 | 18895 |
| 精确检测点 | 2757078 |
| 退休旧passed | 3043 |
| 退休旧split | 35614 |
| 退休queue | 47714 |
| 最终passed / split / queue | 94644 / 56048 / 13764 |
| 剩余尚未测试的eligible源父 | 0 |

测试父盒深度0/1/2分别10,012/14,282/9,744；接受父盒深度0/1/2分别2,871/5,418/6,854。所有新接受margin的有理下舍入最小值为 `1001/1000000000`，严格正。所有34,038个测试源node无重复。

完整新state为 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover066/final/state.json`，691,008,601 bytes，SHA256 `42c67156078cb0eb134910a5248821eb3c826c96981afb3d1d5b2fe3242a871b`。源061运行后SHA完全不变；cache规范SHA始终 `f3b174e2bab3fcaed2c291090b378ddac7414393a6fb7a7bceffedca0ffaccac`。最后检查点SHA `791f4a18c5e6506b724a4580ec52580cae9ec3b9712349e4b522bd06e1f463db`。

## 稀疏保存和重建

每块只保存登记、完整测试证书、退休delta、摘要与链式检查点；未逐块复制cache或完整state。最后仅物化一个完整state。`prefix_checkpoint.json`可在真实故障后继续本块；`checkpoint.json`是已完成不可变块边界；后块登记包含前块文件SHA。

完整delta链和完整输入保存在仓库外，代码、登记、统计、重建receipt与本报告精简归档在本目录。源state顶层旧元数据原样保留；仅passed/queue/splits改变，新增task066/status066/fingerprint066/conversion066。后续程序须识别新状态字段，不能交给旧052控制器自动重启。

确定性重建命令（仅从061+已完成delta链恢复，不再测试父盒）：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/homebrew/opt/python@3.15/bin/python3.15 -B /Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_cover066/convert.py --work-dir /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover066 --reconstruct
```

已实际执行：重建PID61975，exit0且OS确认不存在，用时21.505964秒，从061重新载入全部cache并重放17块，产生完全相同完整state SHA。`--resume`支持未完成前缀恢复；所有源/code/checkpoint锁必须匹配。`STOP`文件可在新派发前停止并保留检查点。

## 结论边界

本轮全部可测源父已处理；仍有13,764个queue叶，不能仅凭旧anchor的该证书闭合。下一步方法和任何新核求值需主另审。独立062正在扩审首块和最终接入；作者核验与确定性重建不替代该独立审计。没有提交或推送，没有修改061及其他人的文件。
