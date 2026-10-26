# upper_framecover077：075全frame曲率的有限零核接入

本轮将已严格审过的075全frame罚项 $Q_{full}=I_4/16$ 接入066冻结树，仅复用178,722个已有H/gradient缓存及保存Delta。两个独立有限进程自然完成全部尚未测试的eligible源父，queue从13,764降至684；没有新核调用，未闭合全域g，也未改变正式 $C_{cert}/32$。

## 冻结输入与数学对象

唯一完整输入 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover066/final/state.json`，SHA256 `42c67156078cb0eb134910a5248821eb3c826c96981afb3d1d5b2fe3242a871b`，691,008,601 bytes，原passed94,644/split56,048/queue13,764，roots52,360，inflight/resume均0。所有旧状态与锁只读。

075完整独立报告SHA `3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225`，manifest SHA `1d5e270eb1664d714ec80f0b371420d40028b24c3943196b9ab41a0b3de1f9d1`；启动时实际核对并读取完整报告/manifest，manifest内4个文件也全部核锁。069 final对066全17块及最终接入独审的报告SHA `bda9cc665e71c5004518b5881f5c625010a6d988a68a760c8666af38d183675f` 同步核锁。

模型仍为 $\varepsilon=1/64$、$g_0=1/10000$、$x_1=1$。075的全frame正混合证明给出

$$\log H_\varepsilon(x+u)\ge\log H_\varepsilon(x)+\nabla\log H_\varepsilon(x)\cdot u-\frac{|u|^2}{16\varepsilon^2}.$$

故对尾位移v，使用旧证书的保存下端Delta及有向gradient，安全下界为

$$\Phi(v)=\Delta_L+\sum_{j=1}^4\min(g_{j,L}v_j,g_{j,U}v_j)-256\sum_{j=1}^4v_j^2.$$

两个仿射函数的最小值全局凹，减正二次仍凹，闭盒16个角点严格正即可；本数学对象不要求旧056的81象限点。anchor理论上无需在盒内，本轮仍从源root/path重建盒并实际核其在内，核缓存H区间及gradient比值包络，同时重现原失败16角点值作为来源完整性检查。新full-frame16角点逐一用Fraction有理数计算；无新增log。

## 有限进程、候选与稀疏状态

每个计算程序总wall上限600秒（含hash、解析、计算、核验和落盘），540秒进入正常检查点，590秒alarm兜底；每进程只处理一个至多4096候选/测试父盒的有限块，完成后退出。后块由代理以另一个新进程登记启动，无无限controller循环。

候选只含仍有queue严格后代的未测活源split，按 `(path深度, -当前queue后代数, root索引, path)` 确定排序。每块先登记确切候选前缀、唯一066源SHA、数学/代码锁、父checkpoint SHA及逻辑frontier SHA。通过祖先即跳过全部后代；失败仅记仍未证，继续其子树，绝不将曲率证书失败当作H反例。每16测试保存可恢复前缀。

块之间只保存完整资格证书、退休delta、摘要与链式checkpoint；不复制未变cache。每个退休passed/queue严格后代及父/后代split均给源数组索引和唯一新父。每块检查完整52,360 roots的精确二进制Kraft和1、antichain、分裂子节点与祖先闭合、最长边确定性轴；未退休记录逐字段复用源对象，cache规范SHA逐块恒定。完整最终状态只物化一次；另做确定性重建验证。

## 实际结果

| 项目 | 第1进程 | 第2进程 | 合计/最终 |
|---|---:|---:|---:|
| PID | 67286 | 67350 | 均exit0，OS查无 |
| 登记候选 | 4096 | 1109 | 5205 |
| 实际测试父 | 3081 | 951 | 4032 |
| 新通过父 | 2594 | 571 | 3165 |
| 仍未证父 | 487 | 380 | 867 |
| 新fullQ角点 | 49296 | 15216 | 64512 |
| 进程总wall收据/秒 | 15.947203 | 19.239527 | 均远低600 |
| queue | 1984 | 684 | 684 |

资格计算累计7.440019秒；接受margin的最小有理下舍入为 `10092297/500000000000`，严格正。测试源父无重复，深度0/1/2分别2,982/962/88，接受深度分别2,501/586/78。退休源passed850、source split10,765、queue13,080；最终passed96,959、split45,283、queue684。无剩余未测试eligible源父，过程自然停止。

全部178,722 cache项规范SHA始终 `f3b174e2bab3fcaed2c291090b378ddac7414393a6fb7a7bceffedca0ffaccac`；唯一066源SHA再次核对完全不变。最终完整state `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_framecover077/final/state.json`，674,754,263 bytes，SHA `c55918ec39ef62629d53d20e1aaed7a5faf454d52212c5897bd16fd15056ea9b`。其唯一旧顶层变换为passed/queue/splits，新增task077/status077/fingerprint077/conversion077，其他066顶层值全部保持。

## 剩余684叶与证据边界

剩余684 queue来自197个roots、342个最近失败父，深度2有664叶、深度3有20叶。这些父缓存H_lower全部严格大于目标1/10000；未证是闭盒margin不足，不能据此称H反例。queue叶本身没有新增H求值。

具体剩余清单在 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_framecover077/remaining_queue.json`，按父缓存H_lower、父fullQ margin、root/path确定困难优先。每行保存叶node/root/独立交错二进制重建的物理box、最近失败父源split索引、parent anchor、缓存H上下端、Delta及fullQ失败margin；`remaining_summary.json`给出精确范围。附录整理程序PID67673，1.319752秒，无核调用。

未授权启动070或其它新H链，本轮不启动。剩余684需主另审下一步；本轮不形成全域g候选或新的C5结论。

## 可恢复性和确定性重建

run registration及每块registration给出source/code/math精确锁与命令；`prefix_checkpoint.json`和attempt登记可在真实故障后用 `--resume` 恢复，完成块由不可变checkpoint链追溯。每个后继进程重新载入源并重放已完成delta，核对选择前缀与逻辑frontier。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/homebrew/opt/python@3.15/bin/python3.15 -B /Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_framecover077/convert.py --work-dir /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_framecover077 --reconstruct
```

实际另PID67566执行重建，记录总21.060225秒、exit0且OS查无；从同一066源+2块delta得到完全相同完整state SHA。原运行与重建收据分别保留；没有在重建中再次计算资格或核。

独立078已获得第一块与最终冻结state/schema/hash进行独审；作者核验与确定性重建不替代独审。主保留数学及最终资格责任。只写新077目录和外部scratch077，未提交、推送、使用MATLAB或修改他人文件。
