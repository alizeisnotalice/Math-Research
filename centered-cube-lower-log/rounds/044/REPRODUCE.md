# 第44轮复现与证据范围

需要Python、NumPy、SciPy。所有本轮正式运行使用已有本机Python3.15.0/NumPy2.5.4/SciPy1.18.1；没有安装新软件。MILP由SciPy连接本机HiGHS执行。

```sh
python3 verify_saved.py --output /tmp/cube-round044-check-new.json
```

默认入口只读核对全部payload与旧helper SHA，在临时视图中重算：反射数值扫描、K2/K4低维联合积分、精确尺度预算、Fraction镜像/质量/核/联合谱检查，及MILP所有已保存incumbent的训练覆盖与完整8192点heldout。同时重算快速贪心控制及其冻结结果。不会重求MILP、不会重新训练mirror selector。输出目标须尚不存在。

```sh
python3 verify_saved.py --replay-mirror-test --output /tmp/cube-round044-full-new.json
```

额外重放三个mirror维数的全部冻结测试，并重算短/长窗见证分布，约20秒，无新独立样本。每个FROZEN的质量/plan关系在默认模式已核对。

旧helper不在本轮重复归档。核验器把仓库`rounds/043/round043_binary/probe.py`映射到临时同级`round043_binary`目录，供mirror通用训练器与MILP旧来源引用；把`rounds/042/round042_selector/probe.py`映射到临时同级`round042_selector`目录，供K4联合积分使用。原字节SHA必须匹配。MILP原PLAN中的source_path是实际执行时的历史路径，不作为唯一复现路径。

不要直接在归档文件夹运行会写结果的原入口。完整新训练需复制三个本轮子目录及上述旧helper到一个新临时根，在副本中移除对应PLAN/SUMMARY/FROZEN输出再运行。mirror入口为`probe.py scan`、`probe.py search`；MILP入口为`probe.py --cost`、`probe.py`，再运行`greedy_control.py`及`verify_saved.py`。这些新训练是独立重做；time-limit MILP在不同机器/版本上可能得到不同incumbent，必须另存status/gap，不能保证原解字节一致。

本轮归档保留最终证明、冻结来源/选择器、各项最小结果及模型状态；不保存完整DP特征矩阵、query数组、solver全日志、重复manifest。低维grid积分不是向外区间证书，部分K4粗网格的质量偏差会使积分小于1，不能把它当负χ²。高维guard测试和CP区间只描述浮点分类器，尚未认证真实连续超水平集。

理论主结论来自自含证明：固定共同长窗的特殊计数门控族联合χ²界，以及稀疏长尺度比较。它不是全f强L1误差估计或一般HL上界；短窗、连续密集尺度和C_n正对数增益仍未完成。
