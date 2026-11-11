# 第96轮认证契约

本轮有三类不同证据：

- 正文证明一般物理源的必要信息下界、有限覆盖族精确对偶及双坐标块的上界；没有闭合一般端点。
- `information96.py`用精确有理数给Jensen上界，`logi`的正级数/显式余项认证对数，分数幂上界逐项化为整数幂比较。r=31/20的1024平方格证书给信息系数严格大于0.12324；mu的有理根括号通过log(mu-1)+mu的异号检查。
- `periodic96.py`的MILP只寻可行源。独立`check96.py`不相信求解器小数：按源格分数、矩形原始端点重算质量、严格余量及实际receiver并集，继而裁剪有限周期盒并结算实际Orlicz分母上界。所有solver dual、gap和最优状态都另列数值，不能解释为一般密度下界。

在项目根目录，用已安装Python（mpmath、NumPy、SciPy）和Wolfram15：

```sh
python code/information96.py /tmp/information96
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel wolframscript -file code/information96.wl /tmp/information96
python code/periodic96.py /tmp/periodic96
python code/check96.py --project-root . --output /tmp/validation96-replay.json
```

所有完整批次均通过`code/math_guard.py --reserve-gib 1 --max-rss-gib 6.5 --log PATH --report PATH -- COMMAND`运行；资源记录归档于`resources96.json`。这是采样监控，非操作系统硬保留。Wolfram任务每案256MiB、120秒，完成CloseKernels。Python搜索12个spawn worker，每worker求解器线程1；120s time limit只证明已有incumbent可行，不证明最优。纯网格cell比例用坐标1上的真实薄片实现，矩形包含整receiver cell，边界零测度不影响体积。

有限原源：把周期源截到[0,K]^N，K=1000N，每个周期见证每边至多1；计[1,K-1]^N内的整周期receiver，给因子(1-2/K)^N。输入1_E/p的分母为|E|p^{-1}(1+log(1/p))^{N-1}；对数取定向上端给明确下界。此归一化没有省略原Orlicz价格。

Mathematica符号残差4个均为零，7个分数阶矩以80/160位目标独立积分；有限积分只是对照，不替代Jensen有理上证书。验收同时检查原始有限LP价2/3及Gram松弛价1/2，以及源码、数据、正文和前轮依赖SHA256。
