# 第50轮复现与范围

正式工具为已有Python3.15.0、NumPy2.5.4、SciPy1.18.1；不安装软件。仓库根执行，输出必须尚不存在：

```sh
PYTHONDONTWRITEBYTECODE=1 /opt/homebrew/opt/python@3.15/bin/python3.15 centered-cube-lower-log/rounds/050/verify_saved.py --output /tmp/round050-new.json
```

入口核对payload和旧轮参考SHA，在临时副本映射48轮manifest。实际重算9组joint预算（含失败项）、6对真实同端点rank、mod距离/phase断点组合测试、16个精确有限原子模型、32个归一化Gaussian源门控探索及3个实际截断窗口/边界积分。旧轮和归档不改。

Fraction预算和有限模型是精确计算；Decimal只用于显示。Gaussian门控、截断网格积分是数值探路，无向外证书。门控grid选点没有独立holdout，不能作已确认的HL收益。渐近小量来自明确假设的证明，不应当作当前n对全部poly(v)表数的小常数。

核心比较属于共同域/共同Fourier/joint基条件、固定权重及分离有限增长网格。截断桥只转移同一有限尺度max，不能上界完整Euclidean HL算子。没有证明单表O(√n)上界或新增C_n正对数下界。
