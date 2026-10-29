# 复现与证据说明

普通核验需要Python、NumPy：

```sh
python3 verify_saved.py
```

它验证原件哈希，在临时副本重放：12类PB log-DP扫描、6个Fraction trace例、4个异表Gram/trace例及2个直接L¹例、连续-bin差额、独立909组有理倾斜参数及有理反例、24组coupled Gaussian积分及一组细化。原始结果不覆写，不抽尾MC，不建大缓存。

完整核验另外需要SciPy，当前本机已验证的解释器如下：

```sh
/opt/homebrew/opt/python@3.15/bin/python3.15 verify_saved.py --full --output /tmp/cube-round039-verification-new.json
```

`--full`另重放channel/probe.py的Gaussian、binary及38轮一维真实核诊断，以及Gaussian大R三候选积分。历史source只读取自`../038/round038_highdim/source_tables.py`，在临时目录映射为原路径。需要保留本仓库的038目录。不会建立38轮的大维特征缓存，也不重跑其MC。

原agent：Python3.15.0、NumPy2.5.4、SciPy1.18.1；主线：Python3.12.14、NumPy2.3.5。数值重放忽略运行时间、解释器文字和NumPy版本字段，数字按明示浮点容差比较；Fraction字符串及整数精确比较。它不是完整FP向外证书。`--output`要求目标文件尚不存在。

posterior/probe_initial_linear.py和EXECUTION_NOTE仅用于追溯一次远尾下溢，默认不会执行已知失败的线性算法。已改为log-DP，输入case没有删减。该故障与证明中严格反例分开。

Gaussian模型的CDF/quantile、输出网格和源Gauss积分都只是数值探路。其一般阶数和耦合模型L¹极限由自含解析证明给出。主结果是固定尺度L¹稳定性；本轮没有将Gaussian熵阶移植为C_n的对数改善。
