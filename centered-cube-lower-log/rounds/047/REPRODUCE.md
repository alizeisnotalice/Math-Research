# 第47轮复现

本机正式运行Python3.15.0、NumPy2.5.4、SciPy1.18.1。仓库根执行：

```sh
PYTHONDONTWRITEBYTECODE=1 /opt/homebrew/opt/python@3.15/bin/python3.15 centered-cube-lower-log/rounds/047/verify_saved.py --output /tmp/round047-new.json
```

output须不存在。主入口核对payload及第46轮证明SHA，在临时副本实际重算复倾斜有限和、单表最优值、96配置两seed的重要性实验；核对完整canonical行SHA，重算400个Gaussian密度比检查及全部96选点Hoeffding界，再生6640点profile压缩数组并逐项精确比较。仅回放原独立测试，不增加独立重复样本。不同数值库若不能精确复现，应保留差异，不放宽canonical检查。

verify/package支持省略v1全文；其历史耗时和hash仅从冻结INDEX引用。v1数据不参与v2推断。全部新中间RESULTS/SUMMARY、缓存留临时目录自动删除。

三种报告中的99%界分别收费：96个独立选点、全部profile、96个高阈值。不能把三套各自99%保证合称共同99%；本轮主报告的正增益仅用独立选点那一套。经验SE/ESS不提供正式置信保证。

严格命题见complex/PROOF与profile/NOTE；数值不是向外区间或真实HL尾概率证书。Gaussian代理与实际输入的迁移、稀有最大事件下界仍未证明。原论文/LaTeX不变。
