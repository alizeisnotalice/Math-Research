# 复现入口

在038目录使用已有Python和NumPy。默认核验：

```sh
python3 verify_saved.py
```

该命令检查归档哈希，在临时副本重算12个精确有理小例及分层平方根包，核验已保存PB数据、288个峰端点、比例/CDF误差传播和336项CI。原件不写入；不抽新高维样本、不建立大缓存。它不能从汇总恢复原逐样本RN或paired协方差。

完整重放另选尚不存在的绝对输出目录：

```sh
python3 replay.py highdim --output /tmp/cube-round038-new-highdim
python3 replay.py source --output /tmp/cube-round038-new-source
python3 replay.py conditional --output /tmp/cube-round038-new-conditional
python3 replay.py pb --output /tmp/cube-round038-new-pb
```

highdim模式运行原冻结代码，约需1.6GB RAM；原主运行13.25秒，耗时取决于本机BLAS及内存状态。拒绝覆写已有目录。重放目录复制必要小文件，并映射原程序所需的034/037依赖，保持冻结源码及哈希不变。代码会写新注册、新结果；不要把新注册的计时或样本流称为原次证据。

主运行：Python3.12.14、NumPy2.3.5；部分一维及小例检查原用Python3.15.0、NumPy2.5.4。跨版本随机multinomial序列不承诺逐位相同。PB模式的随机数只用于固定seed生成显式PB参数，不是原高维尾试验的追加样本。

round038_highdim的source、builder、driver及plan/registration/results均保持原始字节。source验证保留一次reference权重修正的原/新脚本与结果，背景模型未改。conditional/read_saved_stats.py与pb/check_bounds.py、highdim_small_checks.py是独立审查入口。

普通浮点计算尚无完整区间认证。精确Fraction证据仅覆盖所述有理模型/峰部分，不能扩大到全部三角函数核、BLAS或高维CI覆盖率。
