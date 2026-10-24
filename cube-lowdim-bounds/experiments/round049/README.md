# 第49轮复现材料

来源： https://github.com/tejasvi/square-maximal-weak-type-v4 ，固定提交 bca281a706bbf6375b04d601ad19ba4ad0bf28b6 的 lowerbound。两源文件及配置原样保存，许可沿用 ../round020/upstream-v4/LICENSE。Python 3.12；第二程序需要 NumPy。共同使用 Python 整数、Fraction、isqrt；面积分别使用线段树扫描及坐标压缩差分，不能称算术原语完全独立。

复现命令（在本目录执行；第一程序产生 config.check2.json）：

```sh
python -I check_exact2.py config.json --Dr 1099511627776
python -I lb_check2.py config.json 1099511627776 482167717377842991882283179477574801392641890014291050323/283174368139841919615520431313000857449005056000000000000
```

本轮实际各执行一次，180秒独立硬预算；未超时。日志保留实际运行时本地路径，不把路径当作可移植依赖。verified.json保存准确质量、面积及一致性检查结果。完整矩形dump可再生，留本地而不重复上传。六例并集检查覆盖空、相邻、交叠、包含、重复和零面积矩形，小整数单位格穷举作第三参考。

union_cases.py是实际六例测试的可移植副本，仅改源文件定位和输出文件名；归档后未重复运行。运行 `python union_cases.py` 可另存测试结果。
