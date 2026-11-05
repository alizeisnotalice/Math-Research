# 保护算子与限定平面积分的实际证据

```sh
wolframscript -file plane_exact.wl
wolframscript -file h2_spec.wl
wolframscript -file i2_branches.wl
wolframscript -file i2_check.wl
python3 i2_fraction.py
wolframscript -file tonelli_shortcut.wl
wolframscript -file tail_exponents.wl
```

`i2_branches`加载`i2_common.wl`，生成`i2_results.wl`供`i2_check`读取。实际Wolfram15 ARM批次有15个不同子核；十二个矩和三个接点是不同任务。子任务CPU累计10.767秒，WL MaxMemoryUsed峰112557304B；这不是外部RSS，未采得运行中CPU/RSS，不据核数声称持续满载。

I1/I2、完整近零min分支及空间指数是限定平面证据。独立Fraction只核指定三个有理参数的四分支积分；一般参数恒等式及真实几何由解析证明承担。没有运行整套区间恢复或用这些日志证明C_n全阶。

`v2_fibre_patch.json`保留单段坐标明确化的old/new及原/v2哈希，能从最终公开证明精确逆恢复原冻结。原日志编码伪影照实保留。`manifest.json`绑定必要文件及最终源/独审。
