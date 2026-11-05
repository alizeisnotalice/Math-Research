# 任意有限阶条件传递：实际计算证据

这些文件核验有限代数，不证明真实单线性任意阶渐近。

```sh
wolframscript -file recurrences.wl
wolframscript -file composition_batch.wl
wolframscript -file ell_wl.wl
python3 composition_fraction.py
python3 ell_fraction.py
python3 independent_fraction.py
```

`common.wl` 是前两份 Wolfram 脚本的共同定义。实际成功版本为 Wolfram 15.0 ARM；`composition_batch.log/json` 保留 15 个实际子核 ID/PID、54 组和 2700 项结果。批次过短，未采得运行中 CPU/RSS；不据核数推断利用率。初次受保护变量警告已修正并实际重跑，失败日志保留在仓库外。

`composition_fraction.py` 用广义二项式；独审的 `independent_fraction.py` 改用形式 exp/log 和 w 变量直接复合，核验 7088 项复合、144 项一般幂和200项基线。`ell_wl` 及 `ell_fraction` 另核105项基线。q=0 只是辅助连续边界，并非一个实际无限重数模型。

完整有限范围、解析证明边界和共同小 o 论证见本轮证明与独审。`manifest.json` 绑定发布文件与两份冻结稿；保存输出为真实运行结果，没有为利用率重复旧批次。
