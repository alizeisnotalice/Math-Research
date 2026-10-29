# 第48轮复现

本轮实际使用已有Python3.15.0、NumPy2.5.4、SciPy1.18.1及/usr/bin/clang，不安装软件。仓库根运行，output必须尚不存在：

```sh
PYTHONDONTWRITEBYTECODE=1 /opt/homebrew/opt/python@3.15/bin/python3.15 centered-cube-lower-log/rounds/048/verify_saved.py --output /tmp/round048-new.json
```

主入口核对payload与旧轮引用SHA，在临时副本映射第46轮C源码和域manifest。实际重放冻结四表搜索和显式basis再验，再用第46轮独立NumPy实现对四表各自完整验证约210万乘积/419万比率。重算所有秩、Fraction协方差及认证区间，检查9位十进制端点确实向外包含其有理端点。所有旧文件只读。

随后重跑48配置、两个seed各12000条OU代理路径，独立密度检查及全部profile数值逐项精确对照，重算局部粗糙度小增量和主代理稠密协方差求逆检查。它们是原独立holdout的确定性回放，不是新增独立实验。不同数值库若产生canonical行hash差异，应保留差异，不悄悄放宽为新实验。

数学范围：banks/PROOF.md给实际来源在所列尺度上的HH精确式和HB/BB有理包络；roughness/PROOF.md是固定Gaussian OU模型的局部恒等式。process的经验Bernstein解释假定理想iid及精确被积函数，浮点/PRNG未向外认证。全部内容均不提供新的C_n对数下界或连续尺度概率迁移。

不提交可再生全RESULTS、编译二进制、缓存、完整周期标签或产品/比率数组。PROFILE_MIN为最小按阈值统计，非原路径缓存。个体搜索尝试信息可由固定seed重放，manifest保留最终资格统计。原论文/LaTeX未修改。
