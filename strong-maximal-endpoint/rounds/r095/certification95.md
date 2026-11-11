# 第95轮有限认证契约

严格高维结论见[核心](core95.md)和[独立审查](audit95.md)。软件只检验有限概率算术，不把单点覆盖模型称为完整高维物理源，也不声称有限生效维数。

Python用`fractions.Fraction`枚举三色各两副本的64个覆盖状态；独立Wolfram精确求和给同样的实际存活质量、损失及有序成对交叠上界。推广取m=1…16、k=1…8、p=1/31或1/127，共256案，核对

\[
 mkp-m\{1-(1-p)^k\}(1-p)^{k(m-1)}\le mk(mk-1)p^2.
\]

只保留规范行哈希及极值。c=1/1024、eta=1/128、delta=1/1024的有理余量证明各相交好事件的概率下界为正；守护沿长边的最低删除比例3eta/2严格超过eta。周期密度上下界的对数因子只核对符号顺序，真实双向等价依赖正文裁剪与随机平移证明。

项目根目录，使用已有Python和Wolfram15：

```sh
python code/experiment95.py rounds/r095/python --workers 2
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel wolframscript -file code/verify95.wl rounds/r095 1
python code/check95.py rounds/r095 rounds/r095/validation95.json --project-root .
```

完整计算通过`code/math_guard.py --reserve-gib 1 --max-rss-gib 6.5 --log PATH --report PATH -- COMMAND`执行。Python资源报告归档于`python/resources95.json`，Wolfram于`wolfram/resources95.json`；主代理最后在最终code路径重放。Wolfram的这组有理关系使用主核，不额外启动空并行核；进程结束退出。本轮的长时间12核负载由独立周期网格搜索承担，不以反复计算同一常数制造CPU负载。MATLAB未启动。

验收按独立有理表达式、Python/Wolfram规范行哈希、边界余量和资源读数核对。源码、输入、正文和前轮依赖SHA256记录于[validation95.json](validation95.json)。其中哈希只表示文件一致性，存在性反例和完整体积结算必须阅读数学证明。
