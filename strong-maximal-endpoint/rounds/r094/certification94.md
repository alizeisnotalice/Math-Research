# 第94轮认证契约与复现

严格数学证明见[核心](core94.md)与[独立审查](audit94.md)。有限软件证据分别登记：斜带证书为精确有理数；4个驻点 Hessian/分解残差为符号恒等式；4225曲率点、24倾斜积分、315正混合CDF及30848随机任务为有限数值检查。没有给普遍定理的有限起始维、一般指标上界或新的一般下界。

Wolfram最小记录保留完整24积分以及每64个指数网格选一个，共67曲率点，并登记全4225原始输出SHA256。Python独立重放这些67点及24积分。全量随机及正混合行留在仓库外，Git只保留参数、极值、精选行、资源和哈希。独立验收不能从一个哈希推断未归档原始行的数学性质；需要时按最终源码完整重算。

在项目根目录，使用已安装的Python环境（mpmath、NumPy、SciPy）与Wolfram15：

```sh
python code/laplace94.py
python code/normal_probe94.py
python code/tilted_probe94.py /tmp/r094-tilted-replay
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel wolframscript -file code/verify94.wl /tmp/wolfram94 12
python code/reduce94.py --output rounds/r094
python code/check94.py --project-root . --output /tmp/validation94-replay.json
```

完整批次需用[资源守护](../../code/math_guard.py)执行，分别把报告写到`/tmp/laplace94_resources.json`、`/tmp/normal_probe94_resources.json`、`/tmp/r094-tilted-replay-resources.json`、`/tmp/verify94_resources.json`，归档器据此结算。例：

```sh
python code/math_guard.py --reserve-gib 1 --max-rss-gib 6.5 --log /tmp/verify94.log --report /tmp/verify94_resources.json -- env WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel wolframscript -file code/verify94.wl /tmp/wolfram94 12
```

`laplace94.py`默认原始目录`/tmp/laplace94`；`normal_probe94.py`默认`/tmp/normal_probe94`。每任务BLAS线程1，外层12进程；Wolfram12核每任务限256MiB、120秒，结束CloseKernels。守护只终止自己的进程树；不操作其他项目。315案有最大维42案先行试点，恒等分支仅`a=1,xi=0`；相同权重而不同行动仍使用真实两分支CDF。SciPy截尾基于CDF/SF的浮点评估，不包装为区间证书。随机对照的数值优化器状态、容差、种子与概率账本均保留。

当前`check94.py`核对源、输入、前轮依赖与证明/报告SHA256，精确重算有理平均余量及完整halo体积差下界；使用独立mpmath部分分式表达式检验Wolfram倾斜积分。此验收证明这些有限关系和文件的一致性，渐近全法向结论仍依赖正文统一估计。
