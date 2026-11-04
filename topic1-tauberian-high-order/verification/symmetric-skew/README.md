# R11：全对称斜晶格公式的有限参数证书

本包核验给定根方程与坏区上界的算术后果，不独立证明实际几何见证或坏区包含。这里采用精简公式

\[
B_{\rm single}=2A(2\delta)^{n-2},\qquad
B_n=\binom n2B_{\rm single}+\mathbf1_{n\ge3}\,8\binom n3(e+h_M)^3(2\delta)^{n-3}.
\]

删除单活跃对的第三坐标面层依赖单独的几何论证，不由本包有限计算建立。公式输入以及计算细节见 `compute.py`、`verify.wl`。

环境：Python 3.15.0（仅标准库）；Mathematica 15.0.0 for macOS ARM。以下从本目录运行，Python 生成时必须显式指定输出目录：

```sh
python3 compute.py generate --output-dir results
python3 compute.py verify --input results/certificates.json
wolframscript -file verify.wl "$PWD/results/certificates.json" "$PWD/results/wolfram-crosscheck.txt"
```

输入为 n=2,...,16、κ=16n/3、u=2^(-j)，j=8,10,12,14,16,18,20，共105组。与R10存在6个相同(n,κ,u)元组：n=2,3,4各有u=1/256和1/1024；其中n=3,4的行列式模型不同，n=2的两组是更高精度复核。

最终运行15个不同Python worker，约1.01秒。全部显式守卫与区间次序检查通过；按精确G区间判断，89组严格为正、16组严格为负、0组未决。负值均保留，不能将守卫通过误说成G均为正。最终记录为 `results/certificates.json`。

每组的 h* 与 hM 各经400次Fraction有理二分。根记录中十六进制整数 p 和 p+1 对应端点 `(δ−e)p/2^400`、`(δ−e)(p+1)/2^400`。正系数多项式及其积分在正半轴递增；按各根括界形成A、B_n和G的精确上下界，再对G向外舍入至35位有效数字。`verify`读取存档的两个根括界而不重搜，精确重建全式并检查守卫、端点顺序、符号和十进制向外包围性；105组实跑通过，见 `results-verify.log`。

Mathematica 对全部105组以180位工作精度分别求双根、以150位数值积分面积；d先保持精确有理数，避免直接以浮点计算1−d。最终15个子内核均启动；双根括界、G区间与符号各105/105吻合，最大区间中点绝对差约1.44×10⁻²³。逐组状态为 `results/wolfram-crosscheck.txt`，`results-wolfram.log`无消息。此过程是数值交叉检查，不替代Fraction区间证书。

计算没有设置时间限制。105组的通过与符号结果不证明任意参数或一般维几何，也不识别真实极大算子的下一阶系数。

其余必要核验从本目录运行：

```sh
wolframscript -file algebra.wl
python3 combinatorics.py
python3 witnesses.py --input results/certificates.json --output results/witnesses.json
```

`algebra.wl` 实跑以15子内核核对n=2,...,16的全对称矩阵行列式、二次系数与平衡指数计数，并核验一般参数下的优化及阶次恒等式；归档输出 `results/algebra.txt` 全通过。`combinatorics.py` 的n=2,...,9实际组合检查输出为 `results/combinatorics.txt`，一般维数分类仍由正文证明。

`witnesses.py` 读取同一精确根证书，选n=2,...,16、κ=16n/3、u=2^-20、N=u^-9的15个明确有限输入。精确有限比值下界全部超过L_n；孔中心与近角点共30个实际半径的平均值和严格裕量已用Fraction复算。孔中心另用Laurent系数递推与平衡指数公式交叉检查，n≤6额外直接求和完整小整数块。所有坐标、半径和精确分数保存在 `results/witnesses.json`。有限块之外的排除依赖正文几何，不冒充无限穷举。

主代理从归档路径实际重跑 `compute.py verify`、`algebra.wl`、全部105组Wolfram交叉检查和 `witnesses.py`，结果通过。Wolfram输出记录实际内核PID。未为保持CPU占用而重做已完成的批次。
