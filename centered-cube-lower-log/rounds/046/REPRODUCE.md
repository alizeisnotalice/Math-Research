# 第46轮复现与边界

正式工具为既有Python3.15.0、NumPy2.5.4、mpmath1.3.0与/usr/bin/clang。无需SciPy或MATLAB；本轮没有安装软件。

在仓库根目录执行，output必须尚不存在：

```sh
PYTHONDONTWRITEBYTECODE=1 /opt/homebrew/opt/python@3.15/bin/python3.15 centered-cube-lower-log/rounds/046/verify_saved.py --mpmath-vendor /Users/zhengzhihao/Desktop/T/output/hl_cube_paper_20261013/numerics/vendor --output /tmp/round046-verification-new.json
```

若Python正常安装mpmath，可省略--mpmath-vendor。若使用其他本机既有vendor根，替换该路径，要求其下有mpmath/；依赖不会复制入仓库。窗口程序也支持CUBE_MPMATH_VENDOR环境变量，并以正常import优先。所有source数据显式传入，不依赖原论文目录。

主入口先检查payload及第45轮accepted/paper.tex、source_d10.npz的引用哈希。随后在新临时目录编译C源码，实际重新执行冻结seed搜索（首候选成功）、显式basis再验、Python/C64组算术、原小域负例和成本probe、trace矩阵，以及独立NumPy全量约210万乘积/419万比率核验。搜索输出的时间字段不要求相同；其他manifest字段及basis均对照。之后在临时目录恢复原冻结manifest字节，用固定输入哈希进行窗口/安全节点回放。不会覆盖归档。

其他实际重算包括：1040个安全尺度节点、8184个Fraction clip点、直接周期峰积分、零均值分段积分、非零Fourier交项、大域有限点值与trace masks、b1至11 Walsh扫描、36个数值Gauss和及8个小域联合log-score MGF有限和。大n的4行误差是解析界求值，未模拟n个坐标。

窗口adapter的显式入口：

```sh
CUBE_MPMATH_VENDOR=/path/to/existing/vendor python3 validate_oracle.py --labels /path/to/source_d10.npz
CUBE_MPMATH_VENDOR=/path/to/existing/vendor python3 validate_oracle.py --manifest /path/to/MANIFEST.json
```

输出写到脚本同目录，故这些单独入口应在临时副本执行。支持CUBE_SOURCE_LABELS/CUBE_GF_MANIFEST环境变量。主入口已自动使用临时副本。

归档不含编译后的二进制、mpmath源码副本、全部K标签、products/ratios数组、worker重复SHA清单或可替代的旧结果。GF冻结计划成本倍数笔误已在REPORT明确勘误；仍保留原计划以维持可追溯性。

精确整数和Fraction断言属于精确计算核验。70位mpmath、Gauss积分与NumPy有限和不属于向外区间证书。联合矩的严格结论来自自含证明；没有从数值一致性替代证明。已接受主定理未重新核验。全增长尺度最大事件的下界和新的log因子仍未完成。
