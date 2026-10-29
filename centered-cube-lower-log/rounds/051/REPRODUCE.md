# 第51轮复现

依赖：Python 3、NumPy、SciPy、mpmath；本机实际使用Python3.15、NumPy2.5.4、SciPy1.18.1、mpmath1.3.0。没有调用MATLAB，本轮全部数值由本机Python执行。

在本轮目录执行：

```sh
PYTHONDONTWRITEBYTECODE=1 python3 verify_saved.py --output /tmp/cube051-verification.json
```

它先验证payload及第45/46轮引用哈希，在临时目录复制必要文件，再实际运行六个命令并比对全部六个JSON结果摘要；不改归档结果。mpmath正常安装时无需附加设置；使用已有vendor时设置CUBE_MPMATH_VENDOR指向包含mpmath/的目录。本机该路径为`/Users/zhengzhihao/Desktop/T/output/hl_cube_paper_20261013/numerics/vendor`。

channel脚本可单独以CUBE_ORACLE_DIR指定第46轮oracle目录，CUBE_INPUT_DIR指定本轮channel/inputs。高维pilot须传`--high-pilot --manifest ../046/round046_gf/MANIFEST.json`。其几何比较先用Fraction，后以70位mpmath计算；绝不先把全局相位转float。

输入NPZ仅保存不可替代的标签；来源索引/原始labels哈希在channel/inputs中。小源4点和8点分段Gauss积分差、Monte Carlo profile及高维pilot都只是浮点证据。严格证明不依赖这些数值达到某个阈值，亦不以脚本复现替代误差证明。
