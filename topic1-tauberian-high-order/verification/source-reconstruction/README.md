# 第63轮回源重建的有限计算

程序按原字节保存；文件名作了归档重命名。实际输出已分别保存，程序不证明覆盖、全部输入或几何极值。

- `foundation_checks.py/json`：标准库Fraction/comb，n=2…12的ordinary第七系数、n3第五系数零及q=2…24的Young指数。
- `original_checks.py/json`：标准库Fraction，有限中心错误推理的反例、严格等分余量及周期长度恒等式，10项通过。这不是指控原PDF使用错误推理。
- `sixth_checks.py/json`：SymPy1.14，六层消项、Hölder指数、正支凸性、标量极大值、内部质量和周期几何，10项通过。该脚本按原程序写`algebra_results.json`，请在临时目录复现；初次失败输出保留在仓库外。
- `eighth_checks.py/stdout`：标准库Fraction自写有理级数倒数；上界用实际分母，下界用完整斜格行列式，n=2…12及二维微预算。原stdout未经裁剪。

可在临时目录复制相应脚本后运行`python3`；SymPy脚本需已有SymPy环境。实际使用环境和局限在四份解析重建中记录。每份报告的数学证明与有限程序证据分开，未把旧审核判语代入计算。

`manifest.json`绑定所有归档程序/输出及四份冻结重建。PDF的字节哈希、初始提交与原稿/53页版出处另见`evidence/063-snapshot.json`。
