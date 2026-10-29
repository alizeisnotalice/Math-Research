# 第49轮复现

需要已有Python、NumPy、SciPy。本轮主代理正式重放采用Python3.15.0/NumPy2.5.4/SciPy1.18.1；minorant初次子代理运行用已有Python3.14，runtime_python只作运行记录。无需新安装。

仓库根执行，输出路径须不存在：

```sh
PYTHONDONTWRITEBYTECODE=1 /opt/homebrew/opt/python@3.15/bin/python3.15 centered-cube-lower-log/rounds/049/verify_saved.py --output /tmp/round049-new.json
```

主入口核对本轮payload及旧轮引用SHA，在临时目录映射第48轮协方差数据，实际重算7个真实尺度亲和度指数/概率有理帽、30个抽象有限模型、3宽度minorant全频求积和所有s的负性Fraction证书，以及96个多bank有限和案例。GRID中误差分数和也独立重算。旧轮、归档和编辑器不改。

Fraction值、指数下取整和有理概率帽精确核验。minorant的SciPy积分/优化与多bank NumPy有限和是数值证据；严格结论分别来自自含证明。不能把包络预算下界误称为实际CF误差下界，不能以Gaussian校准替换人工dec模型。

本轮严格弱型比较只针对已控制的有限共同尺度族，比较对象为同族各product来源的弱范数，不含对这些范数的维数上界。没有证明连续全尺度比较或新的C_n正对数下界。
