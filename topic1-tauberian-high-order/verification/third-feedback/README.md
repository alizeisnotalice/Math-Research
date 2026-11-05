# 第三反馈与微孔的最小实际计算

```sh
wolframscript -file third_feedback.wl
wolframscript -file cap_l1.wl
wolframscript -file single_corner.wl
python3 fraction_polynomial.py
```

实际成功环境为 Wolfram 15.0 ARM；日志保留 PID 32847、33245、34537。`third_feedback` 核验三次线性多项式算子、241/60与参数归一化；`cap_l1` 精确排除标量 min 违反集并核两条带的系数；`single_corner` 只核解析几何供应后自由体积和参数预算的代数。

`fraction_polynomial` 独立使用标准库有理数重构三次算子及低维系数。有限精确核验不替代实际集合、所有源点、参数稳定和有限外盒的证明。早期partial与冗余有限样本不归档；没有为CPU占用重跑旧批次。

文件哈希及两份冻结证明的绑定见`manifest.json`。几何结论仅为真实有限构造间的比较，不能用这些输出宣布 C_n 的绝对九阶常数。
