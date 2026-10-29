# 第45轮复现

全部正式运行用本机Python3.15.0、NumPy2.5.4、SciPy1.18.1。需要标准库、NumPy及SciPy；不安装新软件。只做参数与独立holdout回放的精简入口可以仅用NumPy，见下文。

在仓库根运行，output必须为尚不存在的路径：

```sh
PYTHONDONTWRITEBYTECODE=1 /opt/homebrew/opt/python@3.15/bin/python3.15 centered-cube-lower-log/rounds/045/verify_saved.py --output /tmp/round045-verification-new.json
```

可换成具备依赖的Python路径；baseline wrapper要求冻结参数结果逐字相同，跨平台三角函数或NumPy改变时可能不满足该强重现条件，应保留差异报告，不能自动忽略。其余浮点复核的明确容差见代码；这些容差不是向外区间误差。

主入口先核对payload、PLAN/FROZEN、最小原来源快照哈希；在临时目录中恢复SUMMARY的三行RESULT并核对原字节SHA。随后实际重算20个Fraction模型、质量陷阱、Decimal80区间及全部profile、22组参数/10万来源样本、低维sampler/积分、144个有理窗口核和通用IS接口、半周期恒等回放、3×32768个冻结holdout的独立核实现。临时输出退出时清理，归档内容不变。

三份NEAR_INDICES原数组未归档。独立回放从seed重建111064项索引，按原UTF-8、json.dumps(indent=2)+实际换行的字节序列计算SHA，与NEAR_SHA256.json冻结值比较；不落盘数组。不需要旧轮代码依赖，也不需要本机原始论文目录。

若只需独立holdout回放：把round045_prior内replay_summary.py、PLAN.json、SUMMARY.json、NEAR_SHA256.json及三个FROZEN放入新临时目录，运行replay_summary.py。仅需NumPy，会生成REPLAY_SUMMARY.json。若同目录有原probe.py，额外核对其PLAN哈希；脚本不import它。这个回放不是新增独立样本，训练未重跑，raw h0舍入差异单列。

若需重新训练，应在新副本保留probe/is_engine、COST，移除PLAN/SUMMARY/FROZEN等已有输出后运行probe.py。那是重新实验，结果必须另存，不覆盖本轮证据。原脚本会输出near数组；正式精简归档使用其hash。

accepted/paper.tex是本地已有接受结果的只读快照，source_d10.npz是参数sampler的最小标签输入；numerical_report.tex仅用于历史资格门的来源对照。快照不是本轮新论文或定理重审。本轮不编译它们，不改原LaTeX编辑器。

验证不能证明：高n匹配GF basis已经构造、连续超水平集已认证、PRNG是精确iid，或C_n具有新的正对数下界。全部这些仍未完成。
