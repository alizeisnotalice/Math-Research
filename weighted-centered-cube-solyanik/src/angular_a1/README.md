# 第28轮：角向A1分类与对数余项

主稿paper/small_homogeneity_weak.tex整合并替代第27轮的有界角向/粗幂余项。零次齐次正轮廓a且0<Ea<∞时，某个τ的完整弱型有限，当且仅当a∈A1(K)，这也等价于所有0<τ<n的弱型均有限。经典A1必要子集检验和本稿的齐次两尺度比较均给完整证明，不把经典A1理论称为新增结果。

SW.1对中尺度按半径比s_j=q^j/2逐壳求和，SW.8改进为统一O_{n,R}(τlog(1/τ))（R控制[a]A1）。q=e^(1/n)还给真实范数的对数双侧界：相对值位于1−τlog(1/τ)−O_n(τ)与1+(e−1)τlog(1/τ)+O_{n,R}(τ)之间；不称精确第二项或锐系数。

## 无界显式例及数值诊断

cube中aγ=(ρ/|x1|)^γ,0<γ<1,n≥2；Ea=(n−γ)/(n(1−γ)),a_*=1，角向A1有限，故τW→(n−γ)/(1−γ)。程序unbounded_angular_cube_batch.py以标准库Decimal120位、15进程计算32唯一(n,γ,τ)，参数见run_summary.json。原子在e1，质量2^n(1+τ)^n，halo是e1+(1+τ)K；按x1切片正积分解析求值。γ+τ=1时显式使用对数分支。签名含源码、环境、模式与参数，旧缓存隔离；本机首次与续跑均已运行，完整表/日志/检查点不重复归档。

unbounded_cube_independent.m在MATLAB R2026a实际独立积分32行；Python/MATLAB最大τR差1.776e−14。unbounded_cube_guards.wl核8组角均值、32行正性及对数分支。数值仅验证这个输入族，不能替代完整算子上界。

angular_shell_guards.wl独立符号核Ea积分、尖峰均值385/384、任意跨零区间不可漏因子2、壳比优化与两个对数极限；另以30个精确有限和核对二项式/几何级数式。失败均非零退出。它核对证明中的恒等式，不认证全域范数最优。

辅助A1控制初稿漏了一维跨零区间因子2：对称区间的1/(1−γ)不能用于任意截面。已改成2*3^abs(n−τ−γ)/(1−γ)，改变源码签名后重新计算，原子数值不变。所有失败与修正保留在本机日志及最小摘要中。

## 复现

```
python3 unbounded_angular_cube_batch.py
wolframscript -file unbounded_cube_guards.wl
wolframscript -file angular_shell_guards.wl
```

MATLAB运行unbounded_cube_independent.m。有限角向均值及正下界仍不保证A1：正文的分离尖峰反例有Ea=385/384、a_*=1而全部弱型范数无穷，证明完全解析，不把有限截断实验当无穷反例。
