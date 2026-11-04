# 第25轮：有限原子研究与可复现审计

主稿paper/single_atom_lower_bound.tex给出经典幂函数系数在本稿外层弱型模型中的单原子测试，并严格证明W1/8>27/4。单原子子问题的最优参数不等于完整Wτ的最优。一般τ的标量最大化已有文献，不作首次声明。

## 有限多原子探索的历史证据

historical_candidates.json保存32个最终原子向量、必要配置/计数及分数；未保存优化轨迹或可再生全部检查点。historical_manifest.json与history/pre_signature_20261105/三份原始源码关联，后者为解释无签名历史不可替代的源码快照。原始完整JSON与运行日志仅保留在本机，manifest中的raw summary/log路径明确不属于精简Git包。

历史三族共7083648次目标评估，实际4 worker；所有记录均达到maxiter，不能当收敛或全域最优证书。τ只为1/4、1/2；这些值未超过当时双原子基准，与新τ1/8单原子下界无冲突。数值盒及种子完整见search_parameters.json。

在本目录运行`python3 audit_search.py`，独立重算32个最终候选的光晕积分与成本，核对每条评估数=人口规模*(maxiter+2)、总数及终止原因；另以200组Fraction小型输入比较连续块与全部子集光晕。不同粒度记录不能相加成新的独立证据。

## 当前签名版重算

需要Python3.15、NumPy2.5.4、SciPy1.18.1（历史实际版本）。三个程序依赖同目录run_provenance.py和multiscale_search.py的公共几何函数；从本目录逐族运行，避免同时启动三套worker池：

```
python3 multiscale_search.py --workers 15 --maxiter 900 --popsize 12 --seeds 2 --ks 6 10 14
python3 asymmetric_multiscale_search.py --workers 15 --maxiter 1000 --popsize 12 --seeds 2 --ks 8 12 16
python3 free_many_atom_search.py --workers 15 --maxiter 1200 --popsize 14 --seeds 2 --ns 12 16
```

这是当前代码的重算入口；历史命令确为4 workers，保留而不改写。当前默认15，实际并发不超过任务数。签名包括源码、导入的公共源码、Python/NumPy/SciPy、任务和参数。无匹配签名的旧缓存移入stale后重算；没有把旧708万记录伪装为新签名结果。

归档入口已实际做三族小规模新算（各2任务）、参数不匹配后重算与匹配缓存内容保持测试；这些只验证计算机制，不新增708万研究次数。完整搜索可不限时运行，日志和检查点在本机保存。

## 单原子验证

进入single/后运行：

```
python3 single_atom_search.py
python3 single_fraction_certificate.py
python3 single_geometry_verify.py
wolframscript -file single_atom_guards.wl
```

15个不同τ的临界点数值搜索使用15实际进程；几何积分独立核25行。Fraction以平方/八次幂的严格有理比较证明固定下界，single_fraction_certificate.json为最小证书。Wolfram四守卫已实际通过；其中曲率守卫覆盖0<τ<1，一般唯一性依赖正文凹性证明。MATLAB执行single/matlab_single_atom_independent.m，独立15行；实际环境和误差见run_summary.json。

本目录weak_atomic_sanity.wl/.m仅复核旧双原子基准与零点小区间成本发散。旧值可比性与新下界不共享未经证明的最优性假设。较弱变质量双原子中间路线未进入正文或Git代码，删改理由在本轮报告。
