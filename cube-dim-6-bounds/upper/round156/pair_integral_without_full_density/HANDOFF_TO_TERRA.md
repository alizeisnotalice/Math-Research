# C156有限诊断交接：仅Terra执行，Sol未预检/求值

初查本clone round156不存在。此目录只完成理论与代码静审，
未运行prepare、数学预检、求根或任何数值。下界150仍USER_PAUSED。

唯一prepare入口（会复制四个源文件到Terra独立clone并冻结源）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round156/pair_integral_without_full_density/prepare_for_terra.py --root /Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round156/pair_integral_without_full_density
```

准备不导入Arb、不生成density、不执行诊断。它冻结130真cone矩阵
与kappa、128准确p区间、146完整真实beta2 cap/目标，绝不从显示
小数复原。旧全6D工作流与149/151 reserved density均不启动。

主与独审必须检查PROOF中导数、erfc Mills、survival log凹性、
normalized submultiplicativity、Cauchy四列seed及真实全域归一化。
parent_release.json需round="156"、authorized/main_approved/
independent_review_approved全true，及registration里的diagnose、
proof、inputs、source_manifest四个hash和registration自身hash。

之后仅运行diagnose.py一次，单进程nice10、四个BLAS环境线程1，
无wall/CPU cap。固定I0及一个sqrt(h)锚点，不搜索连续参数；最多
10416个原joint pairs、每pair35个caps，Cauchy sqrt按准确P缓存
并先作可证明的省略。新scalar calls为exp3、erfc2，复用旧准确p
不新erf；另有sqrt(h)、sqrt(pi)和必要Cauchy sqrt。

结果同时给最优two-column flat family和有限正卷积/两锚点族的
完整group-reduced beta2预算。budget interval的上端是实际mu2
上包，下端仅属于该预算表达式。若其下端仍高于7e-8/8e-8，只
排这个有限预算族；若上端改善旧cap，留主安排独立saved审计与
146消费，不能自行宣布C6。失败/非竞争后不自动加锚点、density
case、kernel网格、worker或候选。
