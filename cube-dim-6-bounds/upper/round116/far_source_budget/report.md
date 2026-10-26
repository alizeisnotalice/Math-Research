# 116 — 真实远区源事件、严格预算障碍与已付大时尾

本轮封存，两个固定配置已取得完整证书。初次执行失败及显式 v2 输出修复均保留；不把初次运行改记为成功。未产生新 C6 上界。

真实区域固定为059的 F={x=(1,y), y∈[0,1]^5, dist(x,S2)≥1/64}，包含原31 unsigned plane距离刻画的完整补集；S2仍是全体实际 body-pair planes。完整111+108账给出的剩余 I 充分额度使用 upstream.json 中原始有理端点，R_L≈4.53710531716404e-6、R_U≈4.81959607749471e-6。全域 g0=1/65000 仍未证明。

## 完整推导及结果

proof.md 给出新正 coarea 包：far fiber 至少三条 |a|>β=ρ/(60√6)，且首坐标之和为1迫一条 |a|≥1/32。按最大绝对值参数的唯一标签分配后，保留其余31条的最大值 M31 和至少两条超过β的事件，得到

    ∫_F p_t ≤ E[c_t(max{1/32,M31}) 1{Nβ≥2}],  全部 t>0。

独立性仅用于解除空间限制后真实31个自由输入的期望，从未假设给定fiber后的独立性。五个自由坐标反射因子1/32作用于完整far measure，与32分支和相消；没有额外无交叠假设。短时此质量包 O(t³)，配完整真实投影点包后可积；大时用108六维峰。详细全域、端点、短时及大时范围见proof.md。

配置1只考察预登记固定窗[1/64,1/32]，无时间格。在窗上，四条真实 full-law 点包的最小值都>1/5，完整face质量包、binomial包、unique-max包和任意合法体积包都>7/1760。因而新显式预算表达式 Ienv 满足

    Ienv^5 > (1/64)^6 (7/1760)^5 / 5。

严格有理比较给 Ienv>1.96053862365116e-5，超过 R_U 约4.067848409倍。**被排除的是本轮指定的源质量＋峰插值预算族，不是 actual I_F、真实Green模型或C6。** 因严格分离，未加密时间、改变ρ/阈值或刷新旧函数。这里丢失的关键是选中最大系数与其余参数和的仿射关联，以及自由参数在五维face上的联合像；计数与最大标签信息仍不足以支付该中间时间窗。

配置2用已接受108的原law六维峰 P6=216/(2^30π^6)，合法替换095远区构造的 det160 basis峰1/(160π^6)，新旧比135/4194304。固定T=1/8，真实完整far区域满足

    ||∫_(1/8)^∞ p_t dt||_(6/5,F) ≤ P6/(5(1/8)^5)
                               <1.37130662983340e-6。

其严格精确上端小于 R_L^(5/6)/20；若单独取q次方约9.21640372589917e-8。结合短时必须先在同域加范数再取q次方，不能把这个单独I尾与短时I直接相加。留给实际0<t<1/8部分的充分范数额度约3.38914059262137e-5，完整原始端点见results.json。这只是已付款partial tail；完整far尚未付款。

## 实际运行与失败保留

原run_prepared.py保持SHA3e987c34...不变。唯一初次run，PID67690，tool chunk6eda4e，exit1：cfg1已完整原子checkpoint；cfg2所有数学断言及Arb crosscheck通过后，将巨大精确残差转字符串时遇Python默认4300位限制。原run_receipt.json保持FAILED/completed1。其CPU采样0.039888s、wall采样0.05741945799672976s；采样在最后收据写入和解释器退出之前，原执行没有独立完整外层wall证据。tool原报wall另保留，不把它当成完整科学耗时。

主明确批准一次输出限制修复后，新resume_tail_v2.py先set_int_max_str_digits(0)，仅读cfg1原checkpoint，调用原cfg_tail一次；cfg_window调用0，无数学参数变化。v2科学PID68172，observer PID68164，tool chunke4b5c2，exit0。独立observer从spawn前计至科学子进程最终receipt、解释器退出及reap后：wall0.08505600001080893s，CPU0.066779s。observer自己的收据写入/退出不在该child范围。原失败不改记成功。

check_saved_v2.py已实际只读核验，tool chunk71b1e6、exit0、PASS_SAVED_EXACT_ONLY_V2_LINEAGE；0新科学配置。check_saved.py是原封预运行checker，仅保留其原始pin证据，不用于续跑后的status；最终核验入口是check_saved_v2.py。

## 核验与范围

使用本机既有数学venv、Fraction精确算术与Arb384新尾交叉核验。所有旧projection/peak常数均从原保存数据读取，0旧peak/native/074-D网格/096重放；cfg1一次完成，cfg2初次未存完整结果后仅一次显式续跑。全部12项来源pins结束时未变。最小保存核验命令：

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B check_saved_v2.py

full far0<t<1/8、全域g0和完整 S=64(g_L/6)^(6/5)-I_U 均仍未闭合。原生成元、全空间正源机制以及001/027的一般L1、全部正共同尺度接口保持；无新C6结果。111及其结果未改，没有新方向、下界、LLM代理、MATLAB、stage或push。
