# 124：rank-peeling 单阈值矩模型的严格认证

固定八个 h 的唯一有限批次及唯一纯保存核验均通过，现封存停止。
未得到123的竞争门槛证书；实际 M、多 h 模型均未被排除。

## 主结论

最佳候选在 h=1/256，给出真正矩上包 M<0.068475443529。
完整 M^6 上端及由123得到的完整 I^5 上端保存在 h_4.json，
不由显示小数构造。该上包涵盖全部正时间和整个正 face，但未通过

    Theta = 36893488147419103232/629200415790081024169921875。

八个固定 h 的固定 B 表 Holder 族均被精确 dual 下证严格排除。
其中最小的族最优 M^6 上预算下界仍超过 (7/4)Theta；完整准确余量
保存在 summary.json。这是指定上预算表达式最优值的下界，不是
actual M、actual I 或 C6 的下界，也不排除其他 h、更锐混合矩或多 h。

|h|目标评估|候选 M 上包展示值|族 dual 下界/Theta 展示值|
|---|---:|---:|---:|
|1/4096|43|0.0736422128825511|2.7201727841|
|1/2048|41|0.0717401228481560|2.3249164255|
|1/1024|38|0.0699887729305563|2.0045149247|
|1/512|31|0.0687388061787004|1.7990758935|
|1/256|25|0.0684754435288341|1.7581087225|
|1/128|24|0.0697424794959182|1.9625610141|
|1/64|30|0.0728027671675408|2.5393735998|
|1/32|41|0.0775609910090300|3.7127807041|

八次 SLSQP 均报告成功。严格排除依赖完整有理 support/gradient 证书，
不依赖优化状态。共273次 objective/解析gradient评估，没有重启、
改 h、改权族或追加探索。

## 引理、原件与严格标量核验

proof.md 给出 Schur/Cauchy–Binet peeling、选中 Y 与真实原 law
投影独立性、Gamma cut 积分、k1..3 完整32/496/4960计数，以及
k4..6 minor 下界和 Phi 的低 N 零分支。

B0 使用已独审002完整 fan 的 R0/pi^3，其中
R0=520066141/48769479475200。原 producer/audit 保存值及 input hash
一致，未重算峰。B1..3 使用033真全 projection 的准确 R，B4..6 只
对独立子集应用108真原 law 峰，未把 binom(N,4) 当独立子集数。

每个 h 的 Arb512 矩 cap 向外取为有理 beta，最终权重在2^-44
simplex网格上严格核验非负、x0>0、和1，故 EdW<=1。全部33个
Binomial态都进入期望；p=erf(sqrt(h)) 区间、264个倒五次根及8个
M六次根端点经准确幂比较。formal阈值只用123完整有理 Theta，
未用显示小数或MC。

凸性和 Euler 齐次关系给 minH>=(6/5)H+min gradient。gradient逐项用
完整概率区间的准确max界支付；gamma重标将beta网格族映回准确Bcap族。
八份严格 positive residual 全部保存。这不使用实际 determinant 样本。

## 实际运行与保存核验

主 release SHA为
5c4fc638a5eb7dcb3b3ccf888dd4dd36e5228e80f79de8d371490f6ee2b65b74。
其后才以 nice -n10 唯一运行。科学PID89699、observer89690，
tool d1054f/session89797启动，cee94d完成exit0。
完整科学child从spawn前到最终receipt、解释器退出及reap后，
wall1.3698894169938285秒，CPU0.916613秒；observer自身收据/退出
不在该child范围。最多一科学CPU，没有使用121的native资源。

唯一纯Fraction保存checker为tool2387db，exit0，
PASS_SAVED_ONLY_PURE_RATIONAL_SERIES：8h、264个倒五次根核验，
128项交错exp/erf级数、384bit整数平方根及原768bit pi区间独立复核
beta/p，同时核W、33态期望、gradient、gamma、support、阈值及六次根。
0优化、Arb或旧生产程序重放。

复核入口：

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B check_saved.py

科学执行和保存checker均无失败。首次封存文档辅助命令因stdin编码
解析失败，未执行任何语句；该行政错误保存在actual_tool_receipt.json，
随后只修正文档输出，未改科学代码/输入/结果，未重跑科学或checker。

## 范围与停止

19项sourcepins结束时全部未变。全部 h_0..h_7 精确权重、cap、root、
gradient及收据保留；0 MC/native、旧fan/peak/tuple程序、MATLAB、
其他frame、下界算子研究、新LLM代理、stage或push。

全域g0仍未证，本批未通过完整I门，完整S/C6未闭合。一般L1、全部
正共同尺度及同一原generator/source接口不变。未来多h九桶multinomial
枚举未生成或运行，需主另行登记；本轮不扩。
