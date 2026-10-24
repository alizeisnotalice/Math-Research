# 第31轮：12个连续小盒的Bernstein充分LP行

状态：只静态实现，未执行。主侧负责唯一分阶段hard包装，builder与search分别300秒/180秒（合计最多480秒，均含导入和读写）；独立审读后才启动。若builder早停，不启动search。没有参数搜索、细分回退或全域checker加深。

## 冻结数学对象

α=β=47/40、R=33/20、h=1/32；沿用第29轮128方向、相同变量盒与完整51574点池。K=cΨα,R+sS1+Σ126项t_jT_j，c,s非负；30个I+J≤9粗H=1/8方向准确细化、冻结96细方向i+j≤48，T支持≤13/8<R。参数有理化后再Fraction展开，禁止完整a独立舍入。B1仍是已证参照，目标3.4不保证可达。

3个父knot固定(7,5)、(8,6)、(10,6)，各局部端点取[0,1/2]或[1/2,1]的2×2笛卡尔积，正好12个闭小盒，不取邻圈，不临时再细分。物理坐标为h(p+x),h(q+y)。全部严格离轴/原点/截断线，端点触knot由已审弦分支处理，不假设全盒四阶导数有界。

## builder与准确充分行

复用pin的independent `Spec`、`AKChecker.AS_model`、`APsi_taylor`。两个模块导入后设置Arb192bit，重建MC3/MC3T（其导入常数原为128bit）。各方向独立Spec/AKChecker，phi表不跨不同mlist共用、D表绝不跨direction共用。每方向先准确展开到同h后建Spec，不改public代码；Ψ模型未乘c，S1是整体一个模型，保留其内部第四差分取消。

每模型输出4×4 power矩阵P及统一绝对余项E；经MC3*P*MC3T得16个Bernstein球，将其向外dyadic端点准确转为Fraction l,u，保存原球。c、s非负使用l−E上端；126个有符号T使用m=(l+u)/2、rad=(u−l)/2，并引入z_j≥±t_j。对每盒每个ell准确添加：

`c*(lPsi_ell−EPsi)+s*(lS1_ell−ES1)+Σ[mj_ell*tj−(radj_ell+Ej)*zj] >=0`。

16行共同保证整闭盒AK≥0；行不是点值采样，也不是将某方向bern_min乘负参数。共192行、254变量（128原参数+126非负辅助）。辅助可固定上界为相应|t|盒最大值；不改变原参数盒。Arb端点用exact `lower()/upper().man_exp()`转换为dyadic Fraction，不将float或显示字符串当准确行。文字对齐注（运行后的文稿说明）：此处以实际冻结源码的man_exp接口为准；本注仅同步接口措辞，不改变数值代码、已存模型或有理行，也不是新增事前预算承诺。

**已批准早停顺序：**先构全部12盒Ψ+S1，共24模型，计算192条c=c1、s=1、t=z=0的准确基线余量。任何负行则builder终态`baseline_enclosure_failed`，保存`cells.json`为`partial_baseline_only`及`build_result.json`；不构1512个T模型、不删行、不增加盒深度、不启动search。此失败只说明充分包围不接受已证基线，不是基线真实违反或核族不可行。

只有192条基线行全非负才继续126×12=1512个T模型，总最多1536项（12 APsi+1524 AS），保存所有power/Bern球、准确端点/余项、盒/方向来源、Spec/hash和完整有理192行。每模型入口先记started/in-flight，模型与有理端点落盘后才记completed；导入、精确细化、Spec构造及转换/保存均计300秒。首次异常/非finite/超时就保存终态，不重试。`cells.json`只有完整+baseline_passed才能供search使用。

## 通过builder后唯一search（先等待builder，不提前运行）

读取旧hash缓存和新完整准确192行，以相同51574×128点矩阵及ε1/2旧基线slack构建254变量模型，加192盒行及252条z≥±t线性行；没有新点、闭式B/L或方向。全部准确行保留，MATLAB用float求解但不冒称准确模型。独占一个batch仅一次linprog20秒soft，总search180秒hard，工具错/非optimal/无finite primal终止，无fallback或第二LP。

浮点解与唯一结构参数有理化共最多2批103148个完整池K/AK值；原/加强/raw门槛沿第29轮固定K1e−9/A1e−8，完整系数/成本交叉不省略。有理化c,s,t_j后，准确取z_j=|t_j|，逐192条Fraction安全行**无容差**核验≥0。若失败不修补；只导出至多1个精确结构candidate。准确成本由主侧独立算，严格低于B1后才安排既定全域认证；12盒仅局部证据，不能替代其他区域、轴/原点/截断线和无穷远覆盖。

文件ownership：本代理PLAN.md/build_cells.py/search.py/batch_solver.m；主侧run_stage.py/enclose_cost.py与归档；audit独立审读。builder输出cells.json/build_result.json，search输出result.json/candidate.json，旧/未完整证据拒绝覆盖或续跑。没有运行前不可将计划预算写成实际次数。
