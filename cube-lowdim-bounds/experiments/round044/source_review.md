# 第44轮builder运行前静态审查

最终build SHA `a03f2f13909d95b140505c2c870cbfe901c2755e93d3c1666de715a603501f8f`，结论GO，无静态阻断。只读审查，未运行模型/数学。

- pin固定B1、round24完整cover及独立通过报告、12方向/old_column顺序、旧rows与12204旧模型。615/661/663及proof_file来源按实际hash核验；历史partial不要求存在。
- 591树叶对应1278region，whole直接用、四片以原node/piece标识全部展开；准确c/s行沿用20448项，不按16L误计。失败分裂父不进入接受cover。
- 每方向独立Spec/AKChecker，使用已pin方向coefficients（coarse已准确refine至h=1/32）；D4半质量、参数盒和T支撑≤3/4核对。设192bit后重建MC3/MC3T，局部坐标与父(52,1)保持。
- 7092次新AS，每个整叶P/E保存后才限制。lower/upper分别准确de Casteljau；整叶E不变；每region每系数m=(l+u)/2，penalty=(u−l)/2+E，供tj*m−zj*penalty，不丢signed误差。
- 固定15336方向区域记录及20448新rows。每dir流式保存、source model行SHA与文件SHA形成来源链；12方向行流按region/index核对。solve读取schema与实际directions_file/legacy_rows_file/new_rows_file及hash一致。
- 300秒外层hard、7092新AS上限，状态/在途键和partial保存，无重试/新Psi/const/独立点值/LP；未完整不得启动solve。末尾bytes/RSS补丁仅记录真实文件大小和平台ru_maxrss单位，已局部复读。
- 独立保存值核验另限90秒：只核7092模型→15336region→20448行，用幂基推导的限制矩阵，与生成器de Casteljau不同；旧16272仅pin复用。不重算Arb/AS或追加模型。
- 审查边界：使用已核解析引理/已有区间原语；静态GO不等于模型数值通过、全域新核认证或C2<3.4。
