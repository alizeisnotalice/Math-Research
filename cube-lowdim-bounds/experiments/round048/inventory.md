# 二维下界重启：只读盘点

2026-10-24。已停止上界round29；本次只读取文件、日志和ZIP中央目录，没有运行数学核验、优化、编译或Git。以下历史结果不算本次重新核验。

## 现有证据及缺件

|对象|准确依据/位置|状态|
|---|---|---|
|repo二维Lin邻近有理构造|`results/lin_rational_certificate.json`，ratio=`167658289715427138435952689653/99470362309541620000000000000`，显示1.6855099933555249|已有准确证书；不是最新公开基准|
|公开Tomar50原子|`tomar-audit/upstream-v4/lowerbound/c2lower_1p70272.json`，SHA `e2ca1616be95afa6680e3ccf157cc16fdd39026ad0f5fa79ebd32e5214779888`|完整有理配置实际存在；本次未运行核验|
|新md所称旧内部50原子|`16979690190120226811729625299082623770069528937/9971651045727228079122845811332874240000000000`，约1.702796268366800390|只有说明，未定位原JSON/准确运行日志|
|新md所称58原子|`33977370204434849545366297954842256187719839393/19950077409778055009062309259909267456000000000`，约1.7031197176097999；md给candidate SHA `11d4e18fc9688580108ad4e2c19e6ecf0a5f346122e7cd32f573a4a2159b2879`|candidate_final.json、verified_result.json及原始两准确日志未到手，未确认/未接纳|
|三维R3|`results/baseline_algebraic_certificates.json` dimension=3，准确多项式及隔离区间；历史日志给1.79393322779477<R3<1.79393322779478|既有严格结果，复用无需重跑|
|四维R4|同文件dimension=4，1.90174261110095<R4<1.90174261110096|既有严格结果，复用无需重跑|

repo根：`/Users/zhengzhihao/Desktop/T/Math-Research-lowdim-agent/cube-lowdim-bounds`。scratch根：`/Users/zhengzhihao/Desktop/T/cube-lowdim-continuation`。

repo账本仍记`project_lower=Phi`；公开下界字段`17027237335/10000000000`是展示截断，且标明未在round20复跑，不能替代公开准确证书。公开`NOTE_lower_bound.md`给准确值为
`482167717377842991882283179477574801392641890014291050323/283174368139841919615520431313000857449005056000000000000`，约1.7027237335962937。本次只读原说明，不声称已重算该分数。

在本次定位的有效低维材料中未发现强于R3/R4的完整证书；不是全磁盘或全部内部项目的穷尽结论。`finite_certificates.json`的n=3/4/5是**二维原子数**，不是空间维数；各比值约1.21，均弱于Lin。代数证书dimension=2项`3/4+3sqrt(6)/8`约1.66856，也不是Lin Phi。

## 可直接复用的准确文件和实际日志

- `results/baseline_algebraic_certificates.json` SHA `2c148ff95f9dcb8f5094c9fb2fa96e248d1f07810518789cbfb725c3c52c5ab2`：含R3/R4全部根定义、参数、多项式系数、见证counts和有理隔离区间。
- `results/baseline_certificates.json` SHA `d41105e016387f33ece5c221738fa60dbf4c5a2e80d4a9687b6bdb991dda65fa`：稍弱的三四维有理坐标证书，准确比值分别为`112120814808410326284319361162979567/62500000000000000000000000000000000`及`1901742484796073585488758280083913395469464402651/1000000000000000000000000000000000000000000000000`。
- `proofs/lower_transfer.md`：有限共同半径见证、完整质量、周期块截断/边界损失、原子平滑到一般L1。坐标闭包S→T极值不变、质量增加，不能把同补丁子集改名当新机制。
- 本地`/Users/zhengzhihao/Desktop/T/cube-lowdim-stage-20261023/logs/`内`baseline_integer.log`、`rerun_verify_bounds.log`明确打印上述有理比值及threshold/independent union PASS；`baseline_algebraic.log`、`rerun_algebraic_bounds.log`打印R3/R4准确区间和多项式；`audit.log`记录结构匹配及105组整数并集测试all_passed=true。这些历史日志已实际存在，本次未重跑。
- repo第6轮每维200权重、第11轮每维80位置实验均未改善R3/R4。scratch `round006/search.log`、`verification.log`和repo相应JSON可复用；独立体积算法共享盒生成，不能夸称全流程独立。

## 现有工具

|位置|能力与注意点|
|---|---|
|repo `baseline/verify_bounds.py`、`algebraic_bounds.py`|标准库准确整数/代数证书；不是通用任意原子优化器|
|repo `baseline/uniform.py`、`weighted2.py`、`search.py`、`src/search_stage.py`|均匀格子/parity权重/自由有限原子数值发现|
|repo `experiments/round005/009/015`|错位双行、卫星、h/V/w联合；已完成有限失败记录，不仅换seed重复|
|repo `experiments/round006/011`|三四维独立权重/非乘积位置扰动，保持共半径和完整周期质量|
|repo `experiments/round019`|四原子超周期，642目标+2独立几何核对完成，无改进；scratch round020九变量准备已停止未运行|
|scratch `tomar-audit/upstream-v4/lowerbound/check_exact2.py`、`independent/lb_check2.py`|公开两套准确矩形枚举/并集入口，完整50原子配置键为pts；本次未运行|
|同目录 `search/eval.c`、`hleval.py`、`search.py/search2.py`、`seg_search.py`、`sym_hf.py`、有理化脚本|C浮点评价、周期质量、真实union、自由/反射/线段族；hleval导入需libeval.so，本次未导入/编译|

新md列举的34参数续接工具`search_followup.py/probe_orbit.py/add_shared_orbit.py/joint58.py`、best JSON及`audit_saved_exact.py`尚未定位。不能用公开等间距线段工具冒充取得了新工具/准确构造。

## 实际查找范围与下一步

新两md在WeChat本月同一临时接收子目录，名为`Pro_lower_and_finite_dims_20261024(1).md`和`cube2_lower_construction_20261024(1).md`。后者声称58原子两次核验及20989/1685/236/10000搜索调用、三四维体心16组失败；缺原文件时这些只属说明声明。

1. 用rg文件清单在Desktop、Downloads、Documents及指定WeChat RWTemp/2026-10查candidate_final.json/verified_result.json精确文件名，未命中实际文件。
2. 在Desktop/T的Math-Research副本、output及该WeChat月目录对新准确数值、旧内部完整分子和candidate_final.json作固定字符串检索，有效命中只有两md，未定位构造或原日志。
3. 对Desktop/T、Downloads、Documents及该WeChat月目录的665个ZIP中央目录完成只读检查（0错误、约1.37秒），查candidate_final.json/verified_result.json/best_joint58.json/best_continued.json/check_exact2.py/joint58.py，未命中。没有解压/执行，也没有继续扩大搜索。

范围不包括全用户盘、隐藏/忽略文件、其他WeChat历史月份、7z/tar或ZIP内嵌ZIP，不能写“全磁盘不存在”。主侧另行网络精确检索/远端树检查未定位，是主侧证据，不冒充本盘点的网络核验。

当前先登记缺件：需取得SHA为11d4e18f…的58原子有理JSON、verified_result、两原始准确日志及搜索源码；旧内部50原子也需原件。不得从展示小数重建。等待期间可由主侧另设预算核验实际存在的公开50原子；既有R3/R4直接承接。本盘点没有新数学结论。
