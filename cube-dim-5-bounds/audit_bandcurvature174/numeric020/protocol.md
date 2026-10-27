# 174：实际173v3的唯一新有限数值审计准备

固定实际输出manifest `28568080f4779a79b134c45b9e707ceb9d9bb2a78b64f17cff1b0a587995ce29`；v3准备SHA `cceaa29e847719bc90cbbc1ef7005f8cf2c62c1924772d942003740233b773eb`。复用已通过的455来源门、原151 jet/projection与168/旧71pair metadata资格；不重跑455来源审核、旧H/log/梯度/旧q/旧角点。prepare只读SHA、AST和已有projection标识，不调用数学函数。

待主授权的新算术恰为455个q、455×49=22295个正band乘积（含完整无穷尾）及实际34170组×16=546720角点；全部35872事件必须审核，包括1702 firstPASS skip及旧71pair严格q降guard。若实际schema/count或任一准确等式不合，不增加算术预算，保存prefix并停。

新math174.py与已读future_math.py的数学路径相同：物理band[a,b]换z=1/t，写B=1/b（无穷尾B=0）、d=1/a−B，用正二项式积分P_n=Σbinom(n,k)B^(n−k)d^k/(k+1)，权P5/P3。不同于173差幂公式。source_q只新增保存已经形成的m*w值供逐项比较，避免第二次形成band乘积；并将0<gL≤gU、所有端点正有序、同cutoff交叠及0<L32<original_gU写为强制前置条件。K正上包的分母为16*max(gL,dL)，153负项仍除original_gU。min只在原qualifiedq、band、512/3、153 split之间取；不放宽原来源(point/epsilon/gamma)资格。

角点先把anchor、Delta、四gamma区间、q及盒端点清到同一正公分母D；16个整数余量为A D²+DΣmin(BL X,BU X)−QΣX²，最后除D³。每组的准确最小值、lex平局bits和passed状态必须与生产相同。同一anchor全16角点；允许盒外anchor及负Delta。事件去重、firstPASS顺序、旧pair guard与失败pair改善best均重建；使用已独立核准的新q，角点阶段不再求q。

15个显式worker分两阶段运行；audit shard=(producer shard的逆7置换)，不产生新点。每阶段先真实pilot1，成功后自动推进所有固定输入。每个成功case fsync ACK，source每个case和corner至多128事件保存checkpoint。成功prefix复用ACK/hash及metadata回放，不重算q/角点；手动续跑还必须有主授权的准确prefix SHA。worker异常立即stop，主不等齐15份summary，进入保存故障后的5秒grace、2秒terminate、2秒kill及实际join。正常没有wallcap，不自动重启。

数值完成后从原155/已接纳158 margin字段与独立1120 receiver best重建完整111201叶的新globalmin，不重算任何旧corner。逐根核exact Kraft，并另做自底向上二叉兄弟合并闭包，完整52360根；完整candidate、updates、pending均逐字段与生产比较。base158和全部实际测试的失败新source改善都必须保留。最终代码/source/protected SHA全过才写PASS结果。

原H73端点、H64及四物理导数区间仍共享017/DD/Arb及080/094/102/113/121、151/155/157/180既有真值资格；167是transport而非新H enclosure。新独审若通过，只认证此次新q/角点/完整消费，不自行更新I、g、C5或替代主的全尺度证明责任。

运行必须精确传`--plan-sha`、主授权文件及其SHA，approval字符串为`audit_bandcurvature174-numeric020-approved-after-review`。主授权JSON字段须含plan_sha256、sources=455、band_products=22295、pair_events=35872、actual_corner_cases=34170、workers=15。准备阶段不运行此CLI。
