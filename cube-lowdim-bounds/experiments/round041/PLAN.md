# 第41轮：续接外区树，采用有符号单侧 Psi 下界

状态：仅静态准备；源码审读后由主侧唯一运行。不得现在执行，不修改第39/40轮来源或重算旧模型。本阶段不求T方向或LP。

## 继承状态与来源

继承round20完整树（cover SHA `6028e48efe444e27c5259858cf1eeb8b41e0081b6e39f06102439ddb198c4cac`）：360个接受叶、218个pending、唯一深度4失败node661。以round21/signed_model.json（SHA `d8e9e38f0c7d61c9b49ea7d5875ab823c85c02c58ba6ef242c8315471f7c0cdd`）的准确16行补证661；其旧AS/T1来源是round20 node661模型（SHA `8dd22932cefd30806acf819ffae37f10d9497e50e2e798361e10f0d9e207db26`）。初始361接受叶，原218pending栈次序和全部node/root/parent ID、端点及原根深度均保留。

node615已经有round18 AS与round19 Psi双来源，继续原样引用、核SHA；node661再增加round21 signed补充路径/SHA与previous_status，不回写来源失败原件。旧最终模型强制核SHA；历史partial路径/原SHA仅作可选元数据，缺失不触发重算。旧446个AS、第37轮全部const、第38轮node615 Psi、第39轮61Psi以及第40轮2Hint均不重算。

固定B1 SHA `21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172`，α=47/40、R=33/20、h=1/32、N=80、c1>0，仅父格(52,1)，所有盒正坐标且u0+v0>R。跨坐标R适用性使用已审H零延拓C1、分片有界Hessian及单调性引理，不假设全盒C2。

## 新盒固定模型

每个新唯一盒，只一次independent AS_model；导入全部模块后重设192bit、重建MC3/MC3T，才构造两套对象。AS保存4×4 Bernstein准确端点和非负余项ES。

新Psi bundle不调用旧完整绝对余项模型，只用修正first `_Hint`（固定hlcert SHA `4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05`）：

1. 中心(uc,vc)处共6Hint：F0=H(uc,vc)+H(vc,uc)，Fu=Ha(uc,vc)+Hb(vc,uc)，Fv=Hb(uc,vc)+Ha(vc,uc)。构造一阶T1的2×2 Bernstein并升至4×4双次数(3,3)，未乘c。
2. 近角(u0,v0)与交换处各一次Hab，保存两球准确端点及绝对上端；Muv为两个绝对上端之和。
3. Dψ=Muv*du*dv为准确有理的**单侧下界扣项**。由Fuu,Fvv,Fuv≥0及线段积分Taylor，有F≥T1−Dψ；不能称其双侧绝对误差。

新16行 `c_column=T1Bern_lower−Dψ`、`s_column=ASBern_lower−ES`、`baseline_slack=c1*c_column+s_column`。准确Arb lower/upper.man_exp转Fraction，全16项严格>0才接受。否则原depth<4则中点四分；首个depth4仍不接纳即存证并终止，不继续其余pending。0新const、独立点值、T、LP、全域checker；无fallback或旧模型重算。

## 预算、账本与保存

300秒外部hard，包含读取/核SHA/来源树、导入构造、计算、端点转换及持久化。最多1500新started盒、1500新AS和1500 signed-Psi bundle；每bundle至多8Hint与8真实acb.integral，全阶段各≤12000。Hint的L≤0分支可以不积分，分别记录实际次数，不把Hint数等同积分完成数。

保留成熟append-only calls.jsonl；result只保存计数、在途栈/最后事件，不重复序列化增长调用数组。AS、Psi、Hint、真实积分各有started/completed/failed/interrupted，已返回有限数值与已保存准确模型分列。代理仅原样转交回调/参数；不改公共源码或flint全局类。调用前后、盒循环与逐项转换查时钟，外wrapper负责硬终止。

每个新盒完成AS后保存partial checkpoint；完成Psi后保存新单侧模型checkpoint，最终模型不可变。首次非有限、负ES、非法单侧扣项、异常、深度/盒/时间上限均保存并停止，原pending与inflight完整保留，不加深/升精度/换模型或重试。

## 输出与接口

leaf_cover.json继承全部旧树，仅rebasing来源路径及增加补证/新节点元数据。旧const/旧poly接受叶保留旧行；node615 mode=outside_poly_reused，保留round18模型加supplemental_Psi_file/hash指向round19；node661 mode=outside_signed_reused，保留round20模型加supplemental_signed_file/hash指向round21。新叶mode=outside_signed。

新final模型含AS={bernstein:[{index,lower,upper}×16],remainder_upper,power_balls}；Psi={bernstein:[同式×16],one_sided_lower_deduction,Muv_upper,du,dv,mixed_derivatives,center_values,one_sided_only:true,unscaled_by_c:true,Bernstein_degree:[3,3]}，及rows16、来源SHA与准确几何。Dψ不使用remainder_upper名称，避免与旧绝对余项混淆。所有来源引用可跨旧目录，原件不修改。

只有0pending/unresolved/inflight、全部接受16行严格正、准确面积和1，exporter才标geometry_complete；主侧另核完整分区/树/来源，不能用总面积代替无缝分区。继承历史调用与本轮新调用分开。成功仅为父格(52,1)取得完整连续c/s安全行，不构成新C2上界；下一T模型/LP另立预算。
