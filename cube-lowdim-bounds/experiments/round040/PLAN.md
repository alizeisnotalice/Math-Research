# 第40轮：node661 的有符号单侧混合曲率扣项

状态：静态准备，未运行；须经独立数学/来源与源码审读后，由主侧唯一60秒硬包装执行。不改第39轮结果、模型或源，不重试原父格。

## 冻结对象与复用合法性

唯一node661，physical=(6677/4096,3339/2048,5/128,161/4096)，local=(21/128,11/64,1/4,33/128)，depth4。整盒r>R，两坐标均严格<R，当前不使用跨坐标R分支。固定B1 α=47/40、R=33/20、h=1/32、c1不变。

输入pin：round20/leaf_cover.json SHA `6028e48efe444e27c5259858cf1eeb8b41e0081b6e39f06102439ddb198c4cac`；node661最终模型 SHA `8dd22932cefd30806acf819ffae37f10d9497e50e2e798361e10f0d9e207db26`；B1 SHA `21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172`；修正first hlcert SHA `4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05`。核cover节点、模型meta和准确几何一致，缺失则停止而不重算。

固定源码 `_outside_taylor` 先构造中心一阶T1的2×2 Bernstein，再升至4×4；保存node661的Psi Bernstein正是这个T1，未乘c，不含余项。故可直接复用其16项准确端点，不需重算F0、梯度、完整Psi模型。旧AS16项和ES也直接复用，不能从显示幂系数反求。

## 新调用及单侧证明

已审外区公式给F=AΨ的Fuu,Fvv,Fuv≥0，且

`Fuv = Hab(u,v)+Hab(v,u)`。

近角点单调性控制全盒Fuv。仅在近角点(u0,v0)和交换(v0,u0)各调用一次修正 `_Hint(...,'Hab')`，以两球绝对值准确上端之和作为Muv。只为稳健包围取绝对值上端，不依赖数值球下端严格正。

沿线段积分Taylor，正对角项不扣：δuδv≥0时F≥T1；异号时F≥T1−Muv|δuδv|。所以全盒有单侧下界

`F ≥ T1 − Dpsi`，`Dpsi=Muv*du*dv`。

这里Dpsi是**有符号曲率推导的单侧下界扣项**，不是双侧绝对误差。本盒同分支；未来跨R使用还须承接已审C1零延拓/分片Hessian引理，不由本次结果额外验证跨R数值行为。

## 固定预算及输出

192 bit，模块导入后重设再构造对象；唯一60秒hard含核SHA、构造、两积分、16行准确算术和写文件。最多2次Hint、2次实际acb.integral；两坐标均<R，预计均实际积分，仍按真实账本记录。代理原样转交acb构造、原callback/参数及Hint，不改公式/积分选项；首非有限、超时或异常保存停止，不重试。

0新AS、0完整Psi模型、0const、0F0/梯度/独立点值、0细分、0T、0LP、0全域checker。原ES、原Psi绝对余项、T1 Bernstein与旧16行保留作准确对照。

一次16行组合：`new_c_column=oldPsiBern_lower−Dpsi`，`s_column=oldASBern_lower−ES`，`new_slack=c1*new_c_column+s_column`。保存旧slack、新slack、差、原绝对扣项与新单侧扣项；全部严格正才accepted。zero/negative充分下界不当作真实AK反例。

输出result.json（终态、真实Hint/积分账本、时间、准确最小值与SHA）、signed_model.json（两Hint球准确端点、Muv上端、Dpsi、几何与来源、16条新旧行）、calls.jsonl（逐调用事件）。已有证据则拒绝重复运行。

成功仅是固定node661充分行闭合，不是父格完整覆盖、T方向准入或新C2上界。若失败只保存停止，不调网格或精度；下一阶段范围另行确定。
