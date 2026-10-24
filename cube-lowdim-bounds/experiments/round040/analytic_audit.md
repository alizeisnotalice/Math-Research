# node661单侧混合曲率扣项：独立解析审查

结论：所提`Dpsi=Muv*du*dv`为有效的**单侧下界扣项**，旧T1 Bernstein可直接复用；无解析阻断。本审查没有积分、模型、点或符号求值。

令F=AΨ，T1为同一盒中心的一阶Taylor多项式。已审外区积分导数给Fuu、Fvv、Fuv均非负；Fuv=Hab(u,v)+Hab(v,u)，近角点的两个Hab上界之和Muv控制整个盒的Fuv。Hab公式包含正确的两倍归一化，修正积分传播analytic标志。

对中心到任意盒内点的增量δ=(δu,δv)，积分Taylor余项为

`F−T1 = integral_0^1 (1−t)[Fuu*δu²+2*Fuv*δu*δv+Fvv*δv²] dt`。

两个对角项非负。若δuδv>=0，整个被积式非负，故F>=T1。若δuδv<0，利用Fuv<=Muv得

`F−T1 >= 2*Muv*δu*δv*integral_0^1(1−t)dt = −Muv*|δuδv|`。

因|δu|<=du、|δv|<=dv，全盒F>=T1−Muv*du*dv。交叉项中的2与积分的1/2抵消，最终扣项不再除2。此推导只用逐项符号，不产生`|F−T1|<=Dpsi`的双侧误差结论。

两个Hab球可取绝对值的准确向外上端再相加，仍为安全Muv，即使数值球下端未能严格正；真正的非负性来自解析导数。node661整盒正坐标、r>R、u1,v1<R，无跨坐标分支。将来跨R时，已审C1零延拓、连续梯度和分片有界Hessian足以给相同线段积分论证；本轮不实测该扩展。

已静态追踪第39轮保存链：`_outside_taylor`返回中心T1的已升阶4×4 Bernstein数组和独立余项；`psi_record`保存原数组端点，直到构造rows时才扣旧余项。因此旧Psi Bernstein**未扣E且未乘c**。固定旧模型hash、几何、修正源码及基线身份后，直接复用该T1包围合法，不需重算F0/梯度。旧AS准确端点和ES亦保持原样。

新的16条准确下界为`c1*(oldPsiBern_lower-Dpsi)+oldASBern_lower-ES`。c1>0，Bernstein基非负且和为1，全部strict>0就给这个闭盒的AK1>0。旧AS误差仍须扣除。新扣项是否比旧绝对余项更小由实际包围决定，不能先验把数值改善当接受标准；最终只按准确新行判定。

本轮最多2次近角Hab/2次实际积分、192bit、60秒，0新AS/完整Psi/const/F0梯度/点/细分/T/LP/全域checker。成功仅表示固定node661准入，失败亦不提供真实AK负点；不外推父格完整覆盖或新C2上界。
