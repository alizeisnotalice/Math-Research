# 第135轮：删源几何能量与反向尺度软尾比较

一般n^(1/2+o(1))一致弱型上界仍未闭合。

- SC1–4：任意相关来源及任意非乘积删源下，真实来源中心Lebesgue积分的加权score方差为O(n)，对应正score为根阶。独审通过；未将e_t(y)dy偷换为实际β(dy)或接收dx，跨测度费用仍未支付。
- JG135.10–18：以真实硬支持坐标数控制K_L、K_L+对纯软核J的比值，给覆盖Hhat<L+的新RC联合阈值。在明确计数层CJ≤3，保留K0本身及bad来源余额；独审通过。原捕获容量不重复领取，重尾总费用未付。
- 真实核数值改变研究重点：12例中10个基础必要门正余量行全在Hhat<L+，另外2例旧RG16成立但无good来源。加回真实ρδ后仅7行仍有正的原始时间交通，另3行score饱和给零；FA/Ψ/history仍未认证。有限菜单不认证FULL或一般维数率。

## 验证

真实核80/100位12例；主代理120位补核残量、tail及来源核字段，相对差小于10^-65、门一致。新较弱反向尺度上界在同菜单60来源点检查通过，其中30个positive点全在CJ≤3软尾层；仍有未消除余量。新增K0项的更强版本有解析独审，此数值脚本未单独核它。另有12例Fraction有限平移群代数核验，以及12例含非零硬坐标、闭边界与s=1的核比较浮点压力检验。数值独审含少样本独立积分和分支核对，无区间认证；具体范围见归档。

## 归档与复现

[research.zip](research.zip)：187,209字节，16成员，SHA-256 `cf8380b97c098d7692a8d82b154e6697dee6a41bf0ba66565b1642b8521eab76`。含详细证明、独审、LaTeX文本、冻结菜单、必要脚本与完整结果。CRC、资产散列及82处显示数学/公式标签保留检查通过。未编译PDF；只新增本README和ZIP，并更新索引。

使用Python与mpmath==1.3.0。便携复现用replay_real_joint.py与冻结real_joint_numeric.py/JSON；新目录无同名复算输出时执行python3 -I replay_real_joint.py，再执行reverse_scale_numeric.py、score_gate_probe.py及kernel_count_probe.py。有限群channel_numeric.py只需标准库，亦须在无输出JSON目录运行。冻结producer的main读取原本机历史文件作散列，不应作为便携入口；原run_count不包含后来的复算。

剩余主瓶颈：在实际same-x门和真实dx下支付纯软尾的κ_r与A/wg联合交通，并处理计数层外显式维数费用，使总价为εX+C√n polylog(n)Wβ。β−ν、尺度/分层、history及未消费负项全部保留；不以√nX代替目标。
