# 第42轮：node663 的四符号区域准确代数证明试验

状态：静态准备，未运行；源码审读后由主侧唯一60秒硬包装运行。仅使用标准库Fraction。0新Arb、Hint、积分、AS、Psi数值模型、点值、T或LP，不修改旧树。

固定node663来源模型SHA `8b614d16e19e8dddd9cb16793595b99aff3acf680795710c50ced9c3f76c108e`、round22 cover SHA `e78b2bf4e6f5e00974a180cb85dfa9044ed42939930c4add48078e6c4064c8e7`及B1 SHA `21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172`。旧physical盒为[6677/4096,3339/2048]×[161/4096,81/2048]，旧tree depth4。原T1/AS的16项准确Bernstein端点、ES及Muv上端全部读取保存值，不能重新评价或取代。

## 下界与四个证明区域

已审F=AΨ满足非负对角曲率与0≤Fuv≤Muv，所以

`F ≥ T1 + min(0,Muv·δu·δv)`。

按原物理中心(uc,vc)分成四个闭矩形：(--),(+-),(-+),(++)。相邻区域共享边界，内部不重叠，准确并集是旧node663。它们是**新代数proof pieces**，不是旧失败流程偷增tree depth5；原node663/depth4/失败记录不回写。

先对已存T1与AS Bernstein**下端数组**分别做准确de Casteljau二等分限制（先u后v）。这些限制是非负凸组合，仍为对应真实多项式的下系数。然后同号区域加0，异号区域加准确双线性多项式Muv·(u−uc)(v−vc)；先构造该子盒2×2端点Bernstein数组，再以准确degree elevation升为4×4双次数(3,3)。不直接把全盒常数D继续扣在每项上。

AS只减原ES，不新增误差项。每区域16条 `c_column=restricted_T1_lower+bilinear_coefficient`、`s_column=restricted_AS_lower−ES`，基线slack=c1*c_column+s_column。c/s列用于以后非负c、s参数；本次固定基线c1>0、s=1。全64项严格正才accepted，不把负充分下界当作真实AK负点。

## 有限计算与比较

总60秒hard包括读源核SHA、准确限制、升阶、64行与保存；内部各步骤和行循环查时间。恰4个预定proof pieces、最多64条新资格行，做完后保存通过或未通过，不新增第五区域、不自适应细分、不调Muv或ES。来源旧16行只作保存值输入，旧行的准确限制结果在同一64行循环中核对，不另评价新函数。

保留旧统一D=Muv·du·dv作对照：准确限制旧c/s/slack数组，逐项验证限制等式、`new_slack−restricted_old_slack=c1*(D+bilinear_coefficient)≥0`。这既比较符号模型增益，也检查基底/索引方向，不凭浮点近似判定。

输出quadrants.json：旧树身份/来源SHA、中心、每piece准确local/physical几何、sign标签、proof_piece标志、限制后T1/AS数组、双线性提升数组、64条新旧c/s/slack/差；result.json保存真实64行计数、终态、accepted、准确最小值、时间及产物SHA。已有证据拒绝重跑。主侧/audit用独立power substitution或Bernstein restriction公式核查，不能把本脚本自检称独立验证。

成功仅补证固定旧盒的四区域单侧充分条件，不代表整个父格及其216pending闭合，不修改C2或开始T方向。任何后续树续接须另定范围和预算。
