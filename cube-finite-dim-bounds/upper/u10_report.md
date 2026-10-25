# U10报告：联合空间修正的有限连续认证接口

2026-10-25，base fcc2a7dbdc0c26ca28352727b8201c6b701d11c0。只upper/u10_*，不改旧证据/bridge/总账，不提交推送扩招。完成一个联合充分接口即停止，无主LP/µ求积/参数搜索。没有新U3<8，U4不在本批具体参数范围。

## 完整机制与可核验输出

固定n3 frame V={(1,±1,±1)}，真实L=−Σ|v·D|，c=(4log2−1)/(24pi³)，H=(G−c)+、LH=μ−δ0来自026已审U09。新基为b0=(1−|x|²)+²以及六个center±3e_i/4、radius1/2同形bump之和b1。w0,w1≥0，F=H−w0b0−w1b1，K=F/D，D>0。

完整证明u10_joint_certificate.md已锁SHA7d3ce12b4cafc3db7d07f3c8c9359f165a777d35f462a48d048f4b46eefae62e：

- C1,1基的Lb是普通连续有界L1密度，无球面atom，∫Lb=0；L1尺度为R^-1，六小球总质量只付一次。
- 两支持均在r≤5/4，六小球互不相交，b1≤1。Jensen真实G≥5/(16pi³r²)使closed球r≤3/2位于Ω，给support到Ω^c的显式buffer1/4。
- 全域非负只需support内核F≥0；全cube另须F≥D/8，D不得凭未证g=edge设定。
- Ω外源确可由F≥0自动结算：{G<c}为openzero区，Markov跳入为正；{G=c}零Lebesgue，故减去Lb dx不扣该处正测度。此处用support buffer及level零体积，不是因b紧支撑就忽略外尾。
- Ω内必须检查整个signed联合源m−w0Lb0−w1Lb1≥0，不能逐基只付正部。m有非负连续密度(1/pi)Σ∫(c−G(x+sv))+/s²ds，适用于Ω内line点；phase未定cell不能套该式。

## 无穷域的显式有限归约

令W=w0+w1、a=5/4，S为四frame lines，T_a为这些line的radius a投影tube之并。tube外Lb均0。two-large-jump与真实heat bound的准确最小积分给

G(x)≤48√6/[pi³dist(x,S)²]。

对far tube x=t v+z、z⊥v、|z|≤a，选另一frame w，在s∈[t,2t]的目标点远离S；t≥60时G≤c/2，故m≥c/(4pi|t|)。同时LB≤40W/(pi√3 t²)。因此

T=max(60,160W/(√3c))，R_*=√3T+5/4

以外正源自动成立。源认证只剩有界Ω∩T_a∩B_R*。这是解析尾支付，没有靠固定源点外推。

近line亦给显式屏障：x=t v+z，|t|≤A、ρ=|z|<1时

G(x)≥log(1/ρ)/[200pi³(1+(A+2)²)]。

可由此认证line邻域的cube/非负下界，以及h=min(G,c)恒c的small-jump零区；µ积分剩余大参数尾有[0,8c/(pi S)]包络。常数粗但有限，未声称实际证书高效或一定存在。bridge可接更锐区间；没有强迫G∞进入浮点求积。

## 准确质量与全L1/尺度

M0=32pi/105，M1=8pi/35=(3/4)M0。满足全部条件后

C(K)=[∫H−w0M0−w1M1]/D。

此时新正源ν=μ−LBdx为概率，ν0=0，LK=(ν−δ0)/D；K满足U09的全部核输入，故一般f∈L1、所有共同正r使用同一bad set的弱型结论成立。不涉及独立矩形多尺度。恒定Lebesgue实际用于密度、普通卷积/Fubini、球质量、level零体积、源总质量、cube分母和接收体积；不推广非恒定背景。

## 实际运行与证明状态

只运行已登记u10_check.py一次，SHA7df5b62eab7fe2617262f7254a5374b8c424825267999e701f237159c47c9918。SymPy准确核：球间距平方9/8及9/4、半径5/4和buffer1/4、两基质量及缩放、四triple det4、line basis剩余坐标(-1,1,1)、横向segment的距离多项式、G尾常数与source支付因式。全部PASS。

唯一合成算术w0=w1=1/1000只测试闭式，给T≈77.56079088、R≈135.58923048，**不测试或主张该权重可行**。无source/G数值积分、无新grid/LP或优化。

MATLAB脚本u10_matlab.m SHA98e81147594ef5cae6347ca030fa0f356ba87129a706253be8d19fc2e4ee8e9e，upU10前缀，无clear/path，已交主串行；主已串行MATLAB实际exit0，script hash匹配，支持/质量/几何均一致，mass defect2.23e-16、G bound代数defect8.89e-16；synthetic tail slack4.24e-22，未称候选。完整receipt及log路径/hash已入u10_results.json，未改共享会话。

|类别|本批状态|
|---|---|
|严格分析/准确算术|联合充分接口、Ω外自动结算、有限tail和line屏障完整锁稿；准确支持/质量/常数PASS。独审排队，尚未登记通过|
|实际数值|一次固定算术Python PASS；主串行MATLAB实际exit0，最大代数差8.89e-16|
|未认证|任何实际w0,w1,D的全support/全cube/有限Ω源区间覆盖；主样本或warnings不补这些条件|
|未实现|无U3<8，无两个基或更广空间修正的普适障碍；未给n4具体tail接口|

主在其分支的2basis/8basis有限LP和留出测试属于独立数值发现，本分支没有重跑，也不把它们当通过条件。下一步只能为具体权重补连续认证或另登记结构；本有限理论批停止。

## 白名单

仅8文件：u10_inputs.json、u10_preregistration.md、u10_joint_certificate.md、u10_check.py、u10_results.json、u10_matlab.m、本报告、u10_manifest.json。输入来源/哈希在inputs，最终源码/结果hash在manifest；不自记manifest自身hash。026已审旧源未改；新接口须另审，不借026自动登记本轮通过。
