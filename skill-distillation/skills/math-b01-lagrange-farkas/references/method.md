# Lagrange 对偶与 Farkas 证书：方法与验收边界

1. 最小化g≤0取λ≥0，L=f+Σλg+Σνh；先证明对每个原可行x及对偶可行乘子q(λ,ν)=inf_x L≤f(x)。

2. 有限实矩阵Farkas形式明确为：b=Ax,x≥0，或存在y使Aᵀy≥0且bᵀy<0，二者恰一成立；零矩阵另直接处理。用精确运算核验方向及严格负margin。

3. P-a9029e19a7586865用有限生成锥闭性及最近点构造y；证书conditioning依赖dist(b,AR₊ⁿ)，不提供统一正margin。

4. 混合LP来源P-0f4471d96cc82d4f采用A_Ix≥b_I、Aᵀλ=c、λ_I≥0：精确gap=cᵀx−λᵀb=λᵀ(Ax−b)。证书充分性不需rank/Slater；该文顶点与必要乘子证明另需可行、rank(A)=n及目标下有界。退化点只要求存在一个最优working set，不要求每个活动子集乘子皆非负。

5. 有限非负线性系统P-8b6157c4cbb03498：Ax=b,x≥0与Aᵀu≤0,bᵀu=ρ>0恰一可行；最小残差z=b−Ax*给bᵀz=||z||²，非零时ρz/||z||²是归一证书。nullspace参数化x=x̄−Kᵀy的显式双射另需rank(A)=m≤n；维数m≤n本身不足。近零残差及小分母不能作浮点可行性认证。

6. 锥/非线性/无限指标情形先定义对偶锥、拓扑和乘子空间；弱对偶、零间隙及对偶达到分别验收。

7. 锥LP来源P-08ccea72ef78920a用h_c(b)与h_c**(b)区分原/对偶值：有限值下零gap检h_c在b的下半连续，dual达到另需∂h_c(b)非空。保拓扑、连续A及正对偶锥；PDF的bar h_c proper不能误抄成h_c proper。Gale全无限指标的liminf尾与截距inf不能由有限采样认证。

8. DC-composite来源P-18e365da3a6d391f的sup_λ inf_μ必须保留一个共同λ。其Theorem3.6印刷∀μ∃λ不足以推出∃λ∀μ，缺量词交换证明时不调用等价。

9. 输出原与对偶可行候选及gap；浮点残差不冒充精确不可行证书。

10. 非闭锥来源P-5ae169908ed70fca限Hilbert与bounded linear A、bounded closed convex generator K∋0：cl A(coneK)的零support方向测试只证明∀ε>0近似可解；exact A(coneK)另需一个共同C使∀y ⟨b,y⟩≤Cσ_K(A*y)。ε>0 dual唯一不推ε=0 dual达到；general support-face inclusion需同时回验残差/normal，H/E恢复原非凸锥是额外条件。该文Remark1.7 upper semicont错误及Prop4.1平方根缩放错误按卡局部修复后才可调用，均尚待两轮审核。

## 不可省略的限制

有限维Farkas不证明一般锥/非凸强对偶。DC-composite原文存在已定位的量词缺口：逐μ选λ不是共同乘子，须修正或另证minimax/鞍点条件。 该有限LP扰动证明的共同小ε与固定working-set子序列已做局部补足，原文Step1计数须加终端检查；不得将rank假设移植为Farkas矩阵限制或无限维强对偶。

本项目前仅提供操作方法与待核验来源；高级定理及常数须从真实阅读卡回查原文。
