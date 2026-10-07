# 第50轮：深细化预算与真实逆向标签反例

日期：2026-10-07。目标：核查保留原粗尺度的比较费用能否支付真实外侧交叉项。

**证明状态：一般维数无关 weak-(1,1) 上界与发散下界均未证明。** 本轮严格结论为深细化平均费用不超过 9X/8、激活首矩平均不超过 5X/2；另构造一维光滑 L¹ 输入，证明任意小剪切位移仍可产生任意大的逆向标签跳变。后者排除标签稳定性路线，不构成发散下界。

## 0. 假设与记号

沿用第49轮，P 为任意概率测度，α>0，n 为正整数，有限深度 D。
V_j=ω_n R_j^n=8/(α2^j)，u_j=1_{B(0,R_j)}/V_j，g_j=u_j*P。
J 是 g_j>α/8 的首次指标，K_β 是 g_j>β 的首次指标。固定
E={α<MP≤2α, J≥2, K_{α/2}≤D}，X=α|E|。
β 在 I=[α/4,α/2] 均匀平均，E_j=E∩{K_β=j}，h_j=u_j*1_{E_j}。
这些 E_j 两两不交并分割 E，且 β<g_j≤2β。
记 B_β=Σ_{j<k}∫_{E_k}g_k h_j；已有平均 B≤3X/4。
O_β 是 R_cross=Σ_{j<k}∫h_j h_k dP 中 |x−z|≥R_j 的部分。
F_l={K_β=l} 在全空间定义，q=x+z−y。
C_G 是 O 中 j≤K_β(q)=l<k 且 |x−y|<R_l 的部分。
D_l=Σ_{j≤l}2^{j−l}1_{E_j}，并定义

C_D=Σ_{l<k}∫P(dy)∫_{E_k}dz∫_{F_l}dq u_l(q−z)u_k(z−y)D_l(q−z+y)，

T=Σ_{l<k}∫(u_l*D_l)h_k dP。

所有尺度指标至少为2，实际有限计算截断于D；无限和只作为非负Tonelli推广。以下给出完整论证，保留英文证明以便逐式审核。

## 1. C_D is a single selected part of the original source cross energy

Fix beta and original finite-depth E_j as in049. Put R_cross=sum_{j<k} integral P(dy)h_j(y)h_k(y), D_l=sum_{j<=l}2^{j-l}1_{E_j}, and F_l={K_beta=l} globally. Denote C_G the compatible OUTER part defined in049, and C_D its displayed formula after dropping the old exterior indicator.

Expand D_l in C_D and substitute x=q-z+y. On |q-z|<R_l, the factor 2^{j-l}u_l(q-z) equals u_j(x-y) times 1_{|x-y|<R_l}. The remaining original fine factor is u_k(z-y). Thus C_D is EXACTLY R_cross restricted to

K_beta(x+z-y)=l, j<=l<k, and |x-y|<R_l.

Each original (j,k,x,z,y) has at most one such l, because the F_l are disjoint. Therefore there is no new layer multiplicity in C_D, and C_D<=R_cross.

Its outer part |x-z|>=R_j is exactly C_G. The remaining inner part is bounded by B_beta. Indeed

H_jk(x,z)=P(B(x,R_j) intersect B(z,R_k))/(V_jV_k)
 <=g_k(z)/V_j,

and on |x-z|<R_j the right side is u_j(x-z)g_k(z). Dropping the extra compatibility restrictions and integrating x in E_j,z in E_k yields B_beta. Consequently

0<=C_D-C_G<=B_beta,
C_G<=C_D<=C_G+B_beta,
C_D<=R_cross<=B_beta+O_beta.

All boundaries are handled by the disjoint partition |x-z|<R_j versus >=R_j. The kernels use open balls throughout; no atomic P boundary approximation is needed. Equalities in observer-distance spheres also have zero dx dz measure, but this fact is not required for the nonnegative partition.

A uniform C_D budget is therefore equivalent, up to the already paid B_beta, to a uniform budget for this selected actual outer part C_G. It is not a genuinely stronger independent budget. It still leaves incompatible actual outer edges.

## 2. A new general paid part of the relaxed T: all sufficiently deep refinements

Dropping F_l gives

T=sum_{l<k} integral P(dy)(u_l*D_l)(y)h_k(y)
 =sum_{j<=l<k}2^{j-l} integral_{E_j}dx integral_{E_k}dz H_lk(x,z).

The term l=j is exactly R_cross. Other l terms are nonnegative, hence T>=R_cross. Deleting F_l loses uniqueness of l and creates repeated charging of the SAME original source edge.

There is nevertheless an absolute budget for the entire portion with l>=j+n. The sharper first-exit form is

T_deep,beta <=3beta|E|, hence average T_deep<=9X/8.

Proof for arbitrary P: on E_k the actual first exit obeys g_k(z)<=2beta (g_1<=alpha/4<=beta, so every exit has k>=2, and g_k<=2g_{k-1}<=2beta). Hence H_lk<=2beta/V_l. Put t=2^{-1/n}; because k>l, R_k<=tR_l, so H_lk vanishes unless |x-z|<(1+t)R_l. Since 2^{j-l}=V_l/V_j, the weighted kernel is at most 2beta/V_j times this indicator. For fixed l the fine masks E_k,k>l, are disjoint, so

sum_{k>l} integral_{E_k} 1_{|x-z|<(1+t)R_l} dz <=(1+t)^n V_l.

For fixed j integrate x in E_j and sum l>=j+n:

cost <=2beta |E_j| (1+t)^n sum_{l>=j+n}(V_l/V_j)
 =4beta |E_j| ((1+2^{-1/n})/2)^n <=3beta |E_j|.

The last inequality is Jensen for the concave function u^(1/n) at u=1 and1/2: (1+2^{-1/n})/2 <=(3/4)^(1/n). Summing j proves the claim, and uniform beta averaging uses average beta=3alpha/8. Finite truncation only decreases the geometric sum. Infinite depth follows by Tonelli, as an extended nonnegative inequality; whenever X is finite it is an actual finite budget. The only input is actual first-exit overshoot and original disjoint E_k.

For the remaining at most n refinements j<=l<j+n, kernel nesting gives

2^{j-l}(u_l*1_{E_j})(y)<=h_j(y).

Thus, at finite depth and then by Tonelli,

R_cross<=T<=n R_cross+3beta|E|,
average R_cross<=average T<=n average R_cross+9X/8.

All deep refinements are INNER original edges: R_l<=R_j/2 and R_k<=tR_l imply |x-z|<(1+t)R_j/2<R_j. For n=1 every extra l>j is deep. Thus in one dimension 0<=average(T-R_cross)<=9X/8; additionally the OUTER part of T is exactly O_beta (only l=j survives). Since the inner part of R_cross is at most B_beta, in n=1 one also has

average O<=average T<=average O+15X/8.

In general the n nearby refinements remain unpaid; this does not give a dimension-free bound for T or the main theorem. A dimension-free T budget would itself imply a dimension-free Q budget, since Q=d+2R_cross and d<=X. The deep budget removes a genuine repeated inside charge, but it does not pay any new actual outer edge.

### 2.1 Sharper source first-moment budget

The049 bound sum_l integral P(dy)(u_l*D_l)(y)<=4X improves using the actual first exit at the ORIGINAL x. For x in E_j, g_j(x)<=2beta, and ball nesting gives g_l(x)<=2^{l-j}g_j(x) for l>=j. Also the original height cap gives g_l(x)<=2alpha. Consequently

sum_{r>=0}2^{-r}g_{j+r}(x)
 <=2beta+2beta+sum_{r>=2}2alpha 2^{-r}
 =4beta+alpha.

The first two terms use g_j<=2beta and (1/2)g_{j+1}<=g_j. Integrating x over the disjoint E_j and using Tonelli yields

sum_l integral P(dy)(u_l*D_l)(y)
 <=(4beta+alpha)|E|<=3X.

Uniform beta averaging sharpens this to5X/2, because E is fixed and average beta=3alpha/8. Finite truncation only reduces the sum; infinite-depth statements follow as nonnegative Tonelli identities whenever X is finite, or extended inequalities otherwise. This first moment still does not pay the correlated future-occupation multiplier in T.

## 3. Actual one-dimensional shear need not preserve the label

The same-label phenomenon in049 finite tests is not geometrically necessary even in n=1. The following actual original-band family has arbitrarily large backward jumps while its displacement tends to zero.

For an integer j>=5 put s=2^{5-j}, k=j+1, D=j+2, alpha=1, and

P_s=.19 delta_{(-.5+.03s)}
    +.04s delta_{(-.015s)}
    +.06s delta_{(.105s)}
    +.19 delta_{(.5+.02s)}
    +(.62-.1s)delta_10.

It is a probability measure. Choose x=0,z=.13s,y=.105s; then q=x+z-y=.025s. Here R_j=.125s, R_{j+1}=.0625s, R_{j+2}=.03125s, while R_3=.5 and R_2=1. Take any common beta in H=[3/8,379/1000]=[.375,.379].

### 3.1 All original clocks and the actual height band

The near mass is .38+.1s<=.48. Its entire support is inside the R_1 and R_2 balls of x,z,q, while the far atom is outside. Therefore g_1<=.12<1/8 and g_2=(.38+.1s)/2>=.19>1/8. Thus J(x)=J(z)=2 and K_beta is not 1 or2.

The R_3 ball of x contains the left .19 atom and both fine atoms, but excludes the right .19 atom. The R_3 ball of z contains the right atom and both fine atoms, but excludes the left. Thus g_3(x)=g_3(z)=.19+.1s<=.29<beta.

All radii from R_4 through the fine scales exclude the coarse .19 atoms. For x, both fine atoms persist through R_j, giving g_j(x)=.4 and g_{j-1}(x)=.2; earlier fine densities are smaller. The next radius R_{j+1} excludes the .06s atom but contains the .04s atom, giving .32; R_{j+2} gives .64. Hence K_beta(x)=j and K_{1/2}(x)=j+2.

For z, both fine atoms persist through R_{j-1}, but R_j excludes the .04s atom and retains the .06s atom, giving g_j(z)=.24. Then g_{j+1}(z)=.48 and g_{j+2}(z)=.96. Hence K_beta(z)=j+1 and K_{1/2}(z)=j+2.

The maximal functions can be evaluated at the finitely many atom-distance endpoints. At x, the central .04s atom gives 4/3; the two fine atoms together give 10/21. The subsequent coarse-atom cumulative ratios are at most .29/.94 and .48/1, while the far threshold is at most 1/20. Thus MP(x)=4/3. At z, the .06s atom gives 6/5; the two fine atoms give 10/29; each later cumulative ratio is below 1/2. Thus MP(z)=6/5. Therefore x,z are strictly within the ORIGINAL band and eligible at depth D.

The exterior and capture conditions are also strict:

|x-z|=.13s>R_j=.125s,
|x-y|=.105s<R_j,
|z-y|=.025s<R_{j+1}.

### 3.2 The backward label

At q, the R_3 ball contains BOTH .19 atoms. Each lies at distance .5-.005s from q. It also contains both fine atoms, so g_3(q)=.38+.1s>.38>beta, while g_2(q)<=.24<beta. Thus

K_beta(q)=3,
K_beta(x)=j,
q-x=.025s ->0,
j-K_beta(q)=j-3 ->infinity.

Indeed MP(q)=5/8: the two fine atoms give this value; the coarse cumulative ratio is at most .48/.99<1/2, and individual/far ratios are smaller. Thus this particular shear target lies below the original height band. Its backward label is still a genuine GLOBAL first-exit label. It does not satisfy j<=l, so it belongs to the remaining incompatible part rather than C_G or C_D.

### 3.3 Positive-area and actual smooth L1 persistence

Take eps=s/10000 and observer boxes |x|<eps, |z-.13s|<eps. With the same source y, q differs from .025s by at most .0002s. The limiting R_3 capture margin at q is .005s, so both coarse atoms remain strictly captured; all canonical memberships used above remain unchanged. Their other smallest relevant margins are .005s for exterior separation and .00625s for the fine K_{1/2} capture. Consequently all exact canonical clock values remain the displayed ones throughout the boxes, for every beta in H.

The central-atom ratio at x is between .04/(2*.0151) and .04/(2*.0149); the fine-atom ratio at z is between .06/(2*.0251) and .06/(2*.0249). All other cumulative ratios remain below1/2. In particular throughout both observer boxes

11/10<MP_s<3/2,

and the distance to every original atom is at least d*=.0149s. Thus these are actual original-band observers, with the original canonical J/K eligibility.

For genuine smooth L1 inputs, convolve P_s with ANY nonnegative smooth probability kernel supported in (-h,h), h=s/100000. For x,z in these boxes and y anywhere in the convolved .06s source blob, q differs from .025s by at most .00021s. All canonical captures and exclusions remain strict after the additional blob spread h. Therefore the canonical clock values at x,z,q are EXACTLY unchanged; in particular J(x)=J(z)=2, K_beta(x)=j, K_beta(z)=j+1, K_{1/2}(x)=K_{1/2}(z)=j+2, while K_beta(q)=3.

Here is an elementary uniform height comparison. For any observer a at distance at least d* from every original atom,

[d*/(d*+h)] MP_s(a) <= MP_{s,h}(a) <= [d*/(d*-h)] MP_s(a).

For the lower bound, enlarge any original atom-capturing ball radius r>=d* by h, which captures all corresponding blobs; r/(r+h)>=d*/(d*+h). For the upper bound P_{s,h}(B(a,r))<=P_s(B(a,r+h)); a positive numerator forces r>=d*-h, so (r+h)/r<=d*/(d*-h). Take suprema, with arbitrarily enlarged radii if required by open-ball endpoints. Thus throughout the observer boxes,

(11/10)(1490/1491)<MP_{s,h}<(3/2)(1490/1489),

which lies strictly in(1,2). The actual smoothed height band and original J/K are therefore reconstructed, rather than frozen. The fine source blob has mass .06s, and all its points are captured by both original balls and have backward-label shear targets.

Both atomic and smoothed versions have a concrete averaged uncovered outer cost. The x,z box area is4eps^2, the source mass is .06s, the kernel product is 1/(V_jV_{j+1})=32/s^2, and the threshold weight is |H|/|I|=2/125. Hence

average O_uncovered >= (.06s)(4s^2/10^8)(32/s^2)(2/125)
 = (12/9765625000)s >0.

Uncovered means q has l=3<j and hence fails the global compatibility condition. This is a genuine positive source/observer/threshold-measure failure in smooth L1 inputs, but the displayed cost decreases with s and is not a divergent global-budget example.

## 4. Exact remaining scope

New paid general part: all relaxed refinements l>=j+n cost at most3beta|E|, hence at most9X/8 after uniform threshold averaging. In n=1 all extra relaxed refinement cost is paid relative to R_cross. C_D itself uses a unique global label and is simply a compatible subpart of the original cross energy, differing from actual C_G by at most B_beta.

The actual 1D family disproves universal same-label or even uniformly bounded backward-jump rules; arbitrarily small shear displacement can cross arbitrarily many genuine global first-exit levels. It does not produce a divergent weak-type lower bound or prove C_D/T divergence.

Still open here: a dimension-free bound for the n nearby relaxed refinements, a dimension-free T bound, all actual compatible outer costs, and the remaining incompatible outer costs.

## 5. 本轮数值核验与复现

运行 `python3 rounds/050/verify.py`（在项目目录内；Python 标准库），完整可再生结果写临时目录，控制台输出最小核验收据。`verification.json` 保存本次实际执行的摘要及代码 SHA256。

- 一维有限原子输入共12例，分3批4+3+5。使用 Fraction 对空间多边形面积及共同β的三方首次记录交集精确积分。核验 C_G≤C_D≤C_G+B、C_D≤R≤T、T−R≤9X/8、T_outer=原O、平均激活首矩≤5X/2。R和B另与第47轮独立卷积实现比对；11例C_G与第49轮收据比对。
- 逆向标签族核验 j=5,6,8,12,20,32,48,80,128，三个批次。使用有理数检查盒顶点、β端点与平滑位移端点，共288次标签检查。端点检查只核验明确公式；任意j、盒内部和任意平滑核由第3节的严格余量证明覆盖。
- 27个有理首矩标量检查；n=1,2,3,8,32,64,128,1024,4096 的深细化标量因子用50/80位Decimal交叉检查。后者不是区间证书，也不是高维体积核验；全维界来自凹性证明。

最大 C_D/X=512639/16777216，最大 T/X=297151/8388608，均出现在 four_atom_perturbation_s471016。α=1 的链 L=2,4,8,16 中 T/R 约从2.07577到2.14857，未见持续增长。有限数值不能推出一致上界或渐近极限。

## 6. 下一步及归档范围

应直接控制保留来源权重的近层相互作用，以及逆向或无标签的真实外侧余项；不能假设剪切保持标签，也不能仅凭首矩预算控制与未来占用的相关乘积。本轮深细化已付项全部属于原粗尺度内部，尚未解决新增真实外侧费用。

仅归档本报告、两个最终探针、复现入口与最小核验摘要（5个新文件），并更新总览。依赖复用先前轮次的最终代码。草稿、全量再生JSON、ZIP及提交规则附件不纳入本轮；未改动其他项目。
