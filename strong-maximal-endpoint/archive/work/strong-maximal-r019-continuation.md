# 第19轮接续：窄平台实际下界（外部研究笔记，待正式报告/最终审计）
2026-11-06。用户未叫停，goal active。r017/r018已完成严格推导并分别提交；r018push已确认。

## 当前权威成果
专用clone /Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal，project strong-maximal-endpoint/，当前main bf02e0d。
r017 03731aa已推送：真实固定斜率ramp的比例偏移改善g logn，g=(mu−1)/8。
r018 bf02e0da715a9344a8545a56b59986b8dc57c908已提交并推送，远端包含已确认。upper_audit审查通过。
R_ramp(mu)~sqrt(2/pi)e^b n^-1/2 exp[nmu−dn+a z_n²]，a=(3−mu)/6,b=(5mu−6)/6。
qsoft(mu)=nmu−dn+a z_n²+b+o1；I_n~sqrt(pi n/2)。
一般C下界目前repo登记nmu−dn−(mu+3)L/24+oL；上界仍nmu−dn+O(n^3/8 L^1/8)。
FGL原始已核验；r018自包含FourierCDF绝对O1/n证明支持家族乘法渐近，不是一般C等价。

## r019关键新构造
固定斜率1/n，轮廓H_w(u)=0(u<0),1(0≤u≤w),exp((u−w)/n)(u>w)。
对输入h(t)=H_w(t−q)，真实J=e^-q K_n(w)，K=1−e^-w+e^-w I_n，完全独立c。
取w_n=.5 logn+loglogn，则K=1+O1/logn，logK=o1，省去纯ramp分母I~sqrt n！
接收均值B_{c,w}(q)=A_c(q+w)+F_{n,c}(q+w)−F_{n,c}(q)。
其精确q导数B'=-A_c(q+w)/n−f_{n,c}(q)<0，故有唯一root qpl(c,w)。
A为r017纯ramp均值，qsoft是A=1根。

## 直接证明路线（已与multilayer_theory独立核验，待写正式报告）
均匀c=mu+tρ、有界t，r017已有qsoft附近f~L/n、A~1，λ~ρ/sigma。
对q=qsoft+O1，差D=A(q)−B(q)精确为
 (e^{w/n}−1) A(q+w) + E[(e^{(S−q)/n}−1)1_{q≤S<q+w}]。
第一项~w/n；带内项≤O(w/n)P带≤O(w/n)*O(w f)=O(L³/n²)=o(L/n)。
A(q+w)=1+O(L²/n)=1+o1，密度在O(L)带的相对变化exp(Oλ L)=1+o1，可由r016倾斜全局导数/局部Fourier界补齐。
故D~w/n，在qsoft处B−1=−w/n(1+o1)。
在qsoft±常数区间B'=−(L/n)(1+o1)；先比较qsoft−(.5±epsilon)的B符号即可把root夹住，无需假设root近。
得到qpl=qsoft−w/L+o1=qsoft−.5+o1，均匀有界t。
因此fixedmu平台score=exp(qpl)/K~exp(qsoft−.5)，相对ramp评分增益~sqrt(pi n/(2e))。
实际固定mu平台乘法式：exp[b−.5] exp[nmu−dn+a z_n²]。
结合r017 c偏移gL且平台根损失均为O1，一般实际下界应升级为
 logC≥nmu−dn+(a+g)L+oL =nmu−dn+(9−mu)L/24+oL，系数为正！
这仍是下界，不是C/Q最优等式。strict q<qpl、有界截断、固定n物理深层极限依旧需要逐项写清。

## 已实际运行本机计算（勿无理由重复）
/work/sm_plateau_probe.m、plateau_probe.json：MATLAB R2026a固定mu，w=.5L+logL，预先生成原始与正倾斜binomial众数±12sd权重，40步root二分。至1e7共约40秒。不是定向证书。
输出plateau score−qsoft：n256 −1.47227115；1024 −1.41802299；4096 −1.34779934；1e4 −1.3020246；1e5 −1.1949771；1e6 −1.10814924；1e7 −1.03908076，理论−.5。
root−qsoft分别−1.2711224、−1.25300181、−1.20797239、−1.17475519、−1.0917253、−1.0213352、−.964203678，慢收敛需保留所有数据。
初始lo=qsoft−10,hi=qsoft+1；root residual保存在JSON。
下一步将最终脚本适配repo输出、Code Analyzer、Wolfram60/100位独立核验n256/1024原始binomial公式，再写r019+独立审查、逐轮commit/push。
理论外部笔记 /work/strong-maximal-r019-plateau-soft-lower.md （multilayer_theory已接受窄平台机制）。

## 上界线独立审计结果（仅启发式，严禁作为定理）
global_upper_adversary在/work/strong-maximal-global-adversary.md核验formal Bellman diffusion：局部控制修正D=mu(mu−1)/8，Duhamel/Brownian bridge给CDF quantile feedback gain≈g[z²+logn]+O1，fixedc仅g z²。
现有r013FGL精度远不足以证明logn尺度，需要累计o(logn)残差的离散Bellman屏障/adaptive local limit。
本机r011/r013 n256/1024实际反馈增益.29745/.36723，formal预测.22566/.30214，差约.07，不能用两点证明渐近。
后续一条路线应继续一般upper而非永远只做受限family。

## 更远探索（纯线索）
KKT完整最优profile由a(t)=e^t f_S(t)单调决定：h=0 ifa<ν，h=1 ifν<a<nν，h>1时Phi'(h)=a/ν。平台宽~logn，softshape非纯指数，可能进一步提高常数甚至不同尺度，需独立证明。
当前窄平台只优化固定斜率1/n的一类，不宣称全KKT最优。
