# 11O：三维积分的高阶全域有限批（新表示，不扩旧自然box网格）

2026-10-25。11H初始与唯一修订均因signed-cone自然box取消膨胀未闭合，停止旧表达。新机制复用11J已本地成功的准确有理point Hermite DD+Chebyshev/Bernstein连续插值，积分验收改为凸弦与准确多项式mean，而非节点或max粗和。

原bridge独立clone只upper/continuous_bridge/cb11_integral_cheb_*；旧11H先锁manifest，不覆盖。固定τ1/3,Hε,ε=9/5000；128×128方格三角半区8256patch，每维degree32、33roots、rho5/2。每patch半宽1/256，s=h(rho−rho^-1)/2，δ=ε−s/3>0，全纯M约2373.63073，二元插值误差Epoly<7.707e−8（用Arb真实包围）。root候选目标更强g_low=131/25000（11J作者证书支持，033独审进行中）及I_target=29/50000=.00058，条件成本约7.94534839<8，仍需全部认证才改U。

主整face短time预算11M：n3Etime=2^-4/3*(3/4)*[(256/3)/pi³]^(1/3)*[(16/3)/pi]^(2/3)*ε^(4/3)约.00012996912955。故I_H允许上限约.0005638179328648。必须准确计算I≤(I_H^(2/3)+Etime)^(3/2)，不直接用浮点阈值。

每patch取P为真Cheb interpolant，误差≤Epoly；从区间Bernstein取得准确L=max(0,min_coeff_lower−Epoly)、U=max_coeff_upper+Epoly，以及真H的平均M_H≤mean(P)_upper+Epoly。mean(P)按μ_0=1、μ_even=-1/(k²−1)、μ_odd=0准确对Cheb系数积分；可截M_H≤U。q=3/2凸性给整个patch平均H^q≤L^q+(U^q−L^q)(M_H−L)/(U−L)。L=0/U=L单独安全处理；全部值和面积准确包围。此为联合连续积分，不是节点计数。

triangle镜像积分须权重正确：j<i方块权2，j=i整diagonal square权1，乘每块面积1/128²；总权面积严格1。不要像minimum覆盖那样将diagonal加倍。全部真实方向奇线由Hε平滑而在多项式内覆盖，G−Hε由整faceLq预算覆盖，无删帽/漏线。

先4代表patch(127,57),(127,127),(127,0),(0,0)端到端；节点优先384bit和2^-128准确近似，加Cauchy节点位移误差。实际记录时间/point半径/Bernstein转换/mean/chord。通过完整接口后执行固定8256patch，最多12本机workers（按当时CPU内存可减少，不是总计算硬限），明确每patch有限输入，无额外分支/τ搜索。若误差失败可仅同节点精度一次提升512bit，不扩网格。总约8990784nodes，缓存复用4patch，不提交全节点缓存；输出必要patch证书与程序，计算资源按实际负载。全批未过即准确报告I上界/未闭合项，不机械加密。

I作者自检不替代独审。032表示、033minimum、11M整face分析均作为明确依赖，后者独审未完成，不宣称新U3。MATLAB由主串行，主承担最终全L1/allr和精确成本/归档。禁止自行扩招/提交推送/其他目录或总账写入。
