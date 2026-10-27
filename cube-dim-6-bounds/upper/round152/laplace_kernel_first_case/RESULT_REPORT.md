# C152：真实group-reduced单项的唯一Arb证书已执行

Terra在其独立brascamp108 clone中完成唯一scalar case，状态
PASS_ONE_REAL_M2_CLOSED_LAPLACE_CERTIFICATE。h=1/128、lambda30，
Arb512，结果保存完整向外有理端点。原PROOF、single_case及
HANDOFF源均未改；本报告仅追加实际结果。

## 准确对象与尺度

147已证的真实cone0/body pair(16,17)条件联合律为
A=E5/60+E6/30、B=E5/60，密度1800exp(-30(A+B))，0<B<A。
其完整密度质量1，输出并不独立。对称性给
E_C[k_h(A)k_h(B)]=900 I30^2。

乘130因子与alpha0后，此次认证的量是

 (46080/pi^3)*(1/171992678400)*900 I30^2.

它严格指group-reduced21-cone正表达式内cone0的一个pair项，
不是原始固定pair在全R6的事件积分。群作用会置换pair标签，
约化因子不能据此单项解释成固定原始pair的全空间事件因子。
SCOPE_ADDENDUM_v2的hash已被frozen_sources绑定。没有额外乘496，
也没有自动扩展为其他pattern或整个beta2的上端。

## 保存的向外有理证书

certificate.json中的lower_fraction/upper_fraction是所有决定
使用的准确有理数；以下仅摘其ball_display_only：

| 量 | 展示中心值 |
|---|---:|
| p_h | 0.09947644966022578595859125982094048… |
| I30 | 0.002487396649868377484054276845476066… |
| 900 I30^2 | 0.005568427884398784920605257996531015… |
| 46080/pi^3 | 1486.15070668183246161816818788901477… |
| group-reduced cone0单pair贡献 | 4.811555417643914944150628519645254…e-11 |

部分贡献准确下/上端共用分母
230344386280611654799899571593522271174128492457058131614758685317726084130380178960334883128519211886075950624436269830308263222041308144077235977260296389382897664，
分子分别为
11083147797323396823530816906570591659081948125962797383284408446680369036096000190547513569250160535195319558784326086986097075928945512632654685424260897
和
11083147797323396823530816906570591659081948125962797383284408446680369036096000190547513569250160535195319558784326086986097075928945512632654685424260963。

唯一producer通过严格正性、I_upper<p_lower/30、E_upper<p_lower^2
及固定10^-80区间宽验收。本文只读保存值及来源/收据链，没有
追加数值求值或运行完整checker；后续独立saved审计由主安排。

## 精确门控、调用与实际收据

main release与独审final gate同时为true，绑定single_case SHA
6a2f8da6bef12d575fd2b7ba99aa1e00bfa7d2e02c77a88732f77575c083483c。
独审final_manifest SHA为
8739dc94a16be8cd74bfaee92469fcff13b83ae217c0d87b582c9d07681ea86e。
实际frozen_sources SHA为
6a56e84ad57c5bc99a3d4ce9fb95e47e1bd16370fb94ce6e513cd056f44506b6，
release SHA为
793a9c421dab53234e85ed22c4fe8c506e9fef8d0316e696e75e5865f9ea119d。

实际各1次sqrt(h)、erf、erfc、exp、pi常数及closed transform。
逐点k、numerical quadrature、density generation、旧343456
functional、MC/native均为0；lambda0/2分支未另求值。

operator_terminal_receipt记录PID83622、exit0、wall
0.11226091600838117秒，user CPU0.069797秒、system CPU0.023004秒，
stderr为空、无自动重试。producer内部wall0.05151720799040049秒、
CPU0.047554秒，计时边界不同。certificate SHA为
ad1676adffc4c1245d559c4f3cc3e55fd98470ff4edfe223cc499a36d0edfc82。

is_full_beta2_upper、new_M_upper、global_g0_proved及new_C6_bound
全为false。首例验证解析核与真实尺度接口；完整signed二维
积分、所有9961模式正权与公共pi的全域误差仍是后续上界工作，
不能由这个正部分项声称完整beta2上包或新C6界。
