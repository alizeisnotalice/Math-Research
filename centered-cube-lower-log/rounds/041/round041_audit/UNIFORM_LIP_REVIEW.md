# 连续uniform-Lipschitz factor6独审

结论：`../round041_mcd/uniform_PROOF.md`的自含证明通过。该改进确实使用坐标Lipschitz常数和连续uniform隐藏law，并非把旧oscillation-only McDiarmid的2直接换成6。可从已验收的source surrogate及逐尺度均值/variance数据做独立离线分支，不需新source、replay或MC；不得改写原factor2冻结文件。

先实际运行自有 `uniform_lip_checks.py/json`：穷举81个四段、每段宽1/4、斜率−1/0/1的1-Lip函数，405个Fraction stoploss检查、567个分片解析指数积分数值检查通过；14个非零affine情形达到uniform MGF等号。它们独立于生产方72函数测试，未使用tail payload。

## 1. 一维sharp MGF

连续1-Lip g的range中，任意a<b间可取最后离开a、首次到达b的路径；中间取值在(a,b)，长度至少b−a。故中间值集合的Lebesgue质量至少b−a。quantile端点的CDF关系再给q(v)−q(u)≤v−u，包括原子/平坦段，因此非降重排仍1-Lip。

中心化quantile a(u)与b(u)=u−1/2满足a−b非增、积分0，故prefix积分a≥b、suffix积分a≤b。用非降函数的stoploss等于suffix积分上确界，得到中心g(U)被U−1/2作convex-order支配。正铰链组合及紧区间凸函数逼近的方向正确，没有假定原g单调。

由此对任意实t，MGF≤sinh(t/2)/(t/2)。级数比较 `(2k+1)!≥6^k k!`的相邻比是(2k+3)/3≥1，给≤exp(t²/24)。affine的variance1/12说明普适二次系数1/24不能缩小。

## 2. Doob及source surrogate

给定prefix后，独立未来坐标的积分保留当前坐标的Lip；当前坐标条件law仍连续uniform。逐步使用sharp一维MGF，得到exp(t²Σ(Lip_i h_i)²/24)，Markov优化给exp(−6gap²/S)。均值、strict阈值、S0/空项及表间union与旧证明一致。

当前source中v_i,s=Lip(z_i,s)h，真实均值下包l_s和q_i=max(v/l)≤q0给min z_s≥(1−q0)l_s。因此 `Lip(log z_s)h≤q_i/(1−q0)`。log正和的坐标Lip不超过各log核最大Lip。平方求和与旧quarter r_i=q_i²/4的log product比较，确由相同C0以及固定C=4.05覆盖。

关键是**相同确定surrogate**仍有效；旧实际Σosc_i²不一定控制Lipschitz cost，不能拿更小的osc值直接替代。新variance/12对应的r/3也不能不付系数就替代这里旧quarter r。

## 3. 离线扩展的范围

固定产品使用R1/table0，joint/prior先同L组合均值，再用R表union；λ≤均值上包时返回1，不对负gap平方。要保留source contribution cap仍须与原Markov/旧upper取min，之后才作共同尺度union。factor6给同target的合法pointwise upper，按 `METHOD_UNIFORM_CI.md` 不因事后提出多收方法数费用；原1008预算仍可保守保留。

producer旧MCD小例是两点隐藏law，不能验证factor6；本证明的连续uniform前提与当前target的bin内隐藏law匹配，counts提议不改变它。对step函数或uniform域外的端点law则结论失败，反例口径正确。

完整mean/BE/r0及机器FP误差前提未消失；尤其正upper指数下溢成0不能当成认证上包。当前独审通过的是解析数学接口，不能仅由数值bound-chain断言升级为有限概率或Cn证书。

后续实际执行：root已用独立 `uniform_upper.py` 分支计算factor6 FINAL，exp exponent向上floor到−700；原factor2文件及receipt未改。本审最终1008 CI的Decimal复核见 `final_ci_checks.json`，最大端点差约2.16×10^-12，84个同阈值比较仍无正分离。floor只支付局部下溢方向，不代替完整FP认证。
