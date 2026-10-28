# R077 论文冻结前先验审读

日期：2026-10-28。责任：独立数学审读；仅写本目录。全文冻结审核尚未开始，本文不登记新论文为 FINAL ACCEPT。

## 已读材料与本次依据

已全文读 R073 `upper/round073/cantelli_frontier/PROOF.md`（作者 clone `cube-r040-cover-selector`，SHA256 `50d48b85f02a14f50614f0eea5cedbd7ccd68e9b2f7e1564fe8b41fcd2af9c4b`）、本 clone R073 `review.md`，以及冻结局部依赖：

* 相对 J / AC / baseline：`c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1`。
* J1：`67bdf01d9aa65c70bf058df08b527b904e954d19f3327ecbb43e5b5c316b5030`。
* J2：`4d2f4fc2d0f3353ba82e32487dba68804d3021d4c0d477137c11f6c1a004e2b8`。

另读 R068 全文审计记录，并直接核对当前作者 `paper/main.tex` 中原有限选择、壳积分、有限 a 公式、一般 L1/正半径传递及实际测度/能量引理（至该引理证明结束）。该主文件仍在修改，本次未把它当作冻结版本，也未完成其整篇审读。

使用数学论文写作技能的 `claims-and-proofs.md`、`computation-and-certificates.md`，并实际读 `/Users/zhengzhihao/.codex/skills/anti-defensive-writing/SKILL.md`。适用原则为直接陈述结果、保留有效性所需条件、删除无信息的重复免责声明。

## 先验结论

现有 R073 证明与上述局部依赖的参数域和概率律相容；本次逐段检查未识别新的实质缺口。它可整合为对所有整数 n、自然对数 L=log n>=8192 的原中心轴平行立方体算子定理：

    C_n <= F(A_n(t), B_n(t)-(291/1000)n)
        < n log n-n log log n-(47/250)n,
    t=n/[(231/256)log n].

独立冻结全文仍须验证作者写出的数学内容、依赖、比较及范围；旧审计状态不替代这一工作。

## 必须在转写中保持的连接

1. **参数域。** qmax=11569/12800 与 1-sigma=231/256 分别控制源带与标量价格；二者不相等。新 h=1/(768L) 的门槛由新有限 guards 保证。R070 的局部有限 h_c 区间不是本 tuple 的全链门槛依据。active source 才使用 R>=1；inactive source 的列为零。

2. **真实相对域与 baseline。** R=sq-h、a0=(1-epsilon)q、g=R-a0，g-R/3=epsilon q-2h/3。完整匹配区间为 [(r-h)/s,(r+h)/s]，baseline 在右端 q_b=(r+h)/s，满足 R=r。所需条件是 g/R<=1/3，不能换回 g<=1/3。

3. **三种 profile。** J 的精确闭式保留中央删除 -m；J1、J2 的中央 phi=0，闭式无 -m。J2 是单坐标积分中 phi²/ell 的 profile。prefactor R/rho 必须随 q 微分；kink 用闭区间绝对连续性连接 a.e. 导数。实际强导数为 3/rho、1/(2rho)、1/(10rho)，对应总桥 6Lh、xh、xh/5。

4. **target law。** 固定完整 signed target face 上有 n-1 个独立 Uniform[-r,r] 自由坐标及一个确定 pin。Fubini 仅在 baseline 以 target capture fraction 取消一次 inverse ell。它给 target profile 均值等于 fullsource phi、phi² 均值。

5. **二阶均值下界。** J2 的正下界使用 a-5/6=h/(3R)-epsilon q/(2R)∈[0,h/(3R)] 与 R/q=s-h/q。相对上界 g/R<=1/3 单独不能推出正下界。n-1 个自由坐标与确定 target pin 分别记账；pin 的 1/(250rho_b) 半值只属于真正的 target pin。

6. **六事件。** baseline log、两个端点 count、target D、J1、双侧 J2 共六项，每项失败概率 <=epsA/8。J2 两个尾共用一次 Chebyshev 事件。事件交集用并集界；无需六事件相互独立。target cap 的 (3/4)^(n-1) 是每张 positive pinned face 上界；negative face 已在 D。

7. **全匹配计数。** ellpin'<0，Nsum 随 q 增大，故 Nsum/(2n ellpin) 增大。右端 mean<=1、左端短缺<=2h+1/n。逐坐标 J 的误差与 n-fold product 的 log 桥是两项估计；6Lh=6/768 必须进入 P0。

8. **固定 source-face 条件方差。** 给定 target w、source q 及 captured signed source face j，自由 inverse 坐标独立，均值 J1_i/J_i、raw second moment J2_i/J_i，pin phi 确定。正确方差为 Σ_(i≠j)[J2_i/J_i-(J1_i/J_i)^2]。方差上包可加回省略的非负 J2_j；不能以 target pin 半值代替任意该坐标。

9. **Cantelli 的顺序。** 每 face 用均值上包 p+m_*x、方差上包 v_*x 分别应用 Cantelli，再混合概率。retention 的全尾性由 (a x-p)^2/x 导数 a²-p²/x²>0 给出。以上概率先在 fullsource inverse law 下完成，随后处理 soft 权与 D 删除。实际增量 mu_a 不需要该概率律。

10. **同一个核的完整 row/column。** arad=r/rhi 与 bsoft 均在 [0,1]。row 使用对数源宽 log((r+h)/(r-h))，column 使用普通 target 宽<=2h。column 保留所有 target sign，并用 T_i/(2ell_i(r))<=1 与 product ell_i(r)<=G。k=nh theta/rhi；rhi 不能从归一化里消失。

11. **两种尾。** 源 cap 删除按未加权环境质量支付，源带质量为 C_log*n，其中 C_log=log(qmax/qmin)。cap/k 保留 (1025/1024) C_log^+ rhi*768L，指数率 cS=13627383/3790929920；同一全局 tau/k 保留 rhi*768L^(3/2)，指数率 cT=1447/327680。两率不同。有限 n 的 exp(eta*x/n)<=1025/1024 明确支付。S64 完整正前缀及相应导数不等式建立整个 L 尾；sqrt(8192)<91 只服务基点。

12. **同一有限历史。** a>=8n^4、source rho>4、实际 target item radius<3，保证 included source 真正先于同一 target。中心位移余量 zeta-(beta+2)/(a+1)>0 使 included pair 满足原资格测试。实际行由原 item weights 支付，每块 increment 只支付实际包括的部分。

13. **单次缺额。** D 内 lambda_a=c_a lambda；energy、degree、target mass 分别缩放 c_a²、c_a、c_a。A 与 P 是点集不交，因而 nu_a(A)+nu_a(P)<=Delta_a。能量恒等式将两项缺失质量合计一次；该条件不能只以 A×P 的形式支撑替代。

14. **原算子及完整根式。** 原 cover factor、捕获预算 tau、terminal self term a^-n 全部保留。固定 n 后只令有效数值不等式的 a→∞，无须 selector 收敛。一般 L1 采用正的实际见证半径、紧集有限子覆盖和内正则性。选择 tau=sqrt(A(B-gamma*n))，最终是完整 F(A,B-gamma*n)，含正交叉项。

15. **价格。** 对数价格是 -log(1-sigma)<257/2500，gamma-price-47/250=1/5000 支付 1/n、A/n 及根式交叉项。不得只以 gamma 或其展示小数宣称最终成本。

## 与旧结果的比较

R073 在 L>=8192 给系数 47/250；R068 的既有 profile 结果在 L>=2^30 给系数 23/100。前者提前起效，后者在其域内给更大线性扣除。若两条保留，论文应据此陈述门槛和系数的取舍。完整成本排序仍需针对各自 t 与 gamma 的 F，而非独立比较两个 gamma 或线性展示系数。

## 去防御性文字与冻结审核方式

论文可在定义处直接给出所用 probability law，在证明处写明先逐 face 估计再混合，并在实际测度桥处使用确定的 row 资格条件；这样能删除后文重复的“不是实际 mu 的随机假设”提醒。定理域、baseline、pin、strict level、全尺度和 L1 传递仍承担数学职责。

全文审核时将检查摘要、引言与主定理是否准确反映每个常数自己的域。纯语言修改过的旧冻结块将比较差异并审其数学范围；不会宣称它们仍逐字保留。待收到 source 路径及 SHA256 后，全文从头至尾读审，补查全部交叉引用、数学比较和必要附件可得性。

## 实际新运行

本次仅运行 `check_adaptation.py`：Python 本机 SymPy 1.14.0，PID 7280，12 项准确符号/有理恒等式全部通过。检查相对 guard、finite-h 结构、pin 负导数、三条强桥系数、两种准确指数率、有限 a 能量缩放余项、完整根式恒等式、Cantelli 单调导数及价格余量。

这些恒等式用于新稿转写适配；未重跑作者 53-gate 程序、R073 旧独立证书或旧局部核验。未做参数搜索、新研究、MATLAB 运行、编译、提交或推送。实际 source hash 与结果保存在 `adaptation_identity_receipt.json`。
