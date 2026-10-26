# E096：联合 middle 重排预算的独立数学审计

## 裁决

**PASS，限于作者准确限定的两个分别求解的抽象松弛预算。** 新的 mass/peak/volume
全时间重排不等式成立，comonotone 族在该单独信息类中达到上预算。六个保存的
抽象 profile 给出两个预算各自的严格下界；两个下界都超过 062 的固定充分余量，
故预算的 min 也超过。准确的条件性 S_plan 上端严格负。

这个结论没有给实际 Green 范数或模型成本的下界，不排除同时施加全部 tube
交叠约束后的更强优化，不排除实际 body32 模型，也不是新的 C_n 界。
g0=1/65000 仍未全局证明，far 仍未付款。最终证明与集成责任归主线程。

固定对象为原等权 body32、同一 semigroup 时间 t、Lebesgue 正面、共同半径、
q=6/5、rho=1/64、31 个 mean_B 在 (1/4,3/4] 的 middle tubes。继承上估计
安全包含空间端点 1/4，不产生漏测度或新的模型。没有将其他 frame 或不同
时间参数引入本轮。

## 封存与独立执行边界

作者目录为 `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round096/joint_middle_budget/`。
先登记再核验：

|证据|SHA-256|
|---|---|
|096 manifest|4381aff5a80fe53e3aa6c692e718fa288402e52604e4a7bcfe7c521aaec2dd68|
|096 results|e60253fa88be251156ca24c1bf58498585d17ebb60b1329a2b30a07dc4ebbbfb|
|096 theorem|3915219aa91a38a839af62f15290e26a4df7fe12fcc030b1b72df251723655d8|
|096 evaluate.py，仅静审|106738ef9c033c84ef0c19aacf4208854579a49825f36657fd1251255e0a3752|
|092 lower_results|fb342512ab2198b7a104ef78777d8e1856e15cca582d3e9b21d3b976ef2b4179|
|074 theorem|a2f7ee97758117244af0588db616721620c76749ebd02aacd856189b454e5717|
|062 results|ffb4c1ac8336fe662f6040dc2c6117eabf58f9d401b8631e8c1aacf548d241aa|

沿用本人独立 no-checkout clone `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit`，
有独立 `.git`/index；本次只写 `audit/round096/joint_middle_review/`。
097 封存产物保持不动。主 AGENTS 已在先前本会话独审时读取，本次再次将其 hash
纳入门禁。作者全部 manifest 文件、11 个冻结依赖及主 AGENTS 共 25 个源文件，
运行前后 hash 全部相同。没有执行或导入作者 evaluate/run_guard 或任何旧科学程序。

理论来源实际阅读范围为 096 theorem/report/evaluate、074 theorem、092 supplement
theorem、043/049/051 的模型与边际接口、059 区域总账，以及 062/092 的指定保存
标量。074 的 F/L/V 和 092 的统一下系数作为已有证明依赖使用；本轮审阅其适用
推导与准确数值传递，没有重新生成旧 F、L 值或重开其 86 格时间表。

## 全时间重排定理的独立证明检查

设 X 是体积至多 V 的非原子空间，可补零到体积 V。f_t>=0 联合可测，
f_t<=F(t)>0，且 integral_X f_t<=M(t)。定义

    theta(t)=min(V,M(t)/F(t)),
    h(s)=integral_(0,infinity) F(t)*1_{s<theta(t)} dt, 0<s<V,
    H(x)=integral_(0,infinity) f_t(x)dt.

对任意可测 A、a=|A|<=V，Tonelli 与峰值/质量约束给出

    integral_A H <= integral min(aF(t),M(t))dt
                 = integral_0^a h(s)ds.

等号中的 min(aF,M)=F*min(a,theta) 用到了 a<=V。h 随 s 单调不增。
对 lambda>=0，正部的集合上确界表示给出

    integral_X (H-lambda)_+
      =sup_A(integral_A H-lambda|A|)
      <=sup_(0<=a<=V) integral_0^a(h-lambda)
      =integral_0^V(h-lambda)_+.

再用 z^q=q(q-1) integral_0^infinity (z-lambda)_+lambda^(q-2)dlambda，
以 Tonelli 积分，得到 integral_X H^q<=integral_0^V h^q。截断与单调收敛
允许扩展值，不需要时间独立、空间独立或事件微分。

在抽象空间 (0,V) 取 f_t^*(s)=F(t)1_{s<theta(t)}，其峰值满足 F，质量为
min(VF,M)<=M，势恰为 h，故上预算确为此信息类的最大值。可测性来自 F、M
的可测性。原 Lebesgue 正面非原子，适用空间假设成立；sharp 性没有增加
semigroup、空间几何或所有 tube 的同时交叠约束。

对 f_t^* 再用 Minkowski，重排预算的 q-root 不超过它自身单支逐时 Lq norm
积分。这个比较仅对同一 F/M/V 单支成立，不能推断它优于 074 在时间上切换
三种估计分支的整体 min。作者保持了这个边界。

## 两个应用、全部分支与可积性

074 的实际分解 P_t+B_t=p_t 是正的真实事件分解；其 F_k=P_k+B_k 同时支付
两类项，不把扩大用于上估计的事件误认成互斥。其正时间项来自原模型、同一 t。
049/051 的真实二维边际保留一个固定 B 坐标，其余四个正坐标有 1/16 的质量，
这个因子只用于完整 p_t，没有赋给单个参数事件。

在每个 tube 上应用上定理于 (F_k,L_m,V_k)。对同一完整势 G，空间正并集
给出 integral_D G^q<=sum_k n_k J_k，n=(1,5,10,10,5)，和为 31。
这里无需也没有假设 tubes 互斥，先加 I 再开 q-root 的次序正确。

对整个 D，F_*=max_k F_k 是点态上界：每点至少属于一个对应 tube。
V_*=sum n_k V_k 是体积上界。每点有一个自由物理坐标落在 [15/64,49/64]，
所以正并集支付的质量是 5L_star；另一个完整正面质量界为
t/[pi(1+1024t^2)]。两者同时有效，质量上界可取它们的 min，计数为 5，
没有误用 31、B-block 大小或额外反射因子。因此 whole-D 预算也上包同一
integral_D G^q。两个合法上预算的 min 仍是上预算。

此 min 不是施加所有 tube 交叠约束后的共同最优值。逐管预算可看作分别求解
31 个信息类后正合并；whole-D 预算是另一个单独信息类。它们的抽象极值族无需
在同一个物理 D 上同时实现。

短时 F_k=P_k+B_k 的 P_k 系数严格正，故 d*t^-2<=F_k<=C(t^-2+1)。
各投影区间紧含于 (0,1)，真实 L_m=O(t^2)，whole 质量因 min 也为 O(t^2)。
于是 theta=O(t^4)。当 s<theta(t) 时 t>=const*s^(1/4)，积分 t^-2 得
h(s)=O(1+s^-1/4)。q/4=3/10<1，故 h^q 在 s=0 可积。长时全部 F_k=O(t^-6)，
其时间尾可积；有限类型的 max F 同样成立。此为全时间存在与有限性证明，
不是竞争性数值上预算。没有在积分范围中遗留短时间洞。

## 六个抽象 profile 的逐时可行性

092 补充证明的统一下系数适用于整个 0<t<=b=1/64，而不只适用于它最终
积分的 [1/1000,b]。F=P+B 中丢弃正的 P 后，剩余分母随 t 增大，故其
f_k 下系数在整个 (0,b] 有效。此处 f_k<=F_k 是“上预算函数有一个下界”，
没有把 f_k 当成实际密度的下界。

每个 L_m 的投影区间含 [1/4,3/4]，包括 singleton m=1 的等号情形；whole
投影区间也包含它。令 u=z^2、c0=1024t^2，分母
D=(1+u+c0)^2-4u 在 c0 增，在该矩形对 u 的导数<=-3/8，故最大值
为 z=1/4,t=b 时的 377/256。积分公共半长区间得
L_m,L_star >=16384/(377*pi^2)*t^2。
完整正面质量另满足 t/[pi(1+1024t^2)]>=256/(5*pi)*t^2。
所以 whole 的两个 mass-min 分支都分别下包，c_W 取它们下系数的 min 合法。
独立保存核验还以继承的 pi 有理上端检查 c_J、c_W 的准确不等式方向，
并逐一检查五个已保存投影区间对公共区间的包含；未计算任何 L 函数值。

五管分别使用 (f_k,v_k,c_J)，whole 使用 (max f_k,sum n_k v_k,c_W)。
本次运行实际采用保存的有理 m<=min(b,sqrt(vf/c))，应把其支撑明确写成

    v_tilde=c*m^2/f <= v <= V,
    fhat_t(s)=f*1_{0<t<=b}*1_{0<s<min(v_tilde,c*t^2/f)}.

这是作者所说“用有理下 cutoff 缩小实际体积”的准确函数。保存字段
actual_used_volume 正是 v_tilde；如果仍以原 v 定义 profile，所报积分下端
仍可作其下界，但不能把缩体积分误称为原 v profile 的完全精确积分。
这是一项明确化，不构成排除证书缺口。

逐时间段验证：0<t<=m 时质量=c*t^2，m<t<=b 时质量=c*m^2<=c*t^2；
窗口外为零。支撑始终<=v_tilde，峰值为 f<=F(t)，且 c*t^2<=M(t) 在窗口
每个 t 都成立。因而整个 t>0 满足约束，而非在端点采样验证。六份 profile
均分别通过准确参数、支撑体积和质量恒等式核验。

其势为 f(b-sqrt(fs/c))，0<s<v_tilde；变元 y=sqrt(fs/c) 给出

    W=2*c*f^(q-1) integral_0^m y(b-y)^q dy
     =2*c*f^(1/5)*[b^(16/5)/((q+1)(q+2))
        -b*(b-m)^(11/5)/(q+1)+(b-m)^(16/5)/(q+2)].

负项在下端使用根上端，正项使用根下端；bracket 正下端经核验后，才能
与 f^(1/5) 的正下端相乘。独立 checker 逐幂验证全部 30 个保存根区间、
根的 radicand/degree 与上述公式的绑定，以及六个完整有理积分端点。
没有重算根搜索、旧 F/L 或作者科学配置。

## 准确比较与计划账本

每个 W_k 是相应 sharp 信息类中的可行目标值，所以 J_k>=W_k；whole
同理。两侧先各自下包，再取 min：

    J_tube >=sum n_k W_k^-， J_whole>=W_*^-，
    J_new=min(J_tube,J_whole)>=min(sum n_k W_k^-,W_*^-).

此处不需要六个见证同时满足所有 tube 交叠关系；也正因未包含这种关系，
不能将 W_new^- 当成实际 integral_D G^q 的下界。

|对象|I 下端展示|B=I^(5/6) 下端展示|
|---|---:|---:|
|逐管正加权预算|3.8474459429884434e-5|2.0940020646672805e-4|
|whole-D 预算|3.081644809085581e-5|1.7404123257458573e-4|
|两预算的 min|3.081644809085581e-5|1.7404123257458573e-4|
|062 固定 D+F 充分额度|7.51537854767938e-6|5.3698116236626482e-5|

展示小数不参与判断。独立程序以062原件的完整有理数逐一核对，两个预算
的 I/B 下端均严格超过额度，min 也严格超过。B 下界由非负六次根证书
B_lower^6<=I_lower^5 得出，没有把上下端方向取反。

这里对 C/E/H 使用的是已选保守预算常数 paid，而非实际积分的下界。
定义固定计划 I_plan=paid+J_new，在未证条件 g0=1/65000 下有

    U_plan=(6/g0)^(6/5)*I_plan,
    S_plan=64*(g0/6)^(6/5)-I_plan.

因为 J_new 是完整时间 G 的上预算，不能另加051 short；far=0 只作对这份
预算最乐观的失败诊断。独立检查 target 正项取上端、I_plan 取下端后相减，
得到 S_plan 上端约 -2.3301069543176437e-5<0；U_plan 下端约
183.3035548729532>64。它们拒绝的是此固定计划预算，绝非实际 I 或 U 的下界。
更大的已证 g、重新分配、实际交叠或更强边际仍在本次排除范围之外。

## 真实运行与交付

继续采用本机 `math-l03-rational-and-interval-certificates` 的精确残差方法，
实际运行本机 Python Fraction，仅使用有理恒等式和整数幂比较。没有浮点
决定、无未经认证区间库、无 MATLAB 调用或状态清理。

独立 saved-only 组 1/3 一次执行 exit 0：229 个断言，25 个源 hash，
30 个根证书、6 个 profile、2 个分别预算比较。CPU 0.004358 秒，
外壳含 reap 墙钟 0.038502 秒。0 新参数配置、0 新根搜索、0 作者入口、
0 旧 F/L、0 旧86格、0 kernel。首失败停止策略未触发，无科学失败或重试。
早期只读展示脚本曾因把 profiles 字典当作数组而停止；它只读打印、无数学
核验或写源操作，随后修正为按字典展示，没有重新执行任何科学配置。

登记先于运行：CPU 硬上限20秒、prelimit18秒、alarm165秒、外层timeout170秒，
为每程序含保存/reap/清理总墙钟180秒留余量。`launcher_receipt.json` 记录真实
外壳收据；checker 排他输出，不能覆盖后重跑。

必要文件仅登记、源门禁、执行绑定、独立 checker、review_results、外壳收据、
本文和 manifest。所有作者和他人文件只读；未 stage/commit/push，未递归分派。
没有发现需要作者修改数值结论的缺口。封存后停止，不启动新研究。
