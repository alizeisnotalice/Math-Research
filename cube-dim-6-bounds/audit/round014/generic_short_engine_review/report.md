# F014：generic 短势、梯度与密度引擎独立审计

**在已验证 fan、规范整数输入、全部准确节点非零和已锁定 ABI 的范围内通过。** 有限 T 的必要 `f′` 对角项、非零重复的六阶 jets、copied-node 梯度、密度复乘积及联合 Euler 符号均正确。一次独立保存数据/标量复核通过，未执行 binary、完整核、完整参考或重新编译。

审核源码 SHA256 为

`247a1f6c6967ad52f34320b999818214488700651b46cb8e5feb3eac7064953b`，

已核 binary SHA256 为

`67119aa799041ed75b66f043444640a2a62556f58517cbdf3bd6f7cccd564f04`。

以下仅认证软件算法和本批证据的适用范围，不认证全域最小值或新 C6 界。F013 已冻结的整数、群映射和 ABI 分析沿用，不重复其实际结构计算。

## 1. 分支、归一化和可求导的真实核

对 `T>0,b≠0`，principal 实 arctan 满足

\[
f_T(b)=\arctan(T/b)
 ={\pi\over2}\operatorname{sgn}b-\arctan(b/T).           \tag{1}
\]

负 b 时右侧的 π/2 也取负，不能固定加正 π/2。源码 Node 构造器直接以 Arb 求 `atan(T/x)`，因此分支正确。令六节点为 `b_j=x·r_j/q_j`、`α=|det R|/∏q_j`；由 F013 已审的 Green 与截时公式相减，

\[
K_T(x)=(2\pi)^{-6}\sum_R\alpha_R D_5(f_T;b_R).          \tag{2}
\]

这是实际 `∫₀ᵀp_tdt` 在合法 generic 域的表达。代码直接算(2)，无需数值上先求两个接近的大势再相减。前因子为 `(2π)⁻⁶`，**不是** Green 的 `1/(128π⁵)`；源码第57行将值、梯度、密度一起除 `(2π)⁶`，准确匹配(2)及 Fourier 密度公式。

全部 finite full-group 节点准确非零时，存在共同开邻域使每个节点符号固定并离0正距离。`f_T` 在其各自非零节点附近解析；F013 的分离轮廓/Hermite 证明适用，即使节点重合也能对真实有限和求导。不以旧低维连续性替代该合法域论证。

若任一节点准确为0，full driver及evaluate都先触发 ZeroNode，main返回 exit3、`unsupported_zero_node`、`full_group_completed=false`，整次不输出成功结果。没有 sign(0)=0、drop cone、忽略主部或“只保留 p”的绕过。虽然 p_T 单独在零节点有定义，这个联合 generic 接口仍整体拒绝它。

## 2. 节点偏导必须包含 f′，同号节点不能归零

六节点互异时，令

\[
w_i=\prod_{j\ne i}(b_i-b_j)^{-1},\qquad c_i=f_T(b_i)w_i.
\]

直接微分 barycentric 值给

\[
\partial_{b_l}D_5
 =f'_T(b_l)w_l-\sum_{i\ne l}{c_l+c_i\over b_l-b_i},
\qquad f'_T(b)=-{T\over T^2+b^2}.                      \tag{3}
\]

源码第45行将 partial_i 初始化为 `−re_i*w_i`，其中 `re_i=T/(T²+b_i²)`；第46行再加入 pair 项，恰为(3)。这一对角项是相对 C013 sign 引擎的关键变化，没有遗漏。程序没有“全同号则值/梯度为0”的分支；有限 T 的 f_T 在同号区域仍非恒定。

重复节点时，程序用 Hermite 表和

\[
\partial_{b_j}D_5(f_T;b_0,\ldots,b_5)
 =D_6(f_T;b_0,\ldots,b_5,b_j).                         \tag{4}
\]

这是对一个带标签节点的导数；复制一次即可，物理 chain 负责共同移动多个同值标签的重数。节点排序不改变对称 DD，partial输出仍对应原标签。

核对真实 jets：写 `D=T²+b²`，

\[
{T\over D+2bh+h^2}=\sum_{k\ge0}g_kh^k,
\quad g_0=T/D,\quad
g_k=-{2bg_{k-1}+g_{k-2}\over D},\quad g_{-1}=0.
\]

因为 `∂_h f_T(b+h)=−T/[D+2bh+h²]`，归一化 Taylor jet 为

`jet_k=f_T^(k)(b)/k!=−g_(k−1)/k`，`1≤k≤6`。

第28行正是该递推，保留到6阶，足够覆盖七节点 DD6 全部重合。这里是归一化 jet，不能额外再除一次 factorial。准确相同才用 jet；准确不同而 Arb 差球含0时仍报错停止。

## 3. 密度与联合 Euler 的符号

由(1)有

\[
\partial_T f_T(b)={b\over T^2+b^2}.
\]

令 `g(b)=(T−ib)⁻¹`，则 `Im g(b)=b/(T²+b²)`。几何级数/有理差商恒等式给

\[
D_5(g;b_0,\ldots,b_5)=i^5\prod_j(T-ib_j)^{-1}
=i\prod_j(T-ib_j)^{-1}.
\]

取虚部即得

\[
D_5\!\left({b\over T^2+b^2};b_0,\ldots,b_5\right)
=\Re\prod_j(T-ib_j)^{-1}.                              \tag{5}
\]

重复情形由有理函数的 Hermite 连续延拓包含。因此 finite fan 的(5)正是实际 `p_T=∂_T K_T`，不只是与 K 有某种数值相关的辅助量。

源码各因子用 `re=T/(T²+b²)`、`im=b/(T²+b²)`，第41行以两对 Arb 实数计算复乘积。先从旧 pr/pi 同时计算 nr/ni，再写回，未混用已更新的实部；没有 binary64 复数。单个 cone 的实部可为负，正确正性来自全 fan 的实际密度，不要求每个 cone 单独正。

缩放 `K_(λT)(λx)=λ⁻⁵K_T(x)` 微分后得到

\[
x\cdot\nabla K_T+5K_T+Tp_T=0.                         \tag{6}
\]

源码累积的是 **加** `T*p_T`，与(6)一致。单个 scalar DD 也满足 `Σb_j partial_j+5D5+T p_cone=0`；本次独立 scalar 与完整点端点均核验该恒等式。

## 4. 物理 chain、群配对及整数/ABI范围

沿用 F013 已核映射 `y_j=s_jx_(perm[j])`，因此

`∂_(x_perm[j])b_m=s_jr_(m,j)/q_m`。

源码第55行先在 y 坐标求和，再按 sign/permutation映回物理x，不再次除point denominator。f_T是奇函数，五阶DD在全局 `b→−b` 下为偶；节点偏导和Jacobian各反一次号，故物理梯度配对相同。p复乘积在全局反号下共轭，其实部同样相同。alpha2 的精确乘2因此对 K、全部梯度和 p 均成立。483840计算代表967680逻辑项，不能改称少算一半数学对象。

Rat归约、排序、相等、节点符号及 `__int128` 范围结构与已审013一致。本轮新增 T 的两个正int64整数只在 Arb 中形成有理值，不引入新的未检查整数乘积。准确不同节点的差球含0、零节点或超范围值均不能成为成功证书。

匹配头archive、FLINT3.6.0/GMP6.3.0/MPFR4.2.2 runtime和本014 binary 的实际哈希均重新核对。源码使用同一已验 generic C ABI 和Arb结构；没有引入新的 native Acb 结构。编译时大小检查与运行时FLINT版本检查仍在。此结论仅限已锁定环境，不是跨平台或任意库升级认证。

## 5. 一次独立保存数据与六配置复核

本目录 `verify_saved.py` 不导入 C 的求值函数、不运行其 binary或脚本，使用同六个已登记配置，T固定4/125。五个非零配置包含互异混号、非零重复、六重正节点、三重/二重混号及重复有理节点，第六个含零配置只核拒绝。

为区别于 C 的 Newton DD/Acb参考，本次使用：

- SymPy 准确有理微分构造正阶 jets；零阶以1536bit Arb的 `(π/2)sign−atan(b/T)` 表达。
- 准确 Hermite 线性系统构造足够阶插值多项式，再由 monomial DD 的 complete homogeneous polynomial 恒等式得到值及六偏导，避免复用 C 的差商递推。
- 对有理函数 `b/(T²+b²)`，用另一六条件 Hermite 插值的最高阶系数求准确 DD5，并与纯 Fraction 的复有理乘积实部准确比较。

五组 native512 完整端点均严格包含本次1536bit值/全部偏导外包；五个 density 有理结果逐项相等且位于native p-cone完整区间内。原存档640bit参考的全部 dyadic 端点也被 native 包含。每组 scalar Euler 同时通过高精度参考和 native Fraction 区间检查。没有添加第七组或新T，也没有为了零拒绝重启binary。

完整点仍是同一旧中心。脚本重新核对四份来源的点分子/分母、完整运行命令、时间、fan/源码/binary哈希和群项数，并由大整数重建 `G=Q/π⁵`、`∂₂G=Q₂/π⁵`。随后读取已存 H32 及其 AD 梯度的完整384bit证书，构造

`Kref=G−H32`，`Dref=∂₂G−∂₂H32`。

此次方向与013不同：

\[
\text{native K interval}\subset\operatorname{int}(Kref),\qquad
\text{native gradient}_2\subset\operatorname{int}(Dref).             \tag{7}
\]

由于 H32/其AD输入是较宽的旧区间，(7)是严格的一致性交叉检查，**不能单凭它把旧宽证书缩成 native 新证书**。native严格外包的依据仍是本报告审核的实际数学算法、准确输入和Arb区间运算。报告没有把(7)说成 native包含一个更精确的完整点参考，也没有从展示数重造区间。

完整点的 p 没有独立第二次全群参考；它的依据是已证公式(5)、源码及五组准确标量锚，另由Euler交叉检查。其他五个 full 梯度也没有各自新fmpq完整参考。脚本用 Fraction 对完整 K、六梯度、p 的dyadic端点重建(6)，准确含0；native存档Euler也含0，且完整点p下端严格为正。这些检查没有被夸大成全域或连续箱验证。

实际一次运行：PID32487，Python3.15.0 / SymPy1.14.0，脚本计时 `0.07796120797866024` 秒，工具墙时 `0.135765291` 秒，exit0。零fullkernel、零binary launch、零编译、零MC/MATLAB。

## 6. Trusted-input 与 release 边界

014继承013的parser限制：`fan()` 在读首字段前没有验证文件流，空/不存在文件可能读取未初始化n；`std::stoll`不检查尾部文本；21cells/数值范围检查也不证明fan的几何合法性。外层必须检查文件可读及准确哈希、完整token格式、规范整数参数，并在运行后检查exit0、正确status、完整JSON和群项数。任何部分输出或exit2/3均应拒绝。

“已验证 fan”必须包括：全空间覆盖、cone内symbol线性、准确正q与det、以及使用同一canonical rates进行signed-permutation复制的真实对称性。任意一般frame不自动具备此对称性。冻结body32 fan在此域内；若改用 augmented frame，需引用其独立fan/rate/对称性证明，不可只改一个文件名便继承旧模型数值。本报告没有重跑或单独数值验证augmented模型。

本报告审核的原始014源码保持冻结，未代修parser。上述trusted wrapper条件足以隔离当前有效批次与已知错误输入风险；它们不使已冻结结果失效，也不能省略。精度512只是固定实现选择；若球太宽或无法分离节点，应停止或在另行登记批中改精度，不允许丢项。

软件提供单点 K/gradient/p 外包，既不是全cube下包，也不是能自动处理全部零节点的全域求值器。其一次3.882秒性能记录不构成所有点吞吐保证。当前 C6 上界64未变。

## 7. 冻结交付

本目录只保留登记、独立纯复核脚本、结果、报告及manifest。完整source hashes在results中；本次文件与来源哈希在manifest。F013/C013/C014及其他目录均未改写，没有新缓存、后台进程或git操作。
