# C016：含零与重节点的正时间带 Arb512 引擎

本批完成软件有限验证：在唯一登记的旧中心，直接计算
\(B=H_{1/250}-H_{6/125}=\int_{1/250}^{6/125}p_t\,dt\)、全部六个空间导数及两个端点热核。完整群值严格位于旧完整 Arb384 差证书中；六组标量差商与全部六个节点偏导均通过独立 Acb640 核验。没有新的空间点、数学上界或全域覆盖结论。

## 合法域与算法

固定既有 body32 Cauchy 模型和准确 fan，归一化为
\[
B(x)=(2\pi)^{-6}\sum_C\alpha_C[b_0,\ldots,b_5]f,\qquad
f(b)=\arctan(b/a)-\arctan(b/T),\quad0<a<T.
\]
此处 \(b_j=(r_j\cdot y)/q_j\)，\(\alpha_C=\det(r_0,\ldots,r_5)/\prod q_j\)；\(y_j=s_jx_{\sigma(j)}\)。21 个 cells、720 个排列、64 个符号总计 967680 项。只用准确的全局负号配对：实际求值 483840 项，权重 2；没有对节点另外作 canonical 负反射。

实数 \(b\) 上 \(f\) 在零点也解析。程序使用等价且少消减的主支表达式
\[
f(b)=\arctan\frac{(T-a)b}{aT+b^2}.
\]
分母正，原两角之差位于 \((-\pi/2,\pi/2)\)，所以没有支路跳跃。其一阶导数为 \(a/(a^2+b^2)-T/(T^2+b^2)\)。若 \(g_u(z)=u/(u^2+(b+z)^2)=\sum c_kz^k\)，则
\(c_0=u/(u^2+b^2)\)、\(c_k=-(2bc_{k-1}+c_{k-2})/(u^2+b^2)\)，缺失的 \(c_{-1}=0\)。程序使用 \(f^{(k)}(b)/k!=(c_{k-1,a}-c_{k-1,T})/k\)，直到六阶，零点不例外。

任一重复节点均采用排序 Hermite 递推，节点偏导为复制相应节点的 DD6。不同节点采用重心公式；若 \(w_i=\prod_{j\ne i}(b_i-b_j)^{-1}\)、\(c_i=f(b_i)w_i\)，则
\[
\partial_{b_l}DD5=f'(b_l)w_l-\sum_{i\ne l}\frac{c_l+c_i}{b_l-b_i}.
\]
这和复制节点 DD6 相等。节点偏导依次除 \(q_j\)、乘准确 ray，再乘原符号并按排列送回物理坐标。有限解析和可逐项求导；此程序没有 principal-part 或 sign(0) 约定，没有删除任何零项。

端点热核每 cone 使用两份 Arb 实虚乘法计算
\(\Re\prod_{j=0}^5(u-i b_j)^{-1}\)。独立标量参考核对了
\([b_0,\ldots,b_5]\{b/(u^2+b^2)\}=\Re\prod_j(u-i b_j)^{-1}\)，包括 Hermite 重节点。它来自有理 resolvent 差商公式，DD5 的 \(i^5\) 与取虚部恰好给实部。

所有数学运算均为 Arb512；二进制浮点仅用于墙钟。输入为规范十进制有符号 int64 字符串，分母正；点、ray 点积与分母乘积经 int128 中间量和 int64 范围门控。stream/open/count/额外 token 均检查，修复了旧工具数字前缀和空 fan 读取风险。DLL 版本/结构大小也检查。对非重的准确节点若512位球差不能排除零，整体退出，绝不跳项；故接口允许所有零/重复类型，但不保证任意 int64 输入都能在固定精度成功。可信 wrapper 还必须验证 fan SHA；任意用户提供的21-cell文件不自动成为真实 fan。

## 登记与实际运行

输入登记 SHA 为 `902d934c53248c056b22e0c7d067d1be90d962503fa03e47820629b1501cc3b7`。先编译一次，再运行六组标量配置，均固定 \(a=1/250,T=6/125\)：全零；对称且两重零；四重零；非零互异；六重非零；两簇非零重复。原始输出、独立参考和全部完整 dyadic 端点见 `anchors_native.json` 与 `anchors_validation.json`。六个 DD5、36 个复制节点 DD6 均通过；12 个端点密度用准确 Fraction Hermite 与另一条准确复数乘积相等，并被 native 区间包住。全零六个偏导准确为零。解析和 parser selftests 均通过。

唯一完整群输入为旧中心
\[
x=(1000000007,998836847,892738545,807464540,467291939,96968636)/1000000007.
\]
实际 967680 逻辑项／483840 配对求值，16931 个准确节点缓存；此旧中心的重节点项和含零项计数均为 0，故完整点不额外检验零节点的群级链式汇总。零/重节点支持由解析构造和上述标量锚检验。本批没有增加特殊完整点。

内核 4.35577 秒，总墙钟 4.363948082987918 秒，退出 0；编译一次，bugfix 0，完整群执行 1。B 约为 `1.4864781404334374e-5`，仅供阅读。**全部证书以 `native_result.json` 的完整整数 mantissa 与共同二进制 exponent 为准，绝不从短 radius 展示重造。**

旧 A007 的完整 \(H_{1/250}\) Arb384 减旧 B011 的完整 \(H_{6/125}\) Arb384 作为独立宽参考，准确来源/点/分母/项数已核对；native512 band 严格在此宽参考内。这个包含方向检验实现一致性，不把宽参考包含方向反过来声称更精确认证。

两个端点密度下端均正。下式由 native Arb 运算和另一次保存端点的 Fraction 区间重建分别包含 0：
\[
x\cdot\nabla B+5B+T p_T-a p_a=0.
\]
原因是同时缩放 \((x,a,T)\) 时 B 为 -5 次齐次，且 \(\partial_TB=p_T,\partial_aB=-p_a\)。

**验证缺口：没有现存 \(T=.048\) 的完整群梯度原件，故完整 band 的任何梯度分量均未做独立全群差证书交叉。**依本任务明确允许的验收退让，本批保留 `full_gradient_crosschecked=false`；当前证据为全部标量 DD6、明确 chain 公式和 Euler 恒等式。Euler 不能单独证明六个分量全部正确。本批未进行新的 full gradient 参考运行，也不宣称该缺口已闭合。

## 依赖、复核与冻结

FLINT3.6.0 / python-flint0.9.0 / GMP6.3.0 / MPFR4.2.2。只读复用 C010 匹配的13个头文件归档和 wheel dylib；没有安装或修改系统路径。准确库、头文件、source、binary、输入哈希与编译命令分别见 `build_record.json`、`inputs.json` 和 `manifest.json`。C010/013/014 冻结源码与 manifest 的当前哈希纳入外部白名单；本批没有写旧目录。

`build_engine.py` 可在本目录 scratch 构建；`check_anchors.py` 重新做六组标量参考；`validate_full.py` 遇已有执行记录拒绝重跑。冻结后常规核验只运行：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round016/band_engine/verify_saved.py
```

该核验只读保存文件与旧来源，检查白名单哈希、标量包含、完整 band 的宽证书包含和准确区间 Euler；0 次核函数求值、0 编译、0 binary launch、0 写文件。scratch（含 binary、解压头文件）与 pycache 不在归档白名单，不复制 binary 到归档。没有 MATLAB、MC、下界计算、扩招、提交或 push。

本交付是未来 Taylor cover 的候选数学计算工具；尚无盒区间输入、二阶导数界或全域覆盖证书。没有 C6 新界；主协调者保留最终证明、release 和集成责任。
