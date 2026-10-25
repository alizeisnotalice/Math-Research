# C026：固定点正时间势的准确 PF 积分工具

有限验证 **PASS**。本工具从已认证 B023 的 ordinary real principal table 直接积分，在固定
\(x_0=(50,50,50,45,27,0)/50\)、唯一正时间 \(T=1/25\) 得到 Arb640 的 \(K_T\) 和 \(p_T\)。\(K_T\) 严格位于冻结准确 Green 原件减唯一一次 C010 Arb512 完整 \(H_T\) 的宽参考内。没有空间梯度、盒 oracle、根区间或 C6 新界。

## 来源、归一化和零节点路径

冻结 B023 表给出
\[
p_T(x_0)=\frac1{64\pi^6}\sum_{(s,j)}\frac{d_{s,j}}{(T^2+s)^j}.
\]
共有417非零主部项、359个正参数，最大实际阶数3。程序逐个准确检查每个最简有理 \(s>0\) 的分子、分母均为整数平方，由 `isqrt` 还原正有理 \(b\)，并核对 \(b^2=s\)。此处不能用近似浮点平方根。

这些系数已包含原群 stabilizer12 和全局负号配对；本工具对每项只使用一次，没有再乘12或2。B023 主表 SHA 为
`8446b3fe84136c005d66501a64a7978e1c1adc452bc1993acd25545ddc8ff299`。
F023 已独立完成准确 Hermite Green 全群审核，报告 SHA 为
`2851b03ea55e4ce3e95e7bbbae0b97daa9724a01f0de051aadc586255427544e`。
其准确 Γ 与 B023 原件逐字符相同，使用
\(G(x_0)=\Gamma/(128\pi^5)\)，没有新 Green tuple 遍历或全群调用。

B023 在完整汇总后才证明所有偶数阶零极点通道准确消失，奇数阶零通道的实部恒为零。程序核对保留的偶通道确实为0；这不是在原 Fourier sum 中跳过零节点。本轮 native H 另走冻结 C010 的完整含零 Hermite 路径：零节点值为 \(\arctan(0/T)=0\)，重复零节点采用准确的非零奇数阶 normalized jets，绝不以 sign(0)=0 删除导数。

## 积分递推和数学合法性

令 \(I_j(T;s)=\int_0^T(t^2+s)^{-j}\,dt\)。\(s=b^2>0\)，于是
\[
I_1(T;s)=\frac{\arctan(T/b)}b,
\qquad
I_j(T;s)=\frac{T}{2s(j-1)(T^2+s)^{j-1}}
+\frac{2j-3}{2s(j-1)}I_{j-1}(T;s).
\]
对有理项求导即得递推；或准确核对右侧导数为 \((T^2+s)^{-j}\)，且其 \(T\downarrow0\) 极限为0，便由微积分基本定理唯一确定该积分。所有359个 \(s\) 均正，有限417项在任意有限闭时间段上连续，故可以逐项积分。由此
\[
K_T(x_0)=\frac1{64\pi^6}\sum_{s,j}d_{s,j}I_j(T;s),
\qquad \frac{d}{dT}K_T(x_0)=p_T(x_0).
\]
这是一条固定空间点的一元闭式计算公式，不含空间或盒导数。

实际运行8个 FORMAL scalar checks：对 \(j=1,2,3,4\) 各做一次符号导数残差和一次零时间极限，全部准确为0。这些检查未代入任何正时间。实际表只有 \(j\le3\)；第四阶检查用于递推模板。`scalar_checks.json` 记录实际 SymPy1.14.0 结果。

## 唯一正时间的实际结果

输入预登记与来源哈希见 `registration.json`，登记 input SHA 为
`89d4a7d72b5530ea26049920be79536a13a0af874bdd7ab46ac8c13c543a5e7f`。
全部计算仅使用 \(T=1/25\)。PF 主部积分及准确有理密度求和耗时0.013183792005293071秒。密度先准确汇总为一个 fmpq 有理数，再除 \(64\pi^6\)；其完整有理原件和 Arb640 端点保存于 `pf_result.json`。同时通过一次保存表的准确全时积分代数检查得到同一 Γ，这不是新的完整 tuple Green pass。

便于阅读的数值为

|量|显示近似值|
|---|---:|
|\(K_{1/25}(x_0)\)|9.255798686659386e-6|
|\(p_{1/25}(x_0)\)|2.6583119036138666e-4|
|\(G(x_0)\)|1.7540724584574725e-5|

显示小数不用于证书。所有准确外包区间均保存完整 lower/upper integer mantissa 与 common binary exponent；没有从短 radius 展示重造端点。

随后且仅一次调用冻结 C010 H native512 binary，SHA
`429cdccca45a69c1de6c04725231297e00ee8d3951ccbb1a6141feb934f13bb8`。
fan SHA 为 `7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812`。wrapper 只传规范整数字面量并先检查 fan/source/runtime 哈希，限制旧 CLI parser 的可信输入范围；没有改动或重新编译旧源码。

|native 完整群项统计|实际|
|---|---:|
|完整逻辑项|967680|
|实际 sign-paired DD evaluations|483840|
|含零节点逻辑项|88320|
|重复节点逻辑项|201408|
|缓存准确节点值|1023|

零／重项和 C022 冻结预检统计完全吻合。native 内核1.14364秒，退出0；包括8个符号检查和本批两个实现的主要运行共1.299768666998716秒。完整 command、墙钟与退出记录在 `execution_record.json`，native 完整 dyadic 原始值在 `native_H.json`。

以冻结 Γ 的 Arb640 包围减 native H512 构造独立宽参考，PF640 \(K_T\) 严格在其内；`comparison.json` 保存完整参考端点。纯只读核验还用 Fraction 对保存 G 与 H 的准确 dyadic 端点直接作 Minkowski 差，并得到同样的严格包含。区间方向是“较窄 PF 结果在宽参考内”，作为实现交叉一致性；没有把此包含倒过来当作更强外包定理。

密度通过准确 PF 有理求和及 FORMAL 导数恒等式得到；本轮没有独立完整密度 kernel 交叉。没有从旧全域模型结果推得当前 C6 界。原 PF 表的完整密度身份依赖 B023 理论与 F023 来源审核；本轮只验证此固定点积分实现及 H 交叉，并不重新证明全部417个系数。

## 复核、资源与边界

使用已有 Python3.15.0 / python-flint0.9.0 / FLINT3.6.0 / GMP6.3.0 / MPFR4.2.2；正式 PF 运算 Arb640、独立 H 运算 Arb512，无 binary64 数学替代。旧 C010 header/runtime/source/deps SHA 均只读引用于登记和 manifest。没有编译、安装、改系统路径或复制旧 binary 到归档。

纯只读验证命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round026/integrated_pf_potential/verify_saved.py
```

该程序仅检查保存哈希、formal receipts、正有理平方根、准确端点差和 Γ／density 单个准确原件的高精度转换；不重新求 PF 积分和，不启动 native，不做全群或 tuple 遍历，不写文件。`run_once.py` 遇已有执行记录会拒绝任何再次运行。白名单不含 pycache，也不含任何旧 binary／头文件副本。

本轮20分钟有限预算内结束：8个 FORMAL cases、唯一正时间、唯一 native H、0新Gfull、0 root操作、0空间扫描、0梯度、0盒oracle、0 MATLAB、0 MC、0扩招、0提交／push。没有读取或自动承接 B024 root 值。若未来使用不同T或根 bracket，必须先另行登记；当前源码固定 `T=1/25`。这是一项固定点软件/符号闭式工具，不是 C6bound，主代理保留后续证明与集成责任。
