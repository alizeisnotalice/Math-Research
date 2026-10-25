# C013：非零 Fourier 节点的严格 Arb512 Green 与六分量梯度引擎

结论：一次编译版本通过本批全部验收。5组非零标量锚的值/节点偏导与精确fmpq Hermite-DD6一致；第6组含零节点，严格返回unsupported。唯一既有完整generic点上，G被旧完整fmpq准确Q/π5证书核实；x₂梯度被另一个全群fmpq DD6准确证书核实；Euler恒等式也通过。

原生G+全部6个梯度外墙1.208953708秒，内核1.19402秒。旧G单值fmpq记录23.328086167秒；独立单分量梯度fmpq本次42.251118166秒。记录之比分别19.30、34.95，但工作量和算法不同，不能将这些比值冒称同工作量或所有点的统一吞吐保证。

总full evaluator executions为2（native一次、独立Python一次），均同一旧空间点；没有新点、网格、MC、MATLAB、下界工作、扩招、stage/commit/push。现已停止并冻结013，不启动014。没有新数学上界或C6突破。

## 1. 模型和严格合法域

本工具回到原六维等权body32 **Cauchy** frame，ψ(ξ)=Σ_v|v·ξ|，V={(1,s):s∈{±1}⁵}。它不是α=3/2替代模型，也不是H_T梯度。

对现有准确21-cell fan、全部坐标置换和符号群，令α_R=det_R/∏q_j，b_j=x·r_j/q_j，使用已有公式

\[
 G(x)={1\over128\pi^5}\sum_R\alpha_R
 D_5(\operatorname{sgn};b_0,\ldots,b_5).
\]

**合法域要求全部full-group Fourier节点都准确非零。** 任何一个节点的准确整数分子为0，整个求值立即以exit3和status=unsupported_zero_node终止，不输出G或梯度。没有skip零节点，没有sign(0)=0约定，也没有用主部抵消替代这个门控。unsupported不等于该点真实G必然发散；零节点及高阶principal parts仍由旧路径负责。

非零重复节点允许。对任意全部节点离0正距离的点，有限个线性form在某个开邻域的符号全固定；sign在每个节点附近的1阶及以上jet为0。因此非零重复节点的Hermite扩展和节点导数合法。这也可从atan(b/ε)在非零节点的全部正阶jets趋0直接看出；这里没有把旧零节点主部论证迁移到非零重合节点。

程序仅接受正int64分母、非LLONG_MIN的int64点分子，fan有21cells，q∈[1,10^6]、ray分量绝对值≤10^6、det∈[1,10^9]。6项点积先在__int128累加，再逐一检查节点分子与den·q能否严格放入int64，不能放入就exit2。约分、排序、重合与符号判定全用整数/__int128；不round任何超范围输入继承证书。本批没有运行B010超范围point2。

## 2. 值和节点偏导算法

当六节点不同且非零时，记

\[
 w_i={\operatorname{sgn}(b_i)\over\prod_{k\ne i}(b_i-b_k)}.
\]

值为Σw_i。由于各节点局部符号固定，f'_i=0，直接对barycentric公式求导得

\[
 {\partial D_5\over\partial b_l}
 =-\sum_{i\ne l}{w_l+w_i\over b_l-b_i}.
\]

程序按无序pair计算一次差，分别对两节点偏导加相反贡献。它没有在热路径中重复6次整张DD表，且不需要atan或H_T计算。

若存在非零重复，则值用Hermite Newton表：端点准确相等的k≥1差商为0；否则以Arb差商递推。每个**带标签的原节点**偏导使用

\[
 {\partial D_5(f;b_0,\ldots,b_5)\over\partial b_j}
 =D_6(f;b_0,\ldots,b_5,b_j).
\]

对多重节点，移动一个原标签只复制一次；物理坐标同时移动多个标签时，再由链式求和提供重数。所有同号且非零节点时，函数局部为常数，值和所有偏导准确为0；这个分支在逐一检查六节点非零之后才允许触发。

标量锚覆盖：未排序的六个不同节点；混号重复节点；六重正节点；混号三重/二重；重复有理节点；一组含0的拒绝测试。native不同节点走barycentric偏导，独立Python参考始终走精确fmpq Hermite DD6，所以并非只用同一浮点实现互相比对。5组支持点的值及各6偏导均准确包含对应有理数，第6组无数值输出。

## 3. 原坐标链、全局负号和归一化

每个置换/符号项写y_j=s_j x_(perm[j])，因此

\[
 \partial_{x_k}b_m
 ={s_j r_{m,j}\over q_m},\qquad\operatorname{perm}[j]=k.
\]

这里对的是物理坐标x，而不是输入整数分子；chain中没有再除point denominator。程序先以partial_m/q_m合成y坐标梯度，再按s_j和perm映回六个原坐标。

程序**没有canonical整体负反射**。它只采用已证明的全局符号配对：D5(sign;−b)=D5(sign;b)。两项的节点偏导互为负，而线性form的Jacobian也互为负，故最终原坐标梯度相等；值和梯度都可用weight严格乘2。720×32×21次实际计算代表全部967680逻辑群项，没有遗漏permutation或cone。

重复节点fallback中的节点排序不改变标签偏导：对原第j节点添加其值后求对称DD6，输出仍放回原j槽。零ray chain系数产生的准确零贡献可不乘加；这与删除零Fourier节点完全不同，每个b的非零门控始终先执行。

G及梯度都除128π5。所有模型运算为Arb512；binary64仅用于chrono墙钟输出。输出均由arb_get_interval_fmpz_2exp给出完整整数端点与公共2进指数，不从短radius展示重造证书。

## 4. 唯一既有整点与独立精确证书

使用旧A007中心，分母D=1000000007，分子

`(1000000007,998836847,892738545,807464540,467291939,96968636)`。

native一次完整求值：逻辑项967680，配对计算483840，所有节点非零，full-group重复节点为0。特殊重数只由标量锚验证，未假称此整点本身测试了退化情况。

G的准确参考直接读取旧完整fmpq记录的Q/π5，不重算旧G基线。将准确有理Q除1024bit向外π5区间后，被native512完整区间严格包含。

独立x₂梯度参考则新写为Python fmpq：遍历720×64×21全部967680项，**不使用全局负号配对**，每个原节点偏导都用copied-node DD6，再乘原ray/q和实际sign/permutation chain。它没有调用native barycentric公式，也没有用Euler倒推出该分量。

参考精确识别646646个全同号非零项的导数为0；其余项实际调用1471024次精确DD6。最后得到准确有理Q₂/π5，完整有理分子/分母保存在gradient_exact_reference.json（约140KB）。1024bit外包被native的x₂分量严格包含。该步骤耗42.2511秒，未触及600秒硬限。

Euler还独立复核两次：native累积x·∇G+5G的完整区间含0；Python仅用存档的七个完整native区间和准确物理x重算，仍含0。它是额外身份检查，未被当成独立梯度证书的替代。

便于阅读的梯度约为

\[
 (-3.8229870,-3.8050647,-2.4369462,-1.5461905,-0.6999007,+1.2794268)\,10^{-5}.
\]

这些展示数不是证书端点。native_result.json保留G、六梯度和Euler的完整dyadic端点；validation_results.json只把短显示明确标为display_only。

## 5. 构建、运行和复核

复用冻结C010的13个匹配头archive和现有wheel runtime：FLINT3.6.0、GMP6.3.0、MPFR4.2.2，python-flint0.9.0。只解包到本013 scratch，未修改010。编译期检查ulong/fmpz各8字节、arb_struct48字节，运行期检查flint_version3.6.0。源码C++只编译一个版本，无实现bugfix重编译。

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/upper/round013/generic_green_engine/build_engine.py
```

脚本验证旧匹配头archive/runtime SHA，编译-O3 -std=c++20；仅修补本轮binary的构建机LC_ID到@rpath，不修改系统PATH、dylib或环境。完整命令和依赖SHA在build_record.json。

现成CLI：

```sh
cube-dim-6-bounds/upper/round013/generic_green_engine/scratch/generic_green full \
  /Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round003/direct_kernel/fan_integer_input.txt \
  1000000007 1000000007 998836847 892738545 807464540 467291939 96968636
```

这是复现用法，不表示又执行一次。本批的两个full executions已经完成并停止。run_native.py是有门控的旧点单次验证driver；exact_reference.py是独立整群参考；finalize_validation.py只读存档作最终交叉，不调用kernel。

## 6. 证据、限度和冻结状态

实际：一次编译exit0；6 scalar cases通过；native整点exit0；独立exact整群exit0；最终存档复核PASS（工具chunk25d000）。没有sign或endpoint不一致，没有放宽端点后重试，没有第二个空间点。独立reference的wall-limit guard在锚校验后从1000秒收紧为600秒，sign_dd数学函数未改；对应前后源hash在锚记录及manifest明确保留。

原生G+六梯度1.209秒与旧G单值23.328秒、Python单分量42.251秒不是同工作量比较。它支持当前点的诊断可用性，但本批没有测全组重复节点点、不同near-zero条件数点或其它空间点的速度。工具不输出连续箱证书，也不覆盖零节点、principal parts或H_T/K_T导数。

输入fan SHA为7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812；旧完整G证书SHA为d31cc5ba69a442e10d81ff2f32094152cf8e974a02ace68768901133f7efda60。算法来源为既有fan及主clone upper/round002/kernel/002_kernel_interface.md中的(B)，非零重复扩展与节点导数在本文单独说明。C++源码、准确梯度参考、依赖、binary及全部必要文件hash见manifest.json。

本轮仅完成新实现验证，不登记数学进展或任何C6界。无存活计算句柄，013交付冻结，等待另行明确有限登记任务。
