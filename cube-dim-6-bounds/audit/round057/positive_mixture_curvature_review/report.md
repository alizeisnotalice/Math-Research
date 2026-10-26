# F057 独立审计：正混合曲率与整盒接口

结论：**ACCEPT_THEORY_AND_CONDITIONAL_BOX_INTERFACE**。冻结版 §2 的固定
fiber 正混合论证、§4 的共享分母极小化及整盒 32 顶点接口均可接受，未发现需要
阻断的数学或实现错误。主线程保留最终整合裁决。没有新的真实 g 证书、覆盖结果
或 C6 上界。

## 独立证明结论

`proof.md` 独立展开了关键缺口审查。底测度 dt dw 的质量可以无限；真正需要且
已经证明的是其正核积分严格正且有限。先作 y=R_M x+w，再在固定 w、固定全空间
上微分，不存在遗漏的移动积分边界或 Jacobian 项。局部
exp(||Rh||_1/a) 倍的固定核支配同时支付一、二阶微分和无限时间尾部。按 F(x)
归一化后，二阶对数导数等于 log-kernel 曲率均值加 score 方差。

因此，原 unit-weight body32 上 F=H_a 或 B_(a,T)、a>0 的新强结论成立：

    Hess log F >= -I6/(16a^2),
    log F(c+u) >= log F(c)+grad log F(c).u-||u||^2/(32a^2).

R=A^T/32 是先固定的同一右逆，所谓 Loewner 最优仅限这个固定右逆证明族。
六 basis 的 Q6=I5/4-J5/32、可选上曲率 25/4 和曲率成本比 39/5 也成立。

固定 Q 下 exponent 的凹性支付全盒。对每个角点，f=F(c) 同时出现在 prefactor
和分母；其准确区间极小点是 clamp(b,[L,U])。这些结论不依赖 chamber、投影
分离、S3 距离，零/重复投影均在定理范围内。不能用 full G 或其他模型的 jet
替换，也不能把 a 固定后使用空间齐次性。

## 实际执行及新证据

预登记八配置，一次执行全部通过；未运行作者 `check_fixed.py`，没有重复作者
七配置。使用本机 Python 3.15.0、SymPy 1.14.0、python-flint 0.9.0；CPU
0.227374 秒、wall 0.302799 秒，低于 30/900 秒限额。

| 项目 | 实际执行及结果 |
| --- | --- |
| 新矩阵配置 | 改用 determinant -32 的 star basis，精确核对固定 fiber 仿射恒等式、AK=0，以及一个非零 nullspace 扰动的右逆范数分解 |
| 新符号配置 | 三个抽象正质量分量的 quotient/variance 恒等式；无限混合的合法性另由解析证明承担 |
| 新 frame32 band 盒 | B_{3/2,5}、非对称且有退化边的盒，三个 clamp 分支均出现，最坏角点由内部 f=3 给出，严格 1.103<3/e<1.104 |
| basis6 路径 | 同一新盒的另一固定 quadratic，正确返回 UNRESOLVED；四种指数均独立核对 |
| 向外舍入 | 对八个不同角点指数，以 200 阶精确有理 Taylor 加几何尾界独立检查 64-bit Arb 包络；全部包含严格参考区间 |
| 精确阈值 | 点盒下界恰为 2；target=2 通过，target=2+10^-19 返回 UNRESOLVED |
| 拒绝边界 | 八种 inexact/bool、八种 model/window/provenance、四种 shape/geometry 错误，全部拒绝 |

实际 source `certify` 调用 **24** 次：4 次正常接口输入，20 次错误输入拒绝。
其中正常接口的 2 个 PASS 与 2 个 UNRESOLVED 都是预期行为。32 个顶点中的重复
顶点只计为配置内部算术，不是 production points。

Production G/H/p/gradient/Hessian/kernel、full-fan、old-formal、MATLAB 调用
全部 **0**。所有 synthetic jet 都明确标为接口输入，不能作为真实 body32
数值证据。定理是解析证明；有限检查只核对代数与实现。

## 接口的准确承诺

输入端强制 root metadata、center、box、F/gradient endpoints、target 的规定
结构与精确有理解析，并核对 jet 的模型、归一化、window 和中心匹配。window
采用字典相等比较，等值而不同写法的有理字符串可能保守拒绝。`source_proof`
只检查非空字符串；核验该外部 jet 的证明是后续调用方的独立责任。

输出的 `tangent_bound_upper` 是 robust tangent minimum 的上端点，**不是**
原 F 在整盒上 infimum 的上界。PASS 仅表示：若输入 jet 的外部证明成立，则
这一整个盒的 F 下界通过目标。尚未组成 root cover，不接受把它升级为 g
全局证书。真实 jet pilot 仍需另行固定 cutoff/window 和中心，再支付覆盖缝隙。

## 冻结范围与复现

审计冻结作者 theorem SHA-256 为
`68fa917e01548c25258612e13d3f976109c55ea29ccf3456d2e83f71d995fb8a`；
checker SHA-256 为
`891a16137bbb483099d577e1ef7199e7a632f141d2517f2540713a19ad27c279`。
所有 source 文件的前后 hash 一致，完整清单见 `input_manifest.json`。

必要文件为 `registration.md`、`input_manifest.json`、`check_independent.py`、
`results.json`、`proof.md`、本报告及 `manifest.json`。仅写本审计目录，没有
stage/commit/push，没有修改作者或往轮材料。

复现原运行命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B check_independent.py
```

脚本会拒绝覆盖已有 `results.json`，须在没有该输出的干净副本中复现，并保留
manifest 指向的冻结 source。当前批次已封存停止，无追加试验授权。

封存时的行政说明：登记中的 “200 terms” 在运行前固定的代码中具体实现为
Taylor 到 200 阶，含常数共 201 项，随后加严格几何尾界；结果文件按阶数准确
记载为 `Taylor_terms_through_degree=200`。这只是登记用词澄清，没有新增配置
或追加计算。

用户在八配置运行结束后新增单程序墙钟上限：后续每个程序统一 timeout
180 秒，整体研究和思考不受该上限约束。超时只终止并 wait 本人启动的进程，
保存已得检查点，不自动重跑。原先验登记和已结束的 0.303 秒批次不追溯修改；
接收约束时没有存活计算进程。行政封存程序亦设置 180 秒 alarm。
