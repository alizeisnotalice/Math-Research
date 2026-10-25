# 第七轮 E：唯一普通晶格候选的独立理想并集数域认证

状态：**理想有限盒并集准确认证通过；新 C6 下界仍 pending 主代理的整数内缩/有限 L¹ 与 F 独立复核。** 本件不提前更新账本。

只认证主指定的一个有理步长，没有重跑 D 的 3 个方向，也没有调用 D 的双数并集、导数或答案系数。D 的 `positive_candidate.json` 只读取字节计算 provenance SHA-256，未解析其数域答案。

## 候选与准确结果

令正实根 u⁶=2、v⁶=3，取单位质量普通矩形晶格，spacing 为

\[
a=\left(1+\frac{u^4}{2},\quad
\frac{1000001}{1000000}\frac{u^3+u^4v^3}{4},\quad
\frac{1+u}{2},\frac{1+u}{2},\frac{1+u}{2},\frac{1+u}{2}\right).
\]

完整理想 witness ratio R 是 Q[u,v]/(u⁶−2,v⁶−3) 中的一个元素，其 **36 个准确有理系数**全部记录于 `candidate_results.json` 的 `coefficients`。该文件同时记录完整有理包围 `[R_lower,R_upper]`，而不是用显示小数回造端点。

R 的显示值为

\[
R=2.06692413818151301542586764203052093702041290108152\ldots
\]

准确 Fraction 比较全部通过：

- R_lower 严格大于主 ledger 读取时的完整旧 `lower.exact_interval[1]`。
- R_lower 严格大于 phase frozen 的完整旧上端。
- R_lower 严格大于本次待验证目标 `2.06692413818=103346206909/50000000000`。

相对完整旧上端的认证余量约为 `3.12170632105328089358e-7`；相对目标的认证余量约为 `1.51301542586764203052e-12`。上下包围宽度约 `4.27474662058e-74`。这些小数只为阅读，严格比较使用结果文件中的完整分数。

## 独立构造与精确算术

程序仅复用已审 `lower/round001_rectangular_certificate.py` 的 36 系数加、负、乘、根因子分解与符号包围。此次实际读取的是 D 独立 clone 内该文件；其 SHA-256 为

`9072275d47bd66c3ba9ec7de0c435efdba0785c39a02b122ea1ebb960403ec82`。

已只读核对主 clone 同名文件字节 hash 相同。没有导入或调用其旧 baseline `main`。本进程把正根隔离分母设为 10⁸⁰，逐项准确验证 u_lower⁶≤2<u_upper⁶、v_lower⁶≤3<v_upper⁶；全部符号决策由正实嵌入的有理包围完成，区间不能分离时抛错，不采用浮点容差。

重新遍历全部 4⁶=4096 个长度向量 ℓ∈{1,2,3,4}⁶。选定连续源块的中央起点为奇长 q_h=−⌊ℓ_h/2⌋、偶长 q_h=1−ℓ_h/2。每块质量 m=∏ℓ_h，共同半径为 r=m^(1/6)/2；所有这些根都在上述数域中，并逐个缓存核验 `(2r)^6=m`。

每个坐标的准确中心区间为

\[
[\max(0,(q_h+\ell_h-1)a_h-r),\quad
\min(a_h/2,q_ha_h+r)].
\]

仅保留正宽度盒。对每个非空盒，再独立核对两个中心端点与源块两个端点的全部距离均≤r；坐标凸性保证所选块中的全部原子均被同一个中心 cube 包含。每个盒也核验包含于真正半胞元。

实际得到：127 个非空 raw 盒、117 个不同几何盒、96 个包含极大盒，共 3048 次源端点覆盖检查。固定候选下的准确包含裁剪不会改变并集。这里计算的是这个完整所列有限 witness 族的理想并集，不是完整最大函数超水平集的等式或上界。

静态并集程序在每一维用全部准确数域端点切片，条带内活跃投影固定；末维合并区间，递归准确相加/相乘。这是新的静态候选求值，没有使用 D 的 `(F0,F1)` 或其 endpoint-rank 双数递归。

因为当前源是普通单位质量晶格，完整周期质量为 1；六个坐标反射的 64 个半胞元内部互不相交。因此 `R=64*union_volume`。这与 D 使用完整一轴周期后的 `16*union` 几何规范不同，但表示同一普通晶格的周期比值。

## 实际运行与资源

唯一实际候选运行 PID **18010**，exec session **39645**，exit **0**，耗时 **7.368186959 秒**，访问 **3302** 个缓存并集状态。上限是 1 个候选、10⁶ 个并集状态、1200 秒；没有达到上限。进程已结束，没有 MATLAB 调用、子代理、Git 暂存、提交或推送。

所有输入文件 hash、根隔离整数、完整旧门槛、系数、完整上下分数、保留长度向量、严格余量与运行结果均保存于 `candidate_results.json`。本轮源码 hash 为

`9534b7626341d95c59b24618a7e703b20ffcbc66e420faa32d4218483971ddfd`。

结果文件 SHA-256 为

`b5661c0eb23c264286085053b3c59ab9f992918dcb1dc4a79fdb01e1ddcb6b81`。

## 复现与尚待闭合项

从 E clone 根目录复现本次唯一步长：

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  cube-dim-6-bounds/audit/lower007/candidate_number_field/certify_candidate.py \
  --arithmetic /Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/lower/round001_rectangular_certificate.py \
  --main-project /Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds \
  --candidate-provenance /Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-astra/cube-dim-6-bounds/lower/round007/shape_jets/positive_candidate.json
```

归档到主 clone 后，省略前两个参数也可使用项目内默认路径。provenance 参数仅供输入追溯，数学参数固定在源代码中。代码没有读取 D 的导数或答案系数。

主代理独立执行有理内缩、正 cube 余量、明确有限 K 与一般 L¹ witness；F 独立复核理论、数域系数、准确比较与传递。只有这些路径全部通过，才可将本件 R 的完整区间登记为新动态门槛，并宣布完整 C6 新下界。
