# 第49轮first静态放行

结论：**GO，仅放行固定公开50原子配置的一次实际核验**；不是数学运行通过结论。审读未运行checker、优化、编译或Git。审查采用code-review-and-quality技能，正确性优先；未因风格或通用化要求扩大本轮。

## 固定文件

- check_exact2.py SHA：`abb311c73157e35f384ee2aca6aded4772cc13ea8f70650b8499400621a80e23`。
- run.py SHA：`f85949d7f5eb7aecc7a794873e5e8e25c150e469741fc76b060b828140c17f73`。
- c2lower_1p70272.json SHA：`e2ca1616be95afa6680e3ccf157cc16fdd39026ad0f5fa79ebd32e5214779888`。
- 公开期望日志c2lower_1p70272.check2.json SHA：`4f15133b984d33430556cedc3393c5da41b56ea6f60d6a789b93801516bba356`。
- Dr固定2^40=1099511627776，lam固定1；期望分子`482167717377842991882283179477574801392641890014291050323`，分母`283174368139841919615520431313000857449005056000000000000`。展示截断17027237335/10^10不替代该分数。

## first数学路径

1. 配置以Fraction读取；点属于基本胞且质量正，完整胞质量为所有配置质量之和。周期窗口中的同坐标质量在整数前缀和中相加，不丢失重复质量。
2. 质量Ms/Dm与坐标采用共同整数尺度。`k=isqrt(floor(Ms*Dr^2/(4Dm)))-1`，k>0时ρ=k/Dr严格小于sqrt(Ms/(4Dm))，故统一半径中心正方形的平均严格大于1。没有独立坐标半径。
3. 枚举闭坐标盒选择窗口内的真实原子子集，前缀和准确计算其质量。即使有限窗口/端点范围漏掉其他见证，也只减少覆盖，不会夸大下界；本次不声称穷尽所有尺度。
4. 用[xb−ρ,xa+ρ]×[yb−ρ,ya+ρ]产生保守中心矩形，裁剪至完整基本胞；退化/空矩形丢弃。坐标在lcm(Dx,Dr)尺度下完全为整数。
5. 并集扫描线/线段树按相邻x差累计当前覆盖y长度；更新计数表示覆盖重数，不把重叠重复计面积。Aint/L²准确，最终比为A/完整胞质量。浮点只打印面积展示，不进入准确分子分母。
6. 有限类型见证、严格内缩及周期块截断/原子平滑的L1传递沿公开NOTE_lower_bound.md Lemma2和既有lower_transfer接口；本次没有重新构造Lean证明。

## wrapper路径与终态判据

wrapper将原配置复制为first.json，SHA断言后用`python -I check_exact2.py first.json --Dr 2^40`。公开程序的输出first.check2.json写在复制配置旁，未改上游配置或源码。已有对应run/log拒绝再次启动；独立进程session/PGID，单次等待180秒预算、超时SIGKILL并保存timed_out，结束清理本轮owned group。

wrapper记录source/config before/after、日志SHA和退出状态。first程序本身不对期望准确分数作assert，wrapper也不自动把SHA变化改为invalid；因此**exit0仅表示程序结束**。主侧终态接纳必须同时确认：

- source/config前后SHA均等于上述固定值；wrapper SHA仍为冻结值；没有超时。
- first.check2.json确实新产出，Dr=2^40；准确ratio_num/ratio_den与固定公开期望相等（允许以Fraction规范化比较）。
- 第二独立路径终态/准确结果另由audit审读，不把first替代second。

这些是现有产物的终态准确检查，不需要新增kernel评价或重跑。失败/超时如实保留，不转为数学反例，也不追加58原子搜索。本轮预算只有固定公开50原子的两次各180秒核验。
