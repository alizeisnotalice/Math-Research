# 第31轮固定盒 builder 运行前审读

结论：最终 builder、PLAN、线性单元引理与分阶段包装器未发现阻断，可由主侧执行唯一 build。此次只读源码、依赖接口类型声明并核 SHA256；**没有导入研究模块、构造模型、评价点、调用Arb/SymPy或求LP**。search尚未实现，本结论不包括未出现的search代码。

实际冻结 hash：

- `build_cells.py`：`a1c763e2c515370f182611900987ab4b34ba3010d68277b54bca86080f74a718`。
- `run_stage.py`：`ebcad288ccd4e20144a5b93e3a3ce23bdba3052aa9d208c2eb32d9e9f8935911`。

## 精确端点与模型安全性

代码对每个Arb系数先核finite，再取 `lower()/upper()`；`exact_point`检查端点exact，按 `man_exp()` 的整数尾数及二进制指数直接构造Fraction。没有显示字符串解析、double中点或向内取整，正/负指数分支正确；再核lower<=upper。本机已装arb类型接口确有 `is_exact/man_exp/lower/upper`，此次仅静态查看stub，没有运行其方法。PLAN初稿写fmpq，代码采用同样准确的dyadic分解，这一措辞差异已通知作者同步，不是数学阻断。

两个公共模块导入后设置ctx.prec=192，再准确重建MC3/MC3T，避免导入时128位转换矩阵宽度遗留；Spec随后才构造。`checker.Spec is core.Spec`另查模块绑定，来源hash固定。每方向独立AKChecker，未跨不同mlist共享phi行或跨方向复用D_table；性能尚无实测，不能据静态审计预报300秒内完成。

12盒来自三个固定父knot的四个局部半盒，实际调用传 `(p,q,x0,x1,y0,y1)`，物理端点另保存为h(p+x)/h(q+y)，没有混用坐标。全部严格在内区离轴、原点及r=R；半盒接触节点的端点由公共phi_table弦分支处理，不把不存在的全盒高阶导数界用于节点。代码断言局部端点在[0,1]、物理分支合法。

APsi模型未乘c，整体S1模型不含奇异项；余项取准确上端，并要求其下端非负。4×4多项式先变换成16个Bernstein球再导出端点。全行用c/s非负下端−余项；signed T使用系数中点和负的半径+余项，配合z>=abs(t)。因此192×254行的列顺序及符号与 `linear_cell_lemma.md` 一致。引理的凸包/绝对余项证明成立；不能把某方向bern_min乘负t。

## 基线优先、配额与输出

程序先构12个未乘c的APsi与12个整体S1，共24模型，落盘其球、准确端点/余项及来源hash。随后对192行准确计算 `c1*(Lpsi-Epsi)+(LS1-ES1)`，判号无容差。若任一负行，终态为 `baseline_enclosure_failed`、cells为 `partial_baseline_only`、complete/baseline_passed均不放行，**不会构造1512个T模型、删行、细分或启动search**。失败只说明固定充分包围未接纳已证基线，不是基线真实负值。

只有基线全过才读取固定细化定义，构造126个准确共同h方向，各12盒，共1512T模型；规范索引、支撑和方向顺序受固定第29轮result及源hash约束。总上限12APsi+1524AS=1536。每项先记started/in-flight，模型及准确端点文件保存后才记completed；首次异常/非有限/到限停止，终态不伪报完整。完整192条行与全部1536模型完成才将cells.complete设true。

模型文件包含power/Bern显示球、准确端点、中点/半径、余项上端、Spec和模型hash；cells引用这些文件。in-flight失败项与早停部分产物不得供下一阶段求解。未来search还须核builder终态、complete/baseline_passed、模型/行hash和准确候选行判定，不能只看包装器exit0。

包装器分别以build300秒、search180秒创建专属进程组，计时包括imports/读写/启动。拒绝已有阶段证据，单次调用不重试；超时SIGKILL自有组并等待，正常退出后也清理自有残留组，保留源前后及日志hash。它不会替代build通过门槛；本轮尚未运行时，预算不能写成实际次数。

12盒仅是局部充分条件，不覆盖其余区域、轴、原点、截断线和无穷远；即使后续局部行通过，也不能宣布新全域弱型界。
