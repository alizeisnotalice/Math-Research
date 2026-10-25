# audit_I028：I027三点负矩cap资格复核

固定句柄`audit_I028_fixed_001`，实际PID由--run在计算开始前打印并保存结果。一次有限三点批次120秒软预算；不重跑022/024几何、旧Gamma公式、旧g或SIP，不提交。

锁027 manifest `27f69351a186acfe18cf3071943933350a88de912d482e60a3dd41e60b96c867`，核五件白名单与输入hash；读取024已锁Zflag准确有理区间。只核gamma=1/5、3/14、1/4。

解析先审ER d^-gamma=Ztilt^-1 E[T^(5gamma-beta)detA^-gamma]、Harris反向单调和非负截断、0<=k<beta适用区间、qopt对gamma独立、16净幂5gamma/3。区分Ztilt与Zflag。

不同数值路径：保留未合并三块16尺度因子，Ztilt^-1用16 beta Gamma(beta)；负矩用未移位Gamma(gamma/beta)/(beta Gamma(gamma))；正矩用Gamma(-k/beta)/(beta Gamma(-k))，区别作者正参数Gamma比。Arb512bit逐端点包围，并以Fraction验证落在作者三点完整有理区间内及各严格cap；只读取同指数旧Q完整球作比较，不重算其公式。底层Gamma/Arb共享明确。

可携入口--source-root、--model-root；默认只锁输入，--run才该PID批次。必要白名单登记、代码、结果、报告；新矩cap不是I或C5界，不能移植旧majorant族排除而省略新资格审查。
