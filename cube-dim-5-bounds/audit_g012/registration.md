# audit_g012：Astra004全域g候选的有限独立审计

锁定作者manifest SHA256=92398b5b1d5b5e26552fbaadbab4e3627d34a02fa26d9f6ce1f65e5a5d832da4。只审eta=1/8、a=333/500给出的最强g>=7/156250，不重跑旧g3，不搜索新参数。

一次固定批次，120秒软预算：

1. 用Fraction准确验证Young二次minorant的AL、BU端点及cL（Machin有理pi界），独立组装空间多项式。
2. 把完整稀疏face幂系数提至公共整数分母，直接使用tensor闭式公式累加全部2401个Bernstein系数；区别作者逐轴Fraction变换。核对全部系数hash与最低值、7/156250目标。
3. 在临时目录编译锁定作者C++，核对生成原始图矩输入hash及space_moments解码。此步骤是复现，不称独立矩生产。
4. 对9个固定Euler图以Python任意精度Fraction直接逆多元Stieltjes分母级数，独立核对C++矩：单边2、4、12次；三角形各边1次及5次；四圈各边3次；五圈各边3次；K5各边1次；相邻两边各8次。最高16阶，不用C++canonical或整数矩递推。二进制、原始缓存只留仓库外。

解析审阅另外覆盖Young优化方向、W>=-2及Z正性、分数minorant、S5图同构、n<=17/4位编码/小int64范围、partition单项代表、G齐次扩cube及所有边界/退化输入。明确：全空间矩仍依赖作者锁定C++生成器，9个独立矩核对并非3389个矩的独立重算；不能隐藏这一范围。

以--project-root指定作者项目目录，默认该入口只校验manifest；--run才执行一次固定批次。只写本目录，不提交。
