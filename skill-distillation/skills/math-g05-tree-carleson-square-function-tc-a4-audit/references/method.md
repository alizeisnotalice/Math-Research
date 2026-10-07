# TC-A4 Carleson 型平方函数界：方法与验收边界

1. TC-A4的系数归一化、树序、测度与平方函数缺正式定义时，只返回缺口及标准模型，不声称私有结论。

2. 标准二进模型α_I≥0且每个J均Σ_{I⊆J}α_I≤C|J|；所有子树都需验证，根packing不足。

3. P-b7592e844240da1c的作者Theorem1.1：p>1有Σα_I|<f>_I|^p≤(p′)^p C||f||_p^p；p=2给4C，p=1强型不包含。

4. 作者sharp常数与Bellman remodeling源自其具体二进/infinitely-refining过滤模型；本次未独立证明全凹性或锐性，不将任意滤过都称同Bellman。

5. 已检查弱构件：一般测度有限过滤树有互不交孩子、统一mass流及所有子树packing≤Cσ(J)时，极大节点层蛋糕+Doob L²给嵌入≤4C||f||²。此上界不需Lebesgue或小边界；sharpness和私有树适用性另证。

6. 平方函数若定义使积分正好为上述非负嵌入和才用该界；多个算子的交叉项另证正交/几乎正交。

7. 有限深度、零测度节点和无限极限分开；中心立方体选族需先证明真实公共输入、树化和packing，不能直接替代任意中心上确界。

8. 新增复用P-6a781aff4000e334 stochastic常数e限连续平方可积、orthogonal/equal-bracket驱动及generalized CR系统，alpha条件tail mass≤1且其closed martingale在驱动stable subspace可表示；不能当一般Carleson树常数e，更未提供TC-A4。

## 不可省略的限制

目录未定义 TC-A4。标准二进 Carleson 嵌入不能冒充对私有系数规范或平方函数的精确结论。

本项目前仅提供操作方法与待核验来源；高级定理及常数须从真实阅读卡回查原文。
