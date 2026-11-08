# 第77轮来源与依赖

本轮没有新增外部端点上界定理。有限网格归约、最大源合并、对偶证书和两角点反例均在本轮给出直接证明；未登记公开优先权。

- [r076近角点估计](../r076/nearcorner76.md)、[必要源复用](../r076/entropy76.md)与[分离副本](../r076/copies76.md)：本轮处理它们留下的有限分组与源预算问题。维数一致的全局预算仍未证明。
- [r023缺角实际源](../r023/report.md)：本轮notch输入取其cap2范围的同类两层形状，只用于分配诊断，不把旧下界宣称为本轮改进。
- [r073一般振幅归约](../r073/amplitude73.md)：一般输入归约到cap2的模块已闭合。
- [r075固定角点](../r075/corner75.md)：每个真正成功的固定角点源片可作全局Hardy结算；分数权重不能被偷换成降低该源片的成功阈值。
- [Wolfram LinearProgramming官方文档](https://reference.wolfram.com/language/ref/LinearProgramming.html)：核对m.x≥b及非负变量约定、精确有理输入和精确算法。实际脚本还检查返回解全为Integer/Rational并逐条精确检验可行性；主要结论依靠独立的有理原始/对偶证书。

allocation77.md建立模型，review77.md含同角点嵌套KKT基准、分离可加条件及残差对偶审查；integrality77.md证明裸分数分配障碍与perspective修复。软件求解状态不替代对偶证书，有限分组最优不替代整个halo或最优分组。
