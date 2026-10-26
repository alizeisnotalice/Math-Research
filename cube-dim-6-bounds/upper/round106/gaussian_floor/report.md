# 第106轮：Gaussian 全域下包的结构障碍与反射平均替代量

仅原 body32、等权 Cauchy 模型；下界研究停止状态保持。正式 C6 上界仍为64。本轮无全域 g0=1/65000 证书、无完整新I、无正S。

## 结论

1. 正 Gaussian scale mixture 的准确系数为
   G(x)=Gamma(5/2)/(4pi^3) E[det(Sigma)^(−1/2)(x^T Sigma^(−1)x)^(−5/2)]，
   Sigma=sum S_v vv^T，E exp(−lambda S_v)=exp(−sqrt(lambda))。原建议的分母2pi^3需改为4pi^3。条件 Gaussian 的协方差是2t^2Sigma，保持同一真实 Cauchy 时间。
2. 先逐样本取整个方体最小、再取期望是合法正下包，但它本身有严格上限：
   B_Gauss<=3/(256pi^3 6^(5/2))<1/165888<1/65000。
   Arb辅助值约4.286013605856154e−6。因此任意精确期望积分、重尾shape消除、正矩阵事件或此量的进一步下估，都不能让这个特定机制证明g0。该结论不排原G、其它g机制或其它S计划。
3. S_v=1/(4Y_v)可将一个径向变量准确积分掉，留下31维Dirichlet(1/2,...,1/2)形状期望；其边界集中仍须处理。该自归一化本身不绕过第2点。
4. 一个具体不同控制量是先平均64个坐标反射后的逆二次型势、再对共同x取方体最小。其期望仍合法下包真实g。对实际body32正矩阵Sigma=I+100vv^T，v全1，可严格降为5条一维分支；有限核验全部通过，得到H_orbit>1/(20sqrt(601))，严格高于该矩阵的原trace上限6^(−5/2)/sqrt(101)，也超过原单样本方体下包四倍。该例不替代随机期望认证。

完整推导、all-cube合法性、真实S2奇异层兼容、同生成元/source口径，以及反射后的min不得再偷换成corner-max，见theorem.md。主代理已在运行前通读并放行这些归约；最终证明责任及后续独审由主代理保留。

## 实际验证

唯一放行程序check_gaussian_floor.py运行1次，退出0，内部wall=0.374116291秒，tool chunk dcdbb1，0子进程，无常驻session。它做Sympy准确常数化简、整数tight-frame/矩阵逆核验，以及5×64=320个闭区间的192-bit Arb下界。每项均严格>1/20，最小展示值约0.09548241606。展示用浮点选择最小行，不参与验收。

源码SHA256：27ccca5717a0feaf12e3f4dd4dc60b5f678fee317c467a0cee0d6a1582db93e5。
输入SHA256：18c9ea9b6dc36491a3ae25625213f4b18fbda847d1da08b96408b4a01a44bd37。
结果SHA256：62bd1e050bf91e8aacfb5f92a94ecb73c3d621d4526b23cf656c9bb8ffb765a6。

results.json保留320个有限行，约80KB。run_receipt.json保存退出及清理信息；primary_release.json保存运行前静审绑定。原源码自设SIG_DFL177 alarm留到exit；这是本次有限短反馈选择，用户已取消一般程序时限，未将旧上限冒充现行要求。

check_saved.py及mincheck_inputs.json仅准备、尚未运行，待主代理另行静审。它只核保存证据哈希、闭网格覆盖和已保存区间端点，不重跑producer，不充当独立数学证明。

## 范围与来源

限定查重见sources.json。006已有一般alpha Gaussian混合表示，本轮不把表示本身称新；新信息是先min后E的结构上限与反射平均控制量的严格有限例子。

专用clone在此前不存在的Math-Research-cube6-sol-gaussian106路径创建。clone本体建成，但checkout因源仓库promisor缺blob失败；只使用已可读的project文件及独立.git，不声称完整checkout，不重试基线。clone失败收据与初始index异常在registration.json保存。所有新文件仅在本目录；未stage、commit、push，未用MATLAB、Astra、子代理、MC、旧kernel、091粗基线或下界程序。
