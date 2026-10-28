# R092：单侧源帽计数引理与临界缩放诊断

LEMMA_ACCEPTED / CONDITIONAL_MODEL；本轮没有新Cn界。
主全文核对两份原件：若w<a0、|w|≤r且2rho>r+1，则源lower strip
[−rho,−rho+1−q)完整落入捕获区间；不需要旧对称条件|w|<a0。
新增第六个full-target计数事件，失败≤4/(e²n)，n≥2/e时全matching
count至少nu_one n、nu_one=(1−e)/2。六事件各≤e/8仍保相同目标质量。
原J/J1/J2、bounded weights与完整column不变；只替换weightedcap中的count。
其应用到新算子构造和显式门槛仍需组合全链认证，不能凭此单独改账本。

临界缩放s~lambda L/x、x固定的诊断保留原B的精确积分极限
I(x)=exp(x/2)E1(x/2)，DCT的近零可积包围已在文中给出。
在额外典型pin-sign和transition-count校准及整合误差假设下，
模型源帽损失为exp(x/2)E1(x/(2Q))/2。当Q↑1，该损失正好抵消
B内I(x)/2收益；Q<1时logQ加残余项严格为负。
这些额外假设没有作为已完成的n依赖s全链；结论仅限条件模型，非所有覆盖不可能性。

作者预登记仅一个x=2、lambda=3、Q=.75案例，本机运行PID33817，
准确log验证及SciPy积分一致性PASS。naive系数1.109103897…扣除源帽成本
.174872920…后，条件模型系数.934230977…；无新Cn改进。主核了全strip端点、
计数中pin损1、n≥2/e代数、theta消去、column限制与条件积分的方向。

原件/最小复核/实际receipt/来源hash保留。没有重扫网格、旧基线复跑、论文修改，
下界及七维专项不启动。公开与内部上界均不变。下一步可单独认证第六事件与
Bernstein/新带宽组合，以判断是否能降低显式门槛；当前仅有局部引理及条件诊断。
