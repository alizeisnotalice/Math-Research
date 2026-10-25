# 11J：三维minimum高阶连续证书有限批

2026-10-25；承接continuous_minimum首批report SHA8278d9ad85b8ab4029b5bb9ae289d1807e771500b6972b0d8191233eb3f6c941及其manifest，不重跑4095直接区间失败批。原GPT6-astra/high、原独立clone，只upper/continuous_certificate/cc11_cheb_*，不覆盖已锁原件。

固定tau1/3、真实H_{1/50}≤G；triangle face用16×16方形中的136patches覆盖。每坐标degree20 Cheb roots，rho3，已得全纯tube统一误差<2.868314e−7。约59976node。先测4个代表patch含最危险edge、diagonal/confluence和远离奇性的patch，确认无理nodes包围/point evaluator/多项式下界的完整路径；若可验再完成固定136patch，不新增tau或更密均匀网格。

节点误差必须可认证：可用Arb无理node enclosure，或准确有理近似节点加全纯Cauchy导数位移误差。全域验收需严格多项式下界（区间Bernstein或带余项分支界）；仅节点通过无效。若Bernstein转换放大误差，优先提高同节点算术精度，不重新搜参。最多一次明确原因驱动的表示/精度修订；每patch子域分支须先登记有限停止条件，不能无限细分。

先形成实际可运行端到端4patch基准及hash/runtime，随后固定全patch批。跨patch共享缓存，原节点值只保存必要证书或仓库外可再生缓存。成功须每patch严格下界>13/2500并覆盖全triangle，然后主/独审验证。失败保存未覆盖patch与误差项，既不是反例也不结束目标。MATLAB交主串行，Python/Arb主动真实运行；不动g以外J/其他项目/总账，不提交推送/扩招。主负责最终整合与一般L1/allr。
