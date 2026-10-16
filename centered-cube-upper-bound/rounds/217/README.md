# 第217轮：捕获加权时钟带与中间层几何

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮证明全原尺度、任意合法来源下的明确正子支根阶预算：令 d_t=(3L0/K)w_g/w_t，H=(K1+K2)/J−p K0/(κT L0 J)，保原时间独立捕获余权R；在剩余mixed域内，d_t>=1且0<R H_+<=ε的完整费用为16ΛL0ε/(3K) W_A。取L0=K√n、ε=1即为根阶。

证明使用完整mixed与good的共同时间值域、完整Pi两腿、原high/FIRST列及固定J的真实dx归一；不增加尺度、soft、time或j菜单费用，不消费负信用，不假设ν<=β。详细证明及独审在研究包内，未编译PDF。

未付补域仍为d_t<1及R H_+>1。原rho0与核差在尺度可比较分支给严格二次中间区间；这只是剩余几何的必要条件，未给其一般总价格，不宣称完整闭合或新增FULL实例认证。小逆指数比较作为辅助弱工具保留，不用统一正方差替代原几何。

## 精简交付

- [research.zip](research.zip)：16成员，66,184字节。SHA256 `697759f6e62b53aef1ebca476bd52eb5df41cb015f7ac7e77ea0fde82881066f`。
- 包含详细LaTeX文本、三份关键证明附件、独审、数值笔记、三组最小复现代码/JSON、任务合同和manifest；没有PDF、缓存或重复原始文献。
- 本轮仓库新增仅本README与一个ZIP，同时更新项目README；工作材料在仓库外逐字节备份后移除，不改历史或其他项目。

## 实际验证

主代理在外部最终checkpoint上执行reproduce.py，脚本再迁移到新的临时目录运行。15项清单前后均通过；small_soft_checks、clock_certificate、capture_polynomial_check三份**完整JSON逐字节一致**。

- 物理一维tilt三阶段14配置运行：保完整M与R，非FULL、非区间积分证书。
- 六个有限min-product来源clock模型，n=2/4/8/16/32/64，保删除和完整分母；Fraction原函数与40/70dps积分核验。
- 二次中间层三阶段10,368项有理代数检查；不替代一般证明。
- 辅助MATLAB 49点核检查代码与实际显示值在SC9，仅有限浮点证据。

父级实际用时：small_soft_checks.py: 9.670s, clock_certificate.py: 0.093s, capture_polynomial_check.py: 0.083s。
外部最终回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round217_20261016_driver_replay/receipt.json`。
冻结展开备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round217_20261016_final_checkpoint`；独立ZIP备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round217_20261016_research.zip`。

复现需Python和mpmath（本机mpmath1.3.0）。解压后可运行`python reproduce.py --verify-only`；完整复现运行`python reproduce.py --output-dir <新的包外目录> --numeric-python <有mpmath的Python>`。本机使用bundled Python3.12.14，PYTHONPATH指向已有vendor，未新安装依赖。

独审解析结论、算法与证据范围全部通过；Git提交与远端接受由提交历史及推送验证确认，不以本README预报推送成功。一般目标保持未完成。
