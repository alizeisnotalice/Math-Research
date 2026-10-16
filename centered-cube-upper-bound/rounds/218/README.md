# 第218轮：全软度自适应硬面根阶几何

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮证明全软度、任意来源位移、固定任意0<L≤b的加权硬面根阶估计，并将软度选择放到每个面、每个接收点的积分内部：Σfaces∫sup_s γ_s≤C√n。核心机制是Hellinger重叠的exp[−c(n−1)(1−s)²]衰减与n(1−s)面系数共同极值化；自适应参数费用用局部包络与dyadic求和严格支付。

固定s的真实内侧上迹权∫√(J/p+)dμhard≤4√n；自适应常数明确且绝对有限，粗值较大，不声称最优。尺度和任意来源权未获免费上确界。详细LaTeX及独审均保留，不编译PDF。

辅助结果：在原mixed补域，捕获加权浮动带0<R(w_t−h0)≤εL0J/p的完整费用为16ΛL0ε/(3K)WA。它同价包含R217旧支，只替代不重复付款。超过该带的原联合交通仍未付。

桥接核验准确指出辅助t不是物理L：δ_t wvwu=wu²[−∂logt(wv/wu)]_+，来源壳层换元保全部Π、历史和有界Jacobian，但没有生成接收硬面。R170旧核外壳薄极限与hard生产的系数及内外迹亦不同。一般中间层仍需真实来源联合正变差的根阶预算；不能把本轮核工具当作已接入原余项。

## 精简交付

- [research.zip](research.zip)：17成员、101,819字节，SHA256 `84de8ceb737870a122bc02b471bb2414eaa5621de3ab8e3f7d226ce9774c4be4`。
- 包含837行详细LaTeX、主证明与必要接口附件、独审、数值笔记、三组最小复现脚本和JSON、driver及manifest。
- 仓库新增仅本README和ZIP，另更新项目README。17项工作材料在仓库外展开备份和独立ZIP均逐字节核对后移除，没有修改其他项目或历史。

## 实际验证

主代理从最终外部checkpoint运行reproduce.py，driver将脚本再次迁移到新的临时目录运行。16项manifest哈希前后均通过，三份**完整JSON逐字节一致**：

- floating_band_checks.py：18个有限clock/capture配置，Fraction精确检查与40/70位积分，0.195秒。
- affinity_checks.py：10个唯一物理配置、24精度运行，显式无限尾，12.019秒。
- adaptive_face_checks.py：三阶段8个维度至10^6的dyadic和、显式几何尾，以及n=2三个真正积分内sup物理配置，1.915秒。

解析证明覆盖全n；有限计算不拟合阶数，不是有向区间或原FULL空间证书。父级另实核39个LaTeX标签与引用；独审纠正并确认共同容量为3Q−QJ≤3Q，Π分母与绝对Jacobian记号均准确。

外部最终回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round218_20261016_driver_replay/receipt.json`。
冻结展开备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round218_20261016_final_checkpoint`。
独立ZIP：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round218_20261016_research.zip`。

复现需Python和mpmath1.3.0。解压后运行`python reproduce.py --verify-only`；完整复现用`python reproduce.py --output-dir <新的包外目录> --numeric-python <有mpmath的Python>`。本机用已有bundled Python3.12.14和vendor mpmath，未安装依赖。

规模：rounds目录217项，当前git松散对象113.91MiB、pack261.91MiB；本轮ZIP约0.1MB。已查阅[GitHub官方限制](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits)，新增规模远低于单对象及单次推送限制。提交与远端接收以Git历史和推送确认，本文不预报成功。
