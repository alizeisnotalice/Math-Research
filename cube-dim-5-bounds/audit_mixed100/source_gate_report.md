**100 混合截时源门与独立数值路径登记**

这是新值冻结前的源门快照；最终独审结果见本目录 report.md。此快照只完成冻结098初态的迁移、保留、typed cache、几何和缩放资格，当时没有运行新 H128 或宣布全域 g105。新生产 raw 未冻结前不读其活文件。旧090、081及既有锁合计133个文件保持原字节。

初态 SHA `19294a9bd493524fdc23ec2abe2bb3e70577f04b762b08d4fc57e78af5377480`，098准备manifest SHA `fe03bfa691b6947a24f13a2f94885962332a77edfce087e2f9e01bc49a3ace15`。099理论manifest SHA `5df7eaf9511a3b91a5a878edd7d6fb5f45ef1c69bdda4ae6886f21145ba49b3c` 已核对，其报告已完整阅读。它给出逐闭叶 G≥Hε 的共同对象资格，不替代新数值和覆盖检查。

source_gate_v2.py 使用只读081的JSON/整数几何工具，未导入作者函数。逐对象验证100804个原通过证书原样保留且leaf_epsilon=1/64；48638个历史split及52360根原样保留。190个旧queue与4个counter inflight恰迁成194个新pending leaf，全部闭盒边界由root与交织path位独立构造。四个H64障碍的原record与上端严格小于目标的事实完整保留，不重贴H128标签。根前沿prefix-free，每根dyadic整数区间精确铺满单位区间，活split恰为所有严格祖先。

H64 namespace中184796条记录逐对象等于旧180040、089正complete4752及4个旧counter record之并。H128 namespace严格为空，按(epsilon,point) pair绑定。全部194个计划请求及执行任务的epsilon、point、candidate、280投影守卫、原H64 guard、闭盒、point路径一一重建；4个同point的旧H64观测明确不作为H128缓存命中。fingerprint098独立再生一致。准备状态仍partial，无inflight/resume_jobs、无新H128记录。

第一版源门在源/缓存/mapper/树完成后，对非命中任务的显式old_H64_interval=null字段检查过严而停止；没有H/log计算或数值成功前缀。旧程序和登记保留，新v2明确检查该null。最终PID89827、10.100845083秒、exit0，source_gate_result.json SHA `7cebdfeab699d05e9f3cfb5ea09149bba4fb429cffc54a07c5846cff7374b13e`。这是源门通过，仍未通过任何新H128叶。

独立数值路线使用锁定081 arith.py原EPS64实现，但输入为2x。在同alpha1、五维、16个未单位化方向及Fourier(2π)^−5模型下，半群密度有

\[
p_t(2x)=2^{-5}p_{t/2}(x).
\]

将u=t/2代入正时间积分，dt=2du，下端1/64变为1/128，故

\[
H_{1/64}(2x)=2^{-4}H_{1/128}(x),\qquad
H_{1/128}(x)=16H_{1/64}(2x).
\]

对真实物理尾坐标求导，链式法则再乘2，得

\[
\partial_{x_j}H_{1/128}(x)=32\,\partial_{y_j}H_{1/64}(y)\big|_{y=2x},\quad j=2,3,4,5.
\]

这同时改变截时和输入，未将Hε当作同epsilon齐次函数。只读arith允许首坐标2，直接输出真实物理导数，不能省略梯度的第二个因子2。EPS128的Taylor罚项准确为1/(16ε²)=1024。

scale_gate.py实际以Fraction核对全部194个有理输入、2x输入、全部280投影和时间参数的共同缩放比。令D、E分别为x、2x的公共分母，投影整数和时间参数同时乘r=2E/D；准确检查144E/64=r·144D/128。全部新输入首坐标2，互不重复；原面输入首坐标1。该过程没有调用H、梯度或log。PID89747、0.066420583秒、exit0，scale_gate_result.json SHA `db2904086846d8eeeec20c5bfa7747d4b1371030fe585c5673e5f2f3886e0363`。

待raw冻结后另登记有限point manifest、代码及缩放结果锁，按固定有限块运行12个独审worker，先核首点质量，逐点持久保存成功前缀。每个H128独立值为EPS64(2x)的有理区间乘16，四个导数乘32；保存向外dyadic端点并逐分量严格内含原生产球。只消费194个新EPS128 pair，不重算旧180040或旧4756个H64面点。最终还须逐叶检查新Delta、physical-gradient/H、1024罚项全部16个闭角点、typed cache/journal及完整树；只有queue/inflight全部清空且每门通过，才可接083/029的同完整G资格。
