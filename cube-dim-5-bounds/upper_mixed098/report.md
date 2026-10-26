# upper_mixed098：194 个 H128 叶的有限生产与混合 G105 覆盖候选

首批 194 个新 H_{1/128} 锚点全部完成，194 个闭盒均通过本叶曲率罚项 1024 的 16 角点证书；没有新增二分，queue=0、inflight=0、fault=0、H128 目标障碍=0。保留原 100804 个 H_{1/64} 通过叶，得到 **100998 个通过叶、48638 个内部节点、52360 个 ordered roots 的完整 mixed-G105 覆盖候选**。目标为 21/200000，同一 alpha1、16 个未单位化 body 方向、Fourier 因子 (2π)^−5 模型。

这是实际 DD/Arb 生产和作者核验结果，**不是独立数值审计**。主代理及 audit_mixed100 保留最终资格；本目录不宣布新的 C5 或完整积分 I 数值。旧 g100 资格和四个固定 H64 障碍原样保留。

## 数学与来源

075 的任意正 ε 全 frame 估计为

`log Hε(c+u) ≥ log Hε(c)+∇log Hε(c)·u−|u|²/(16ε²)`。

每叶使用同 ε 的 H 和四个物理尾导数、有向区间梯度比及向下 log(H_lower/g)。安全线性项为各区间线性项的最小值，减去正二次式后全局凹，16 个实际闭角点足够。ε=1/128 的罚项为 1024；旧 ε=1/64 叶仍用其原证书和原罚项，未重算。

099 已独立审查逐叶 `G≥Hε` 的异质截时拼接和 G 的 −4 齐次延伸；不能对混合 H 主张同 ε 齐次性，也不混合 H64/H128 梯度。完整理论来源为 audit_mixed099，manifest SHA `5df7eaf9511a3b91a5a878edd7d6fb5f45ef1c69bdda4ae6886f21145ba49b3c`。主核记录位于主 clone 的 `round013/mixed098_primary_review.json`，其哈希写入本次实际 registration。

唯一输入为外置 `initial_mixed_state.json`，SHA `19294a9bd493524fdc23ec2abe2bb3e70577f04b762b08d4fc57e78af5377480`。它恰把 089 的 190 queue 与 4 个保留 counter inflight 转为 194 queue，未丢失未证叶。原准备包 manifest SHA `fe03bfa691b6947a24f13a2f94885962332a77edfce087e2f9e01bc49a3ace15` 不变。

H128 standalone kernel SHA `128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0`。主核已逆应用限定 diff 核对其与 046 原 prefix 逐字一致，变化仅独立说明、ε=1/128 和时间参数 `144D/128`。数值核仍含 280 个 generic 投影资格、3504 个完整锥、同一 DD/Arb 公式和全部 H+4grad。256/384/512 仅用于数值资格梯级；本次全部 194 点在 256 bit 完成，每点一次 segment。

## 实际有限运行

只运行一个固定 194 输入程序，不设置 wall alarm，不自动续批。`run.py` 先拒绝已存在 batch 目录或任何 point checkpoint，登记 source、全部源码/理论/表哈希、实际 master PID 和点集，再建立并登记 12-worker pool，之后才派 H128。任何未来恢复必须先核 ε、kernel、point、precision 和 cone prefix；本次没有重启。

实际 master PID **90096**，workers **90116–90127**；全部正常 close/join、exit 0，并由 OS 确认不存在。完整程序含加载与最终落盘耗时 **25.134462083 秒**。194 个新锚点、194 个完整 H128、776 个物理梯度分量；0 个 H128 cache hit。4 个同点 H64 历史记录仅作跨截时单调性检查，没有作为 H128 值或梯度复用。

缓存严格分离为 `cache_by_epsilon['1/64']` 的 184796 个历史观测（含原四个 H64 障碍）和 `cache_by_epsilon['1/128']` 的 194 个新观测。记录、checkpoint、证书都显式绑定 ε。

派发前已比较 095 冻结 jobs 注册的 H128 时间与完整五坐标有理点：098 的 194 个新请求与 095 的 13406 个新 H128 点**精确交集为零**，不重复其并行 H128 零阶请求。比较收据 `overlap095_verification.json` 保留 095 registration/jobs 哈希。

新 194 点最小 H128 下端约为 **0.00011871452571952959**；最小新角点 log margin 约 **0.11312666511154411**。小数只供阅读，证书保存完整有理端点和余量。

## 作者终核与可携来源

`verify.py` PID **91072**、exit 0、耗时 **7.784022541 秒**，零 H/gradient/kernel/Arb/log 调用。它用纯整数 768 bit 向下固定点、180 项正 atanh 级数检查全部 194 个 Delta 下端，用 Fraction 重算 3104 个完整闭角点，核同 ε 梯度比包围、锚点物理盒与记录/checkpoint/journal 一致性。另直接比较全部 184796 个旧缓存对象、100804 个旧通过叶及全部保留顶层字段，确认原样保留；核 583 个 hash-chain 事件与 194 个 dispatch/receive/merge；从 ordered roots 和 root/path 重建闭盒，核 52360 根精确 Kraft 和反链以及全部 split 父子关系。

冻结外置目录：

`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_mixed098/`

关键冻结对象：

- `batch0001/manifest.json`：`0333c44e9d8bc29f2f1c2ac16ab38bbc935af8ae738f334a767812aed4a61928`。
- `batch0001/state.json`：`a32a7e7ac99058223812554dfaa8f3c8b54fffdb78442d0a94b5d212bb83a07a`。
- `batch0001/journal.jsonl`：`cc0105f13d681fac9b7a7573704fd0dd8242b736e40ef0bebb040ee55bf1cfb9`。
- `batch0001/point_manifest.json`：`293622876e07bb998227bad5b1a2f14af56a9727d17e3242023a43ca77d7d06c`，含全部 194 raw 及每点 checkpoint 哈希。
- `final_state.json`：`73bc873f23f7b68f6f70132e590d28ebbe4cd619012ac9c81980f9fa500314af`，完整物化状态，不依赖活进程。
- `new_H128_leaf_mapper.json`：`6c7078e87762e7d8a3cf818acd65b5fce8867c688fa7af6f620ca9f26dac0653`，完整 source queue→新叶/记录/证书/物理盒映射。

最小可携消费集合为冻结准备输入及来源 mapper、batch registration/initial overlay/journal/state/point manifest/point checkpoints、最终完整 state、作者收据、075/099/090 来源锁及本目录源码；生产 manifest 给出精确路径和 SHA。未重复复制历史大 cache。源码不对 fresh-only 生产提供盲目重算命令；审计可从固定 journal、raw 和 final 进行零 H 重放。新数值独立重核由 audit_mixed100 完成后交主验收。本次没有提交、推送、使用 MATLAB 或改动他人文件。
