# execution135：主已完成唯一标量的源绑定终核

**源绑定通过；没有重新运行135标量、四次根、七条备用比较或成本残差表达式。** 135 冻结初稿、源码、计划与 manifest 保持原字节锁；本执行附录独立保存。

主在自己的 `round015/join135_actual/` 实际运行了一次已审 readonly consumer。原记录 PID **51083**、耗时 **0.00035354201099835336 秒**、case_count=1、pool=null。exit0 由主的可信任务消息给出；本分支不从进程消失推断退出码，也不补跑原程序。

本分支只读终核实际 PID **51975**、exit **0**、约 **0.00295 秒**。不导入135 consumer，不使用 Fraction/SymPy/Arb 或任何数学求值器，只做 SHA 与 JSON 字段一致性检查。

## 准确输入与来源字段

|输入槽位|实际锁定来源字段|终核结果|
|---|---|---|
|Cbar_upper|134 `condition_results.json:Cbar_upper_256dyadic`|输入、actual结果和134 manifest逐字符串一致|
|M_upper|132 `portable_components.json:objects.analytic_stationary.M_upper_portable_cap`|输入、actual结果和132 manifest逐字符串一致|
|IH_upper|014 `primary_integrals_acceptance.json:complete_H512_integral_upper`|输入、actual结果和132 portable `UH`逐字符串一致|
|π区间|076 `results.json:derived_constants.pi_outward_384_dyadic` 的 lower/upper|actual结果两端点逐字符串一致|

M 的实际选择是132已保存的准确有理 portable cap

`488568596392024500738533860719/250000000000000000000000000000`。

它没有从 display 小数构造，也没有被误称为更紧的原 raw `M.upper`。本终核不重新优化 M，不重新验证其数值算术；数学与独立算术资格由原132、136和主承担。

三项资格元数据均明确 ε=1/512、domain=complete_F。Cbar 绑定整个短时区间，M 绑定同 ε 的完整 H512 sup，IH 绑定同 ε 的完整 H512^(5/4) 面积分。原 g105 与029传递文件也与冻结计划/运行注册中的 source SHA一致。

## 单次执行链

已逐项核对：

- 主授权文件的 code/plan/input SHA、唯一 case_count=1、root_bits=192，与 actual注册及冻结文件一致。
- actual注册中的全部 plan sources 与冻结计划相同；每个 source 文件当前 SHA 匹配，actual注册的 qualifications 与准确输入完全一致。
- actual结果的 Cbar/M/IH/ε/PID、registration SHA 与原输入和原注册一致。
- actual结果保留 `global_I_g_Cn_claim=false`、没有实际 J/R/I 的下界主张、`all_sources_unchanged_after=true`。
- 原结果记录 `conditional_cost_upper_below_32=true`，保存的准确有理残差分子为正；本分支只读取该符号，不重新计算残差表达式。
- 本次只读终核前后全部锁定源 SHA 未变，冻结135 manifest仍为 `bc05ffe0e93315c764928f7e9b34dde450cccdc47777a5b223f8cffa41772d54`。

原根上下端点与完整条件 I 上端按原字符串复制进本目录 `results.json`，供来源追溯；它们没有在此重新计算。准确根/残差的独立数值核验属于136及主的既定任务，本文不替代该门。

## 原件锁与归档范围

|原件|SHA-256|
|---|---|
|135准确输入|`c6eb027efee798f140f6cb7362721b246de257a722b58cf65ba57ea1a6ad426c`|
|135 actual registration|`b1ca8c818d5cc0928141f83a72fc21f4f199792e0bc28f8b5d35d1fccf4b1bab`|
|135 actual results|`8b0623db708b493cf7d196e5f8dc0f88724d401154a7c888da95987568502686`|
|冻结 consumer|`80c162f2d1f90d6cc88140c700ce94fe83c776a0157ef279632b0fce17947512`|
|冻结 plan|`98a200a9645f2532252af43ed03ac4871346aa95b2cc028aa3a835578c0df50d`|
|134 condition|`1cf1aa32b44eb580226692bfe37ab40d2ef96deeda4edc155f16a2f0a48cb024`|
|132 portable components|`2ec64713d80f2ecc791a7207a706d035b4c6723c24868cc44383c4d91ff9855b`|
|014完整IH接纳|`14374097c6d164fcbaca56a4b6d003069e92e5e16ffa990a2ebb6cb470d2c219`|

只归档本报告、纯来源绑定脚本、终核注册、结果和 manifest。原实际预算、主七条备用证明和 MATLAB 运行由主 `round015/` 按原来源保存；本分支没有重跑或复制这些大/重复材料。

136与主029/完整分析接纳仍待完成，本分支不宣布完整 I、g 或 Cn 更新。0新H/R/梯度/log/exp、0新M优化、0真实标量重算、0七链重算；不启动下界或Astra，不提交、不推送，不修改他人文件。
