# 原selector实验与最终区间修复的审计索引

质量、DP、prior、周期到L1截断、有限查询格认证见 `PROOF.md`。独立公共checker在 `audit.py`，实际小例结果在 `AUDIT.json`。

`INDEPENDENT_SELECTOR_CHECK.py`核对主n8实验真实盒核和DP，24个真实query×7尺度逐点用Fraction交叉检查，重放训练后selector、阈值、π、training records全部一致；还用独立kernel和字典DP重算全部16384 test query。质量与公式通过，但原严格λ1计数重算失败，最终JSON状态保留为 `STRICT_BOUNDARY_FAILURE`，不能误读成最终区间修复失败。

失败原因是source半周期使fixed/prior在L1.5真实平均恒1，原浮点严格 `>1`计数对1ulp敏感。`BOUNDARY_DIAGNOSTIC.py/JSON`保留n8、balanced n7/n15全部原测试的独立重放、阈值guard及精确L1.5修正结果。每例另有20真实query、全部table/coordinate的Fraction恒等式检查；n7/n15全部旋转type轨道和π=1/4精确验证。guard不是认证。

最终主区间修复由 `INTERVAL_REVIEW.md`审查，独立小checker与实际结果在 `INTERVAL_CHECK.py/JSON`。该修复保留原freeze和query，认证局部source/DP算术，全部18model成功数上下包相同；修复后prior差仍保留，最佳fixed table仍更强。36个CI端点高精度重算通过，CI及连续总体积分仍保留MC/数值资格。未取得C_n新增长阶。

`DP_ERROR_BUDGET.md`给出可选的正DP gamma深度预算；主修复直接逐操作向外，不依赖此替代方案。

复现本目录证据：分别运行 `python3 audit.py`、`python3 INDEPENDENT_SELECTOR_CHECK.py`、`python3 BOUNDARY_DIAGNOSTIC.py`、`python3 INTERVAL_CHECK.py`。这些程序只写本审计目录；原source、frozen selector及原RESULTS保持只读。INDEPENDENT_SELECTOR_CHECK的JSON已在更新conclusion后实际重跑，明确否定原naive计数的资格；不是仅改历史JSON文字。
