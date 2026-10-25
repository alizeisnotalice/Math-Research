# 007 clean replay 最小白名单

要从零重放一个 007 单点证书，归档时只需同目录放置：

```text
portable_replay_007.py
005b_interval_heat.cpp
007_interval_fejer.cpp
007_generate_fejer_inputs.py
007_final_verify.py
body_fan_geometry.py
body_fan_arb_anchor.py
body_fan_h_epsilon_arb.py
body_fan_green_arb.py
```

以及保存本次实际执行证据的 `portable_replay_007_receipt.json`。节点位模式、
热核输出和最终 JSON 均由临时目录重新生成，不应从旧 007 输出复制。

明确排除 `body_fan_check.py`、旧浮点 fan 输出和所有 005b 吞吐基准；它们不是
该可移植重放所需依赖。运行环境须有 Python、python-flint、SymPy 与支持
`-std=c++20 -fno-fast-math -ffp-contract=off` 的 C++ 编译器。
