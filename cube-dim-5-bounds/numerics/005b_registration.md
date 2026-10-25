# 第 005b 批：审计后 production-candidate 区间节点接口

本批不修改锁定的 `005_*`，仅根据审计 008 的五项必要修复建立新实现。

验收对象：

1. JSON 的每个端点以 IEEE-754 binary64 的 64-bit 十六进制位模式导出；十进制不作为端点证据。
2. 编译时 `static_assert` IEC 559、binary radix、53 位尾数和 8-byte double；运行时拒绝非 `FE_TONEAREST`，并读取 ARM64 FPCR 或 x86 MXCSR，拒绝 FTZ/DAZ。
3. `nextafter` 本身若给出非有限端点立即失败。
4. 输出 source SHA-256、完整 build command、compiler、architecture、浮点控制寄存器。
5. 对请求的全部 `count` 点导出输入区间及全部 5 个时间值，而不再只导出前 8 点。

仍声明的数学前提是：在记录的 IEEE binary64 环境与禁用 fast-math/FMA contraction 下，每次标量基本运算正确舍入；每个基本结果均作一 ulp outward。该程序是 production-candidate 节点工具，不是形式化验证的区间库，也不认证积分、Green 函数或全域界。

有限批：先 8 锁点 × 5 时间与 256-bit Arb 比较，再最多 100 点性能/导出测试；软限 120 秒。

## 实际运行（2026-10-25）

源文件 SHA-256 为
`bbc523c6c953e7baea02fd79aa97d5403950720f97342d0ebd692271aaf0c88c`。
实际编译命令为（`CUBE5_BUILD_COMMAND` 以同一命令的自指占位形式写入 JSON）：

```text
TASK_SHA=$(shasum -a 256 005b_interval_heat.cpp | awk '{print $1}')
c++ -O2 -std=c++20 -fno-fast-math -ffp-contract=off \
  -DCUBE5_BUILD_COMMAND='"c++ -O2 -std=c++20 -fno-fast-math -ffp-contract=off -DCUBE5_SOURCE_SHA256=<sha> -DCUBE5_BUILD_COMMAND=<self-recorded> 005b_interval_heat.cpp -o /tmp/cube5-005b_interval_heat-macos-arm64"' \
  -DCUBE5_SOURCE_SHA256='"'"$TASK_SHA"'"' \
  005b_interval_heat.cpp -o /tmp/cube5-005b_interval_heat-macos-arm64
```

其中 build 还以预处理器字符串参数写入上述源 hash 与完整命令；结果 JSON
记录 compiler `Apple LLVM 21.0.0 (clang-2100.3.34.2)`、`aarch64` 和
`aarch64_fpcr_0x0`。运行时 `FE_TONEAREST` 与 ARM FZ/FZ16 检查均通过。

最终 8 个锁点 × 5 时间实际运行耗时 `0.883605 s`。使用输出的**原始 binary64
位模式**（不经过十进制端点解释）对照独立 256-bit Arb，40/40 均包含；最大
区间宽度 `6.751845420651173e-12`，最大相对宽度
`5.2568483198844945e-08`，复核耗时 `25.312013 s`。

100 点 × 5 时间的有限导出测试耗时 `10.829 s`（约 `46.17` 节点/s）；
`005b_interval_benchmark_100.json` 实际写有全部 100 个输入区间与 500 个结果
区间。随后重新运行 8 点并再次完成上述 Arb 交叉检查。所有运行均在 120 秒
软限以内。可执行文件位于 `/tmp/cube5-005b_interval_heat-macos-arm64`，未放入
证据目录。
