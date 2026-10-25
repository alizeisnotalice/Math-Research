# 009 — 005b 区间节点修复复核

本轮只审计 005b 对 008 所列缺口的修复，不重新审计既有 body-fan/热核
公式。审阅源码 SHA-256：
`bbc523c6c953e7baea02fd79aa97d5403950720f97342d0ebd692271aaf0c88c`。

## 结论

通过，状态可升为**可重放、平台条件性的 production-candidate 节点工具**：
其端点格式、已记录的 ARM64 环境和指定编译假设下，8 个锁点的 40 个
节点已有直接 bit-pattern Arb 交叉，100 点/500 节点也已完整导出并独立
重放。它仍不是形式化验证的通用区间库，且不能单独认证积分或 \(C_5\)。

## 实际独立复核

审计者用 Apple clang 21.0.0 在临时目录以
`-O2 -std=c++20 -fno-fast-math -ffp-contract=off -Wall -Wextra -Werror`
重新编译，带源码 SHA 宏和审计用 build-string 宏。结果：

* 8 点运行报告 `aarch64_fpcr_0x0`，五个时间 bit pattern、40 个输入区间
  和 40 个输出区间逐位匹配 005b 生产 artifact。
* 对临时重编输出实际运行独立 256-bit Arb 脚本：40/40 包含，最大相对宽
  `5.2568483198844945e-08`（约 25.6 秒）。比较直接读
  `lo_bits`/`hi_bits`，没有经十进制端点回读。
* 重编程序实际运行 100 点（约 11.0 秒）；100 个输入区间和 500 个结果
  区间逐项逐 bit 匹配 `005b_interval_benchmark_100.json`。JSON 结构也
  实测为 100 行，每行 5 输入与 5 输出。
* 独立以精确二进制有理数检查 8 个锁点的 40 个输入端点，均包含原始
  \(1\) 和 \(\mathrm{TAILS}/1009\)；输出的 5 个时间 bit pattern 正是
  \(1/128,1/32,1/8,1/2,1\)。

## 对 008 修复项的判定

1. **端点协议：通过。** 输出使用 `std::bit_cast<uint64_t>` 的 16 位
   hexadecimal bit pattern，JSON 不再把十进制当数学端点。`005b_arb_crosscheck.py`
   按该 bit pattern 还原 binary64 后直接比较 Arb。
2. **IEEE/环境前提：通过，但为平台条件。** 编译期断言 IEC 559、radix 2、
   53-bit significand 和 8-byte `double`；运行期拒绝非 `FE_TONEAREST`。
   ARM64 读取 FPCR 并拒绝 FZ/FZ16，x86-64 分支拒绝 MXCSR FTZ/DAZ。实际
   artifact 的 FPCR 为 0，故本批无 subnormal flush。
3. **`nextafter` 溢出：通过。** 输入及 `nextafter` 结果均作 finite 检查，
   不再可悄然产生 Inf 端点。
4. **全点导出：通过。** `count` 的每一个点的输入 bit interval 和五个结果
   都写入 JSON；不再有 005 的“算了但未导出”问题。
5. **\(\pi\)、复运算、sum：无新的回归。** 仍保留一 ulp outward 的 hex
   binary64 \(\pi\) 包围、正实部的倒数拒绝、complex interval 乘法及每次
   累和 outward。这些路径已由本轮重建/Arb 和 008 审计交叉。

## 仍须明确的限制/建议

* 接受的 IEEE 前提应准确写成：在该记录的 ARM64/Apple clang 环境，
  `FE_TONEAREST`、FPCR FZ/FZ16 关闭，且所用 `+,-,*,/` 为 binary64
  correctly rounded；`-fno-fast-math -ffp-contract=off` 是实现这一前提的
  构建条件。静态断言和运行检测支持这个范围，但不是对任意编译器、优化器
  或硬件的形式证明。
* 结果的 `source_sha256` 和 `build_command` 是编译时宏注入，程序没有自行
  重新 hash 源文件；而 JSON 内 build command 的 `CUBE5_BUILD_COMMAND=<self-recorded>`
  是自指占位，不是逐字符的 shell transcript。独立重建与 artifact bit-match
  为本批提供强交叉证据；若要更强 provenance，应存一份不可自指的构建脚本
  和其 SHA、编译器版本/二进制 SHA，以及命令 stdout/stderr。
* Arb crosscheck 目前逐项验证 40 个锁点输出，却未在其自身代码中 assert
  JSON 输入端点包含精确 rationals、时间 bit pattern 正确、端点有限且有序；
  本审计已独立完成前两项。建议把这些检查纳入 crosscheck，使证书自足。
* 输出流未检查 `ofstream` 的打开/写入失败。应在正式批量节点生成后检查
  stream 状态，避免磁盘/路径失败时将不完整文件误标成功。

以上限制不妨碍将 005b 用作下一轮的有限节点工具，但使用其值时必须保留
“平台条件性 floating interval”标签，直到改用形式验证区间库或提供等价
的机器验证证明。
