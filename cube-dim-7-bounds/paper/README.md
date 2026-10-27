# 中心立方体极大算子的覆盖估计与一般高维弱型上界

R060 中文修订：`main.tex` 与实际 XeLaTeX 编译的 24 页 `main.pdf`。

本轮包含两个原算子结论，首项系数均为 1：

- 自然 `log n ≥ 2^20` 时，`C_n ≤ Φ56 < n log n − n log log n − n/40`，准确 γ 与完整 Φ 见定理 2。
- 存在尚未显式量化的绝对 N，使所有 `n ≥ N` 满足 `C_n ≤ Φw < n log n − n log log n − n/10`，准确扣减为 `1637/12500`，见定理 3。不把前一条的显式门槛用于此条。

正文按引言、预备知识、主要证明、应用组织。主要证明重新写出扩大局部域的一维准确积分与 AC 桥、来源面数校准、目标面计数和交长的耦合、逐完整来源面的逆概率与截断、角点和资格阈值支付、同一实际历史的有限膨胀参数连接、一般 L¹/全部正共同半径及完整标量代价。宽带证明额外说明有效来源的 R0、第五个 Kw 坏事件及 Csrc 来源质量。显式版本列出全部有限 guard 和 Taylor8 尾；宽带版本列出有限个均匀 eventual guard，严格保留未知共同起点。

应用章证明新旧完整根式函数的准确排序，保留较小门槛的 R034 常数族和 R037 逐点框架范围障碍。旧 R051 长参数族正文由 R054 已发布源码身份恢复，无需重复保留全文或旧 PDF。没有将内部材料表述为已发表成果或提出优先权/全局最优性主张。参考文献沿用已核原始版本。

从仓库根目录编译两遍：

    /Library/TeX/texbin/xelatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp/cube-cover-manuscript-build cube-dim-7-bounds/paper/main.tex

新固定参数的准确转写检查（Python 标准库）：

    python3 cube-dim-7-bounds/upper/round060/manuscript_notes/verify_transcription.py

本机 60 项准确/源码检查 PASS。实际双 pass 编译成功，24 页均已渲染并看图，无公式越界、字体缺失或未定义引用。只有参考文献的两项 underfull hbox。完整正文及唯一换行增量已获独立全文数学审核 FINAL ACCEPT，锁定 TeX SHA256 `a5ffdfac3871d4ecd3d7a8e27ecb484019d62f309ce2d7a6e434754e8e17e445`。审核来源见 `BUILD.json`；作者排版与转写检查不替代独立审核。

本轮没有旧 baseline/n7、生产 grid、参数搜索、共享 MATLAB 或 Git 同步操作。历史七维脚本和证书保持原样且未运行，不作为本轮定理依据。R054 旧源 SHA256 `17e3de43f0d841a285dbfcaf8ff12c57f9670270aab9ce95cf7fe5ce359b9554`，已发布主提交 `4474b8be0cc47c904b23a5327108e37b5e6aea79`。private 恢复副本与阶段稿不在归档白名单中；不重复旧 PDF。

正式主目录的实际双pass编译与逐页渲染记录见BUILD.json；作者run另锁upper/round060/author_build.json，数学源a5ff相同。本文固定R056/R059证明范围，后续已接受成果与有效域以../STATE.md和准确账本为准。旧20页完整源/PDF用source_copy_manifest中的稳定blob恢复。
