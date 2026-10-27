# 中心立方体极大算子的覆盖估计与一般高维弱型上界

当前中文正文为 `main.tex`，24页 `main.pdf`。结构为引言、预备知识、主要证明、应用。核心结论为log n≥2^108（即s≥2^36）的整数 n 下

C_n < n log n − n log log n − 2^(−32)n/(log n)^(1/3)。

正文还保留 log log n≥128 时 C_n<T_n−n sqrt(log log n/log n)/131072 的完整显式上界，分别按适用范围取最小值。

正文保留有限覆盖及原角点公式的完整证明，新增二阶矩列截断与新参数的完整证明、反对角多壳的共同增量能量证明、有限 a 捕获余量、全部14个门槛包络及无限尾、有理基点、较小门槛参数族与真实选择历史的逐点框架限制。后者不是算子常数下界。未开展七维搜索或下界计算。

在仓库根目录编译：

    mkdir -p /tmp/cube-cover-manuscript-build
    /Library/TeX/texbin/xelatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp/cube-cover-manuscript-build cube-dim-7-bounds/paper/main.tex

重复运行至交叉引用稳定。基点代数的独立准确核查：

    python3 cube-dim-7-bounds/upper/round041/manuscript_notes/verify_manuscript_base.py

它只执行 Fraction 有理基点运算，不调用七维旧认证器。计算范围及结果见该目录 `exact_base_receipt.json`；来源与逐项转写核对见 `transcription_checklist.md`、`sources.json`。编译与排版结果见 `BUILD.json`。

原13页稿可从提交 `a3cd1967f4879b7398535691dcfac6dc7d38a724` 恢复，身份记录见 `upper/round041/manuscript_notes/chapter_plan.md`。`constant_certificate.json`、`VERIFY_RECEIPT.json`、`verify_constant.py`、`dimensions_table.tex` 和原 `manuscript_review.md` 是该历史稿的材料，本次未运行或重认证，不是当前一般高维定理的证明附件。公开引文只使用已经归档的 Stein–Strömberg 1983、Aldaz 2017 v2 与 Aldaz 2020 v2 原始版本。

R042二阶矩符号转写核查需要SymPy；本机实际使用1.14.0：

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round041/manuscript_notes/verify_second_moment_transcription.py

其范围仅为单坐标积分、固定面因子和准确常数转写，结果见 `second_moment_transcription_receipt.json`。新主定理已接入最终独审接受的显式s≥2^36门槛，正文逐项给出几何、行乘积、二阶矩补偿、Taylor8及无限尾证明。接入前17页阶段TeX仅保留在作者独立clone中，不重复归档阶段稿。

新显式门槛准确核查（仅标准库，无外部路径依赖）：

    python3 cube-dim-7-bounds/upper/round041/manuscript_notes/verify_second_threshold.py

实际22项全部通过，结果见 `second_threshold_receipt.json`。应用中同时保留新旧参数下完整上界函数及准确差式，不以显示小数或扣减项大小判定全函数排序。

e523全文数学独审后的三处定义已补齐：角点积分首次定义、完整来源带测度与face径向因子、能量引理M>0。独审原源仅保留在作者独立clone；也可在临时副本对最终main.tex反向应用notes/definition_delta.patch恢复，原源hash见sources.json。仅重新编译与查排版，未重跑数学或七维计算。

最终独立增量审阅已锁定94c5正文并给出FINAL ACCEPT，三处定义事项全部关闭；报告hash及路径见notes/sources.json。数学源独审与作者24页实际编译/排版检查分开记录。

主代理整合编译与作者运行证据分别记录在BUILD.json及upper/round041/author_build.json。推荐复核命令使用现有数学venv；作者当时实际SymPy路径保留在原始记录，不把推荐命令写成重跑证据。
