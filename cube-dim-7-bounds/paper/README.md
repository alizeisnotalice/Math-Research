# 中心立方体极大算子的覆盖估计与一般高维弱型上界

R068 中文修订：`main.tex` 与实际双轮 XeLaTeX 编译的 28 页 `main.pdf`。

本文证明两个显式原算子上界，所有对数均为自然对数：

- `log n ≥ 2^20`：`C_n ≤ Φ56 < n log n − n log log n − n/40`。
- `log n ≥ 2^30`：`C_n ≤ Φp < n log n − n log log n − (23/100)n`，其中 `t=1000n/(917log n)`，`Φp=F(A_n(t),B_n(t)−(317/1000)n)`。

首项系数仍为1。第二个显示常数来自同一个完整 Φp 的准确对数价格，并未另换核、gamma或完整代价函数。门槛是充分条件，不声称最小。

正文按引言、预备知识、主要证明、应用组织。R056较早显式域的完整证明保留；R059旧宽带长证明替换为新相对间隙J与连续匹配、目标一阶Jphi剖面、负ellpin导数下端点计数、同一目标剖面控制每张来源面均值、概率先逐面后混合、径向行列measure校准、两个有界权支付原资格与同一实际历史。新显式门槛自足给出所有几何/五坏事件/概率/径向条件、Taylor10 cap及tau全尾与Taylor4完整根式余项。一般L¹和所有正共同半径通过原有限覆盖和真实见证连接。

应用比较完整根式函数和准确正价格，另给不交点支撑缺额相加命题与共用来源的反例；R034低域常数族及R037真实history框架范围完整保留。公开ABC与原创编号分开，参考文献沿用已核版本。正文不使用全局最优性或公开优先权主张。

从仓库根目录编译两遍：

    /Library/TeX/texbin/xelatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp/cube-cover-manuscript-build cube-dim-7-bounds/paper/main.tex

新固定tuple及门槛的独立转写检查（标准库Fraction）：

    python3 cube-dim-7-bounds/upper/round068/manuscript_notes/verify_transcription.py core
    python3 cube-dim-7-bounds/upper/round068/manuscript_notes/verify_transcription.py threshold

本机本轮core23项、threshold41项分别实际PASS；receipt保留各实际运行阶段的源码与脚本身份。最终源码另作引用/文献/保留块扫描，不重复数学检查。没有重跑旧60/128/producer baseline、参数grid、生产计算、七维脚本或共享MATLAB。

最终双pass编译成功，28页全渲染/看图，无越界、缺字或未定义引用，仅文献两项underfull hbox。独立完整数学审核及两处定义范围补全均获 FINAL ACCEPT，锁定源码 SHA256 `8135adfb51470228ed18371c6c5f1a22caac0b35523b0cd49c50b13cfb9b386b`，报告身份见BUILD.json。作者排版与转写检查不替代全文审。

R060旧源SHA256 `a5ffdfac3871d4ecd3d7a8e27ecb484019d62f309ce2d7a6e434754e8e17e445`，已发布主提交`2f393a17f443f10411fa9ea435fd88fd5350c14a`。旧R060notes和历史七维附件未改；private恢复稿不进入归档白名单，不重复旧PDF。来源、转写责任、实际检查、归档范围见R068 manuscript_notes；主负责最终整合、独立复核与归档提交。

正式主目录再次实际双pass编译及全部28页渲染复核，见BUILD.json；作者运行单独锁定upper/round068/author_build.json。新研究第二剖面候选仍为独立轮次，不进入本冻结数学源。旧24页稿由source_copy_manifest的稳定Git blob恢复。
