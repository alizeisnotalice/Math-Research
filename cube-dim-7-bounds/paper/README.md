# 中心立方体极大算子的覆盖估计与一般高维弱型上界

[论文 PDF](main.pdf)及[可编辑 LaTeX 源码](main.tex)，共37页。全文以引言、预备知识、主要证明、应用组织。

对自然 L=log n≥8192，本文证明 Cn≤F(A_n(t),B_n(t)−291n/1000)<nlog n−nloglog n−47n/250，t=256n/(231L)。对 L≥2^30，另一组参数给出更强的完整函数界 Cn<T−23n/100。对象是 Lebesgue 测度下中心轴平行共同半径立方体，允许一般 L1 输入和全部正尺度。

正文给出两阶目标剖面、逐来源面 Cantelli 保留、同一有限选择历史和完整行列缺额、全部显式门槛与无限尾及完整根式代价。按用户指定的 anti-defensive-writing 技能修订表达；数学适用条件和真实依赖均保留。

最终源获独立全文数学审阅；作者与主分别实际完成两遍 XeLaTeX 编译和全部37页渲染图审。无 overfull、未定义引用或缺字；有3项可读的 underfull 警告（正文1项、参考文献2项）。两次独立编译的 PDF 元数据可能不同，源码相同，见 BUILD.json。

准确证书、原始材料和最小复核代码见 [补充索引](../upper/round077/manuscript_notes/supplement/README.md)。脚本会写自己的回执；为保留原回执，可在仓库根目录创建临时副本执行：

```sh
paper_check_dir=$(mktemp -d)
mkdir -p "$paper_check_dir/paper" "$paper_check_dir/upper/round077/manuscript_notes"
cp cube-dim-7-bounds/paper/main.tex "$paper_check_dir/paper/main.tex"
cp cube-dim-7-bounds/upper/round077/manuscript_notes/verify_transcription.py "$paper_check_dir/upper/round077/manuscript_notes/verify_transcription.py"
python3 "$paper_check_dir/upper/round077/manuscript_notes/verify_transcription.py"
```

较早28页论文可通过已发布 Git commit f26db46db9e320aaffca74807ad25a815caad237 恢复；本稿源身份、审核与实际运行记录见 MANIFEST.json 和 upper/round077/REPORT.md。
