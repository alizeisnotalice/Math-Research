# 中心立方体极大算子的显式弱型估计及其维数渐近

扩展修订中文论文，13页。旧8页稿可从Git提交f51eea16恢复。原覆盖定理从头复核通过，新增U_n的强化渐近、前三维端点与n≥4内部全局极小值存在性、维数传递，以及21维准确数值与有限a证书。七维上界没有降低；不声称候选全局最优、真实Cn渐近等价或全球首创。

复现（在项目根cube-dim-7-bounds）：

    python3 revision_v3/numeric/verify.py
    python3 revision_v3/compare_bounds.py
    cd paper
    python3 verify_constant.py
    latexmk -norc -xelatex -interaction=nonstopmode -halt-on-error main.tex

前两个验证器仅用Python标准库；可选SymPy符号检查在revision_v3/asymptotic和endpoint。generate.py需要SciPy，只用于候选搜索，不是验证证书所必需。numeric/table.json保留B、F和有限a/tau全部准确端点；dimension_comparison.csv用于阅读。dimensions_table.tex是表格生成快照，最终main.tex已内联相同21行以避免booktabs输入展开问题。

一般Lebesgue L1、中心轴平行、共同半径、全部正尺度。论文保持四章结构；补充报告、来源及新增独立审稿见revision_v3/REPORT.md与final_audit/。BUILD.json记录实际本机XeLaTeX编译及13页Ghostscript逐页检查；没有将历史运行当本轮计算。原七维verify_constant.py及旧结果收据保留不变，新审稿实际复核运行记录见final_audit/review.md。

本次一次性授权扩展已结束，旧长期goal工具状态仍为paused；下界及其他项目未恢复，不自动启动下一批。
