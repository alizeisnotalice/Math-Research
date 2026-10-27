# 中心立方体极大算子的一个显式弱型估计

中文论文，8页。主定理给一般维数的cube-cap显式界，七维应用为C7<24.554960525965。一般Lebesgue L1、中心轴平行、共同半径、全部正尺度。论文只使用已经通过二次独立审核的主线；原下界工作和研究goal保持暂停。

## 复现

在本目录运行：

    python3 verify_constant.py
    latexmk -norc -xelatex -interaction=nonstopmode -halt-on-error main.tex

需要Python3标准库和含ctex/Fandol字体的TeX Live或MacTeX。verify_constant.py用纯整数/有理数证明一个确切参数的常数比较，不需要NumPy、MATLAB、Arb或任何旧runner。constant_certificate.json是实际输出，VERIFY_RECEIPT记录代码hash与本次运行。

main.tex为可编辑全文，main.pdf已在本机TeX Live 2026/XeLaTeX成功编译；BUILD.json记录编译、源hash、8页Ghostscript渲染检查。内置编辑器在线编译曾因sandbox Perl依赖及Tectonic网络下载超时失败；不影响已成功的本机PDF，可直接打开main.pdf。

数学/数值/文献二审详上级second_review/REPORT.md和各自独立材料。manuscript_review.md是正文转写检查；其可选题名拼写建议未采用，因为SS原刊及正式题录确实印作functios，本文按原刊保留。主文不包含内部轮次或审核背书。

研究保持暂停，不以本稿声称确定最优Cn、连续参数最优或公开文献首创。
