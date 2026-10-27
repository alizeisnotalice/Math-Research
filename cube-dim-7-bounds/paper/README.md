# 中心立方体极大算子的覆盖估计与一般高维弱型上界

R048中文修订稿：main.tex与30页main.pdf。结构为引言、预备知识、主要证明、应用。

新主结果：每个固定0<d≤2^-24有自己的未知N(d)，进入该域后Cn≤F(A(t_d),B(t_d)−nd/65536)≤T_n−nd/131072；固定端点d=2^-24给Cn≤T_n−2^-41n。该新算子起点没有显式认证。R042的log n≥2^108完整显式定理仍保留，各结果按自己的适用域取min。

正文完整给有界逆列权、准确一维J、baseline均值及log方差、a.e.导数与AC连续匹配桥、inverse-row概率截断、双方度数重算、同历史有限a、一般L¹与全正共同半径、完整scalar代价。L≥2^240只用于完整cost标量排序，不是新算子起点。R047数值不作定理证据。

在仓库根目录编译并重复至引用稳定：

    mkdir -p /tmp/cube-cover-manuscript-build
    /Library/TeX/texbin/xelatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp/cube-cover-manuscript-build cube-dim-7-bounds/paper/main.tex

新转写检查只需SymPy（本机1.14.0）及标准库：

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round048/manuscript_notes/verify_transcription.py

实际13项符号/有理检查通过。没有重跑旧baseline、七维、生产程序或优化，没有使用共享MATLAB。编译与排版检查见BUILD.json，数学新全文独审已ACCEPT，见audit/round048/manuscript/review.md；主代理另行实际编译并渲染检查，证据分别记录。

新来源hash、转写清单与归档白名单在upper/round048/manuscript_notes。旧94c5稿身份与private恢复源已记录，不重复归档旧PDF或阶段全文草稿。历史七维verifier、证书和R041旧notes本轮未修改；它们不作为新加权定理的证据。公开引文沿用已核Stein–Strömberg1983、Aldaz2017v2及Aldaz2020v2，不新增未核文献。

推荐数学venv命令是可复核指引，作者实际运行依赖身份在原receipt中保留，未把修改后的命令记为新数学重跑。原24页稿可用upper/round048/source_copy_manifest.json中的稳定Git blob恢复；不重复保存旧PDF。
