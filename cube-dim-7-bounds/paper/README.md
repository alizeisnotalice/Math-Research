# 中心立方体极大算子的覆盖估计与一般高维弱型上界

R054中文修订：main.tex及20页main.pdf。最新主定理为自然log n≥2^27时

    C_n ≤ Φ_simple < n log n − n log log n − 2^-16 n。

raw准确组在log n≥2^32也满足相同display，且完整根式严格小于Φ_simple。正文保留两组完整gamma、sigma、t与Φ定义，逐项量化operator门槛，并自足证明完整函数排序。

主线为R051准确泛budget、exact inverseJ/pinned face/AC连续匹配、有界source权的双方degree、同history actualμ、finite dilation及general L¹/all-positive-common-radius桥。geometry、probability、cap、tau、scalar的全部finite guard与Taylor尾均在正文给出。有限搜索只提示两组输入，最终proof不依赖全遍历最优性。

较旧R042/R036高维长证明移出本轮正文，来源与旧稿可恢复身份保留；R034较小L门槛的完整family与R037 pointwise范围仍在应用章。

仓库根编译（重复至引用稳定）：

    mkdir -p /tmp/cube-cover-manuscript-build
    /Library/TeX/texbin/xelatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp/cube-cover-manuscript-build cube-dim-7-bounds/paper/main.tex

两准确tuple的独立有限常数转写检查（仅Python标准库）：

    python3 cube-dim-7-bounds/upper/round054/manuscript_notes/verify_budget_and_threshold.py

本机实际两组各24项及3项比较全部通过。没有旧baseline/n7/production/参数搜索/MATLAB/Git同步。编译与排版见BUILD.json；新源已完成独立全文数学审核及唯一K定义补全后的增量锁定。

R048已发布源SHA256 ecd7fd7487083b7eb6318a0a94f969f9d291b9d4df443cc3374497a7f103b248，主提交5122f5eba011430a857692896772ddffd02fba5d。其private恢复副本不进入新归档，不重复旧大PDF。历史七维verifier、证书及其他轮notes本次未改，不作为新主定理依据。公开引用沿用已核原始版本，无新优先性或最优性声明。

a05全文数学独审后，按唯一P3建议显式补齐所用K/P与对应row/col定义。最小差异见R054notes/kernel_definition_delta.patch；原a05全文留private，不双份归档。未改估计/门槛，51项准确检查未重跑；只重新编译与看排版，最终增量锁审已返回full ACCEPT，无待修事项。

主代理正式目录独立编译、全部页渲染/排版证据见BUILD.json；作者实际编译单列upper/round054/author_build.json。数学源17e3保持字节不变。建议venv命令不是新数学运行记录。rebase后以source_copy_manifest的稳定blob恢复旧稿，作者当时的原commit身份在copied来源中保留。
