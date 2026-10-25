# 文献表（2026-10-25；持续补核）

|作者/原文|版本、定位|对象和范围|常数、证据状态|
|---|---|---|---|
|J. M. Aldaz, [A remark on the centered n-dimensional Hardy–Littlewood maximal function](https://dml.cz/bitstream/handle/10338.dmlcz/127554/CzechMathJ_50-2000-1_14.pdf)|CMJ 50 (2000),103–112；Prop.1.4、Lemma1.1|中心轴平行方体；Lebesgue；一般L1，经原子传递；全部正尺度；n≥2|A_n公式；原文纸本110–111页公式及证明已读，根区间计算通过|
|J. M. Aldaz, [Kissing numbers and the centered maximal operator](https://arxiv.org/html/1811.10478v2)|v2,2020-03-10；Th.1.2、Cor.3.5；版本页已核|中心∞范数球；所有局部有限Borel背景的一致界；各背景L1；全部尺度；所有n|2^n；已阅读相应定理。可给Lebesgue上界，锐性不可移作Lebesgue最优性|
|J. M. Aldaz, [The weak type (1,1) bounds for the maximal function associated to cubes grow to infinity with the dimension](https://annals.math.princeton.edu/wp-content/uploads/annals-v173-n2-p10-p.pdf)|Annals173(2011),1013–1023；第1013页Theorem后，引AV07 Th.2|中心轴平行方体；Lebesgue；一般L1；全部尺度|C_(n+1)≥C_n，独立传递证明另附；渐近发散不生成固定维新常数|
|Tejasvi Singh Tomar, [Bounds for the weak type (1,1) constant of the centred maximal operator over squares](https://doi.org/10.5281/zenodo.23554628)|Zenodo version 4，2026-10-22；Th.1.3 p.2；Lemma6.1/6.2 p.10；证明pp.10–11；原PDF已读相关页，哈希见tomar-v4-source.json；公开代码固定bca281a…|中心轴平行二维方体；Lebesgue；一般L1；全部正尺度；周期原子输入经严格半径余量和截断平滑传递|c₂≥1.7027237335；准确公开分数见results/public_lower_run.json，本机整数核验已通过。原PDF版本与定理编号缺项已补齐；未重建Lean、未独立核验其上界|
|Yongxi Lin, [Two partial balayage principles for weak-type estimates](https://coolrmal.github.io/articles/two-partial-balayage-principles/)|网页未标版本；2026-10-25访问；§3 Principle2|Lebesgue卷积；一般L1及全尺度；需齐次对称保守Markov生成元与障碍分解|方法资料；本轮只审读假设，不登记新上界；核必须偶，生成元符号非正|

令u=2^(1/n)、a=(1+u)/2：
A_n=(3/2)a^(n−1)−(1/2)[(3a−2(3/2)^(1/n))^(n−1)−((3−u)/2)^(n−1)]。
三维公开比较采用Tomar准确分数（优于A3），四维采用A4（优于该分数）；所有当前内部下界均更大。另以compare_public.py根有理夹逼完成准确比较，结果见public_comparison.json。
