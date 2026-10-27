# R048 作者转写核对

正文为paper/main.tex；本文件是作者自审，不代替新稿独立全文审。源已按数学写作规范保留引言、预备知识、主要证明、应用；前人A/B/C，本文结果数字编号。公开引文沿用已核原始版本，无新公开优先性或最优性主张。

## 新主线

- 固定0<d≤2^-24后统一选择N(d)，不允许d随n变化；端点d=2^-24全部参数重新代入，最终T−2^-41n。未知N(d)不由scalar L≥2^240替代。
- 原有限选择、固定实际增量、λa/μa支配与准确倍率保留。来源因子b(z)=min(1,theta/G)只用于后选择的核，不改变原算法权重或资格。
- normalized目标交长随目标r下降，完整column≤4nhG/R；b乘入核后column≤4k。行须按同一b重算，未保留旧raw energy再单独缩小列分母。
- 一维局部引理完整给精确J：稳定excess Psi、两侧m±、来源端点至多一端删去、ordinary保留长度及inverse excess均展开。
- Psi偏导和q导数分deep/tail控制；q map只有a.e.导数，capture onset可有kink，有限分支+有界one-sided slopes给AC，故可穿过kink积分。
- q0=r−d+h的baseline严格满足R=r；Fubini取消ell得到EJ=1。log均值与log方差均保留，不用小variance替代非零mean的控制。pinned target coordinate保留确定项。
- q0可能略出source band，但完整扩展Q仍满足local lemma域；good target只在自己的r面定义，之后AC给每个matched q的确定性乘积界。未用有限q网格或不可数随机事件交。
- 完整inverse壳系数保留pinned source ell_*与Ni/2面计数；自由坐标才积分成J。U在n/8与8n之间。
- clipping概率使用capture/G inverse-row law；其free coordinate密度normalizer J≥1/2，pin贡献≤d，E_inv[-logG]≤3Ld²，Markov tail≤1/4。没有用raw source概率替换。
- 重算soft角系数3ntheta/32至8ntheta；完整径向行3k/32至32k；b≤1后按未加权corner删除O(nd exp(−Ld/32))，正指数d/32−12d²>d/64确保cap=o(k)。最终row≥3k/64、row≤32k、col≤4k。
- same-history缺额用同一目标一次资格，weighted actual row≤tau无theta^-1因子。target实际半径<3、source>4及κ+2/(a+1)<αβ逐式给出。
- finite-a c_a环境倍率下energy乘c_a²、degree和目标质量乘c_a；缺额≥3nd/65536，所用nd/65536留严格余量。无测度浓度/selector convergence假设。
- 根号B>D与tau/k≤sqrt(L)/h exp(−4d²L)→0；固定n先a→∞，原紧严格水平集见证桥保持一般L¹/全部正共同半径。
- scalar x=(1−32d²)L：B≤T+64nd²+1，A+cross=o_d(n)，64d²≤(d/65536)/4。没有复用旧x的W_n或只比较负项。

## 比较与旧结果

- R042完整显式链保留：log n≥2^108时Cn≤Phi42<T−2^-32 n L^-1/3。旧R036及R034证明也保留，各按自己的域取min。
- 应用以Phi_{sigma,Delta}的完整上下包比较；actual根式预算Delta与display Delta/2区分。
- scalar L0=2^240处E_sigma的准确有理上包及单调powers覆盖整个尾部。端点Phi_*<Phi42和d=2^-30旧展示cost是完整公式排序，仅在新own-N内才兼有Cn意义。
- R037框架范围和fixed-level spike原准确限定保留；不把全weak norm声称趋1。R047局部数值、R043障碍或另一soft主链不作本文新定理证据。

## 本轮实际检查

- 仅执行verify_transcription.py的13项SymPy1.14.0/Fraction局部转写检查：Psi偏导、baseline primitive、精确J与分段A等价、bridge与endpoint及完整cost常数。无旧baseline/生产/优化/n7重跑，无共享MATLAB。
- XeLaTeX实际编译及最终字体/越界/引用检查另见paper/BUILD.json；编译不代替数学审核。
- 94c5原稿私有恢复源在private/main_94c5.tex；私有阶段草稿及旧PDF不进入R048归档白名单。
