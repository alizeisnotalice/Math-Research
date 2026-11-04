# 中心方体 Tauberian 常数的高阶端点研究

研究Lebesgue测度下的中心轴平行立方体极大算子，使用严格水平集。
线性输入为任意有限正测度可测集；共享方体的多线性常数以公共交集归一化，单个输入允许无限测度。

## 当前严格结论

置a=2^nt、v=2t^(1/n)，令
\[
L_n(t)=\left[1-\frac{a}{(1+v)^n}\right]^{-1}.
\]
对每个固定整数n≥2，真实常数满足
\[
C_n(1-t)=L_n(t)+O_n(av^6).
\]
因此直至t^(1+5/n)的全部系数均确定：
\[
C_n(1-t)=1+\sum_{j=0}^{5}b_{n,j}t^{1+j/n}+O_n(t^{1+6/n}),
\quad b_{n,j}=2^{n+j}\sum_{k=1}^{1+\lfloor j/n\rfloor}
(-1)^{j-(k-1)n}\binom{n+j-1}{j-(k-1)n}.
\]
通式包含低维幂次碰撞；三维t^(8/3)的系数为零，不能将所有六个幂次都称为非零项。
上下界覆盖全部允许输入和全部半径；常数及起效阈值依赖固定n，没有维数一致主张。
av⁶尺度的归一化极限及精确系数尚未确定。

此外，对每个固定n≥2，一般近单位矩阵晶格配合孔形与质量重配的实际有限输入给出
\[
\liminf_{t\downarrow0}
\frac{C_n(1-t)-L_n(t)}{2^{n+6}t^{1+6/n}}\ge\binom n2\frac{n^2}{27}>0.
\]
因此L_n在充分小参数下严格低于真实常数。对所有固定n≥2，与全部输入上界合得C_n−L_n与t^(1+6/n)同阶；尚未确定该尺度的极限或精确系数。

多线性超临界区精确约化为
\[
C_n^{(m)}(\lambda)=C_n(1-m(1-\lambda^{1/m})).
\]
阈值为(1−1/m)^m；m≥2时在阈值及以下发散。
令A_m(ε)=1+(m−1)ε/(2m)+(m−1)(2m−1)ε²/(6m²)。将t换成ε=1−λ后，关于全部m≥1一致有
\[
C_n^{(m)}(1-\epsilon)=1+\sum_{j=0}^{5}b_{n,j}[\epsilon A_m(\epsilon)]^{1+j/n}
+O_n(\epsilon^{1+6/n}).
\]
三维ε^(8/3)系数为160(1−1/m)，说明参数替换会填入单线性中消失的幂次。

## 材料与历史

- `paper/tauberian.pdf`：实际编译、检查的32页论文；`main.tex`、`sections/global.tex`、`sections/skew.tex`和`sources.bib`为可编辑项目。
- `rounds/015.md`、`reviews/015-*.md`：一般矩阵构造、质量偏移、联合下界及有限输入审查；`rounds/014.md`、`reviews/014-*.md`：局部dyadic小源回收、精确E+M结算及全维数双边幂阶；第13轮记录孔形重配下界，第12轮记录此前高维上界。
- `literature/2026-11-05.md`、`sources.md`：本日最新针对性文献核对及模型比较，不宣称穷尽或新颖性。
- `proofs/`：仍能排除错误路线的反例和受限构造排除；`verification/`为必要证书及探索复现材料。
- 第4轮闭合二阶；第6轮取得幂余项；第7轮证明固定异常尾桥反例；第8轮以异常回收闭合第三项；第9轮进一步取得n≥3的第四项；第10轮证明周期包络在更细尺度并不精确；第11轮以全对称晶格增强其正差下界；第12轮用全源约束与免费面缝将高维余项推进至av⁶；第13轮以孔形重配给av⁶正下界，确定n≥5的首差幂阶；第14轮消除低维损失，推广到所有固定n≥2；第15轮将正下限增强为所有坐标对的联合值。

正文只保留当前证明链。被覆盖的空间缓冲和双深度主路线可从历史提交恢复；异常尾反例仍有效，新证明没有假设固定阈值尾o(vS)。没有把周期包络的任意阶展开直接认作真实常数展开。

## 复现

```sh
cd paper
latexmk -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

需要XeLaTeX、BibTeX、Songti SC和Times New Roman。主要符号工具为本机Mathematica/Wolfram15 ARM：

```sh
wolframscript -file verification/wolfram/augmented-corners.wl
wolframscript -file verification/wolfram/multidepth.wl
wolframscript -file verification/wolfram/fourth-transfer.wl
wolframscript -file verification/wolfram/young-face.wl
wolframscript -file verification/wolfram/endpoint-six.wl
```

角部和系数批次实际各启动15个子内核；多深度脚本实际18项恒等式/阶次通过。有限维系数与标量检查只承担各自范围，一般维数几何、全部输入和全部重数的一致性由论文解析证明承担。

第10轮最小计算包为`verification/skew-lattice/`：204个有理根括界及向外区间，182个正G、22个负G全部保留；15个明确有限输入和30个实际见证平均值也已精确复算。其README含完整复现命令。数值与证书的覆盖范围由正文一般证明及相应守卫解释，不从有限参数外推。

第11轮当前计算包为`verification/symmetric-skew/`：105组全对称晶格有理证书（89正、16负），15个有限输入和30个严格见证；Mathematica以15子内核独立交叉核验全部105组，并核验行列式和优化恒等式。一般维数证明不依赖有限枚举。第10轮最小包为历史结论保留。

第12轮新增两份短Wolfram脚本及实际输出：Young余项六项精确检查全通过；系数核验使用15子内核检查n=2,...,16并核验全部m的标量表达式。它们只核验代数，一般维几何和全部输入上界由正文证明及独审承担。

第13轮当前最小包为`verification/reshaped-skew/`：120组筛选中90预算可行（87正差、3负差）、29预算不足、1守卫失败；另15个n=2,...,16的有限输入均精确认证。Mathematica各以15子内核交叉核验119个守卫合格筛选组及15个新有限输入。中心平均记录为严格下界，未冒称等式；135条存档根括界重建通过。

第14轮当前最小符号包为`verification/sixth-order/`：最终清洁Wolfram脚本实际使用15子内核核验n=2,...,16、全部m的物理参数回归及新账本恒等式，FAILURES=0。正文以一般系数证明和统一四阶参数余量处理所有维数/重数，有限符号表只作核验。经典强矩形长证明已由独审通过的局部dyadic计数替代。

第15轮最小包为 `verification/matrix-skew/`（22个文件）：Mathematica15子核处理不等互边精确积分及矩阵行列式；Python15进程给出n=2,...,16的实际有限输入有理证书，另用Mathematica15子核独立交叉检查。二维子包通过两个代数根明确指定实际切片孔，给出指定中心见证的实际平均等式和严格裕量；几何已由非作者Sol审查。完整复现命令和边界、预算、精度范围见各README。

该轮的变分最优值仅针对固定矩阵方向、u³扰动及明确的共同见证安全选择类。它给出真实常数的更强下界；类内上界不构成全部可测输入的最优常数上界。归一化极限与精确系数仍是待解目标。

第16轮终止纯q≥3有向环的指定安全重配路线：独立短证明给出其正收益为o(av⁶)，不缩小全输入精确系数缺口。只保留`rounds/016.md`、`proofs/016-directed-cycle-exclusion.md`和三个最小符号核验文件；不向正文加入更晚项的长篇模型优化。

第17轮证明全输入临界分支的加权盒形/面半径紧性，并以真实有限稀有相混合反例排除“每胞一致近整格”桥。`proofs/017-critical-compactness.md`及`reviews/017-compactness.md`保留完整有效证明和LT可见性修补；LT接触计数在本轮仍有条件，原输入⇒LT及真实锐系数未在本轮证明。五个最小计算文件核验代数和横截面。主正文/PDF不变。
