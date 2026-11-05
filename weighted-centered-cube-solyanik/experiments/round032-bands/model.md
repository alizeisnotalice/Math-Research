# 连续薄孔带模型：从实际中心见证到系数 \(J\)

本附件研究 \(\mathbb R^2\) 的平坦权 \(w(x,y)=|x|^{-\gamma}\)，其中
\(0<\gamma<1\)，以及有限条按 \(y\) 方向排列的带。它是寻找近单位下界
构造的有限模型；它不提供完整 halo 上界。附件中的 \(\gamma=1/20\) 两带输入
给出一个严格的 \(J>1\) 下界例，而原 \(\gamma=1/2\) 数值搜索仍只是受限模型的
阴性证据。

令 \(B_i=[\ell_i,u_i]\) 按序排列，只有端点接触而没有正长度重叠；取
\(r_i>0\)、\(H_i\geq0\)，写
\[
 H(s)=\sum_iH_i\mathbf 1_{B_i}(s),\qquad
 H_t=\{(x,y):y\in B_i,\ |x|<tH_i\text{ for some }i\}.
\]
将第 \(i\) 条带的候选输入矩形取为
\[
 D_{i,t}=[-r_i-tH_i,r_i+tH_i]\times(B_i+[-r_i,r_i]),
 \qquad D_t=\bigcup_iD_{i,t},\qquad E_t=D_t\setminus H_t.
\]
假设对每个 \(i\) 及 \(y\in B_i\)，有严格余量约束
\[
 G_i(y):=\int_{y-r_i}^{y+r_i}H(s)\,ds
 \leq 2(1-\eta)r_i^2,\qquad 0<\eta<1. \tag{1}
\]

## 实际的中心方体见证

若 \((x,y)\in H_t\) 属于第 \(i\) 条带，则以该点为中心、半径 \(r_i\) 的
方体 \(Q=(x,y)+[-r_i,r_i]^2\) 包含在 \(D_{i,t}\subset D_t\)。其中横向包含性
来自 \(|x|<tH_i\)，并非将一个固定方体错误地当作所有孔点的中心方体。
该方体中所有孔的面积至多
\[
 2tG_i(y)\leq 4t(1-\eta)r_i^2.
\]
故 \(|Q\cap E_t|/|Q|\geq1-t(1-\eta)>1-t\)，每个孔点均在严格 halo 中。

函数 \(G_i\) 连续且分段仿射。在 \(B_i\) 上，它的所有折点都来自
\(\ell_j\pm r_i,u_j\pm r_i\)，再加 \(\ell_i,u_i\)。所以验证 (1) 只需对有限
节点检查。这正是独立 `band_fraction_verifier.py` 所重新生成并用有理算术
验证的列表。

## 连续极限与 \(J\)

令
\[
 R(y)=\max_i r_i\mathbf 1_{B_i+[-r_i,r_i]}(y),\qquad \beta=1-\gamma.
\]
当 \(t\downarrow0\)，并矩形的截面半宽趋于 \(R(y)\)，因而
\[
 \mu_w(D_t)\longrightarrow \frac2\beta\int R(y)^\beta\,dy,
 \qquad
 \mu_w(H_t)=\frac2\beta t^\beta\sum_i|B_i|H_i^\beta.
\]
由上面的严格中心见证，halo 包含 \(D_t\)（模零集）。于是本模型至多给出下界
\[
 C^{\rm out}_{c,w}(1-t)-1
 \geq \frac{\mu_w(H_t)}{\mu_w(E_t)}
 =Jt^\beta+o(t^\beta),
 \quad
 J=\frac{\sum_i|B_i|H_i^\beta}{\int R(y)^\beta\,dy}. \tag{2}
\]
因此只有一个经严格余量核验的 \(J>1\) 才会反驳候选首系数 1；本轮没有得到
这种构造。反向地，若有限搜索得到 \(J<1\)，它绝不构成完整 \(C^{\rm out}\)
的上界。

## 有理公共半径/高度实例的精确幂证书

对于有理 \(\beta=1-\gamma=p/q\)（既约）且所有 \(r_i=r\)、\(H_i=H\)，设
\(L=\sum_i|B_i|\)，并令 \(U\) 为区间 \(B_i+[-r,r]\) 的并集长度。此时
\[
 J^q=\left(\frac LU\right)^q\left(\frac Hr\right)^p. \tag{3}
\]
这只适用于该同半径、同高度的分支；一般情形并未由小数输出提升为精确比较。
`band_rational_inputs.json` 中两个 \(\gamma=1/2\) 搜索后有理化输入由 (3) 得到
\[
 J^2=\frac{6272}{6561}<1,\qquad
 J^2=\frac{1197900}{1216609}<1.
\]
它们同时满足 (1) 的严格余量。另一方面，\(\gamma=1/20\)、
\(B_1=[0,3/38]\)、\(B_2=[41/38,22/19]\)、\(r=1\)、\(H=25\)、
\(\eta=1/76\) 的解析构造满足 (1)，且
\[
 J^{20}=\frac{5^{18}}{2^{40}}>2,
\]
所以 \(J>2^{1/20}>1\)。其每个孔点的实际半径 1 中心方体见证仍由上节证明；
该正例并非原 \(\gamma=1/2\) DE 搜索所发现或覆盖。
