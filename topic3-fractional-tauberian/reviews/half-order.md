# 半阶固定凸体首项与 Petty 极值的独立复核

## 结论

所给半阶论证正确。对中心对称满维凸体 \(K=-K\)、\(0\in\operatorname{int}K\) 及
\(\beta=n/2\)，固定测试集合 \(E=K\) 的完整切向首项确为
\[
 \lim_{t\downarrow0}\frac{|H^K_{n/2,\mathrm c}(1-t,K)|}{|K|t^n}
 =|K|^{n-1}|\Pi^*K|,
\]
其中 \(\Pi K\) 使用本文已经固定的经典投影体归一化
\(h_{\Pi K}(u)=|u|\,|P_{u^\perp}K|_{n-1}\)。对 \(n\ge2\)，Petty 投影不等式给出
\[
 |K|^{n-1}|\Pi^*K|
 \le\left(\frac{\omega_n}{\omega_{n-1}}\right)^n,
\]
且在中心对称类中等号当且仅当 \(K\) 是以原点为中心的椭球。

这只是“同时改变固定输入体和其锚定位似基”的固定模型极值；它不涉及固定某个基后
对所有可测 \(E\) 取上确界的 Tauberian sharp 常数。

## 半阶变分的核对

令 \(h(y)=\langle y,\nu_K(y)\rangle\)。因 \(0\in\operatorname{int}K\)，有
\(h>0\) 于 \(\mathcal H^{n-1}\)-几乎处处的边界点。把
\(w=\langle v+ay,\nu_K(y)\rangle=v\cdot\nu_K(y)+ah(y)\) 代入
\[
 D_K(v,a)=\int_{\partial K}w_-\,d\mathcal H^{n-1},
 \qquad \int_{\partial K}w\,d\mathcal H^{n-1}=na|K|,
\]
即得
\[
 D_K(v,a)=\frac12\left(\int_{\partial K}|v\cdot\nu_K(y)+ah(y)|
 \,d\mathcal H^{n-1}-na|K|\right).
\]
所以 \(\beta=n/2\) 时
\[
 \Phi_{n/2,K}(v,a)=\frac1{2|K|}
 \int_{\partial K}|v\cdot\nu_K(y)+ah(y)|\,d\mathcal H^{n-1}.
\]

中心对称性把几乎处处的边界点配为 \(y,-y\)，并有
\(\nu_K(-y)=-\nu_K(y)\) 与 \(h(-y)=h(y)\)。对
\(b=v\cdot\nu_K(y)\)，恒等式
\[
 |b+ah|+|-b+ah|=2\max\{|b|,|a|h\}\ge2|b|
\]
表明 \(a=0\) 是一个极小点。因此
\[
 N_{n/2,K}(v)=\frac1{2|K|}\int_{\partial K}|v\cdot\nu_K|,d\mathcal H^{n-1}
 =\frac{h_{\Pi K}(v)}{|K|}.
\]
这里最后一个等式正是 Cauchy 投影公式的经典归一化。应写“\(a=0\) 是极小点”，
而不应写成“唯一极小点”：例如多面体和某些方向可能允许一段非零 \(a\) 同样达到最小值。

从而
\[
 S_{n/2,K}=\{v:h_{\Pi K}(v)<|K|\}=\operatorname{int}(|K|\Pi^*K),
 \qquad |S_{n/2,K}|=|K|^n|\Pi^*K|,
\]
与完整固定模型切向定理相接后即得到首项公式。

## Petty 常数的归一化

已实际查阅 Lutwak--Yang--Zhang 的原始 PDF，第 4--5 页明确说明其
\(\Pi B=B\) 是不同于经典的归一化；同页给出归一化后的 Petty 右端为
\(\omega_n^n\)。本文定义满足经典的 \(\Pi B=\omega_{n-1}B\)，故极体体积相差
\(\omega_{n-1}^{-n}\)。因此本文应使用
\[
 |K|^{n-1}|\Pi^*K|
 \le\omega_n^n\omega_{n-1}^{-n}
 =\left(\frac{\omega_n}{\omega_{n-1}}\right)^n.
\]
原文还说明经典 Petty 不等式对凸体的等号体为椭球；限制到 \(K=-K\) 时即为原点中心
椭球。这个等号条件适用于 \(n\ge2\)。如保留 \(n=1\)，应单独处理，避免把高维
投影体记号和 Petty 定理不加说明地延至零维投影。

## 正文措辞建议

推荐把本结果命名为“半阶固定模型的仿射极值”，并在定理后紧接以下范围声明：

> 此处极值同时改变资料体 \(K\) 与锚定位似基，且仅对资料 \(E=K\) 的切向首项成立；
> 它不比较固定 \(K\) 下任意可测资料的 Tauberian 常数。

这能避免将 Petty 不等式误读为所缺的全体集合统一上界。
