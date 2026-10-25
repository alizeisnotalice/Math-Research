# 奇 Taylor 截断为何是 cosine 的全实下界

令
\[
 P_m(x)=\sum_{k=0}^m(-1)^k\frac{x^{2k}}{(2k)!},\qquad m\text{ 为奇数}.
\]
记 \(R_m=\cos-P_m\)。在 \(m=1\) 时
\[
 R_1(x)=\frac{x^2}{2}+\cos x-1
       =\int_0^x(x-t)(1-\cos t)\,dt\ge0\quad(x\ge0).
\]
对任意奇 \(m\ge3\)，直接比较有限几何级数的 Laplace 变换，得到
\[
 \mathcal L R_m(s)=\frac1{s^{2m+1}(s^2+1)}
 =\mathcal L R_1(s)\,s^{-(2m-2)}.
\]
所以
\[
 R_m(x)=\frac1{(2m-3)!}\int_0^xR_1(t)(x-t)^{2m-3}\,dt\ge0
 \quad(x\ge0).
\]
两边均为偶函数，故对所有实数 \(x\) 成立。这个证明只用于本轮的 \(m=3,5\) 尾积分下界。
