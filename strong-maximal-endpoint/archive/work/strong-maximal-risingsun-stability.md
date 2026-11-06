# 一维 rising-sun 稳定性与坐标复合缺口

日期：2026-11-06。仅处理一维非中心最大算子及其逐坐标复合；不声称高维 Steiner 重排结论。

令 \(M^\pm\) 为左右单侧最大算子，\(H^\pm=\{M^\pm f>s\}\)，
\(H=H^+\cup H^-\)，\(I=H^+\cap H^-\)。rising-sun 等式
\(s|H^\pm|=\int_{H^\pm}f\) 与容斥给出
\[
 s|H|=\int_Hf+\int_I(f-s).
\]
若 \(D=|H|\)、\(F(u)=\int_0^u f^*(v)dv\)、\(A_s=\int(f-s)_+\)，则精确缺口为
\[
 F(D)+A_s-sD=
 \bigl(F(D)-\int_Hf\bigr)+
 \int_{I^c}(f-s)_++\int_I(s-f)_+\ge0.
\]
第一项是质量未向 halo 集中造成的重排损失，后两项是两个单侧 halo 的交集错配。偶对称递减重排使三项均为零，因此在固定分布函数下最大化一维非中心 halo 测度。

对 \(p\ge2\)、\(\|f\|_p=1\)，设 \(c=\|Mf\|_p\)、\(b=\min(M^+f,M^-f)\)、
\(H=c^{p-1}-\int f(Mf)^{p-1}\)。两侧幂等式相加给
\[
 p c^{p-1}+1-(p-1)c^p
 =pH+\int\{(p-1)b^p-pfb^{p-1}+f^p\}.
\]
因 \(b\ge f\)，其积分被
\[
 {p(p-1)\over2}\int f^{p-2}(b-f)^2
\]
下界。令 \(c_p\) 为锐一维范数根，逐坐标正规化 \(F_i=f_{i-1}/\|f_{i-1}\|_p\)、
\(a_i=\|M_iF_i\|_p\) 后，若 \(D_i\) 是上式右端，则
\[
 \log{c_p\over a_i}\ge
 {D_i\over p(p-1)c_p^{p-1}(c_p-1)}.
\]
因此
\[
 n\log c_p-\log{\|M_Sf\|_p\over\|f\|_p}
 \ge\sum_i{D_i\over p(p-1)c_p^{p-1}(c_p-1)}+
 \log{\|M_n\cdots M_1f\|_p\over\|M_Sf\|_p}.
\]
右侧两部分均非负，分别量化每个坐标的 rising-sun 近饱和与“逐坐标适配复合”相对于真实强矩形的几何松弛。

MATLAB R2026a 对标量余项的二次下界在 \(p=2,4,8\)、\(b/f\in[1,6]\) 上核验，最小数值松弛为 0（\(p=2\) 的浮点误差为 \(-3.6\times10^{-15}\)）。解析证明来自二阶导数下界，不依赖该采样。
