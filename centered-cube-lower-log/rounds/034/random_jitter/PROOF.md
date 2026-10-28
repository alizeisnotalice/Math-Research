# 独立格内抖动的统一自相关核

本命题在本机FFT/Fourier探路之后建立并独立审查。它给出另一种一维来源的二点核，不依赖有限域或Sidon标签。它没有证明新的极大算子下界阶，也没有把完整旧概率证明自动迁移。

## 1. 对象与结论

取整数K≥2，令U₀,…,Uₖ₋₁独立均匀分布于[0,1]，并置a_j=(j+U_j)/K。在圆周上定义周期分格原函数
\[
 P_U((j+r)/K)=\mathbf1_{\{U_j\le r\}}-r,\qquad 0\le r<1.
\]
有 \(\|P_U\|_\infty\le1\)、\(\operatorname{TV}(P_U)\le2K\)，分布导数为 \(\sum_j\delta_{a_j}-K\)。令
\[
 R_U(t)=\int_{\mathbb T}P_U(x)P_U(x+t)dx,\quad
 \Phi_K(t)=\tfrac16(1-K\operatorname{dist}_{\mathbb T}(t,0))_+^3,
\]
\[
 \varepsilon_K=4\sqrt{\frac{\log K}{K}}+\frac3K.
\]
则
\[
 ER_U(t)=\Phi_K(t),\qquad
 P\{\sup_t|R_U(t)-\Phi_K(t)|>\varepsilon_K\}\le2K^{-6}.    \tag{1}
\]
固定A>0。在一个与(1)相交后概率至少 \(1/2-2K^{-6}>0\) 的事件上，还可同时使所有M<K、\(|c_m|\le A/m^2\) 的低频补偿满足下述统一交叉帽。因而对每个K都存在一份确定抖动来源满足全部这些界。

## 2. 期望核的准确形状

对均匀U，\(E[\mathbf1_{U\le r}-r]=0\)，且
\[
 E[(\mathbf1_{U\le r}-r)(\mathbf1_{U\le s}-s)]
       =\min(r,s)-rs.
\]
不同格子的U独立，乘积期望为0。若 \(0\le Kt=\tau\le1\)，只有不跨格的部分贡献：
\[
 ER_U(t)=\int_0^{1-\tau}r(1-\tau-r)dr=(1-\tau)^3/6.
\]
在圆周另一端由偶性得到同样公式；距0至少1/K时所有贡献来自不同格子，期望为0。K≥2排除了单格自环的特殊情形。这给出(1)的期望恒等式。

## 3. 有界差分与全连续位移

改变一个U_j只在第j格内改变P，改变量绝对值至多1，故 \(\|\Delta P\|_1\le1/K\)。因P和修改后的P均绝对值不超过1，对每个固定t有 \(|\Delta R_U(t)|\le2/K\)。

为展开所用概率界：均值为0、取值范围长度c的变量X，其对数矩母函数的二阶导数是倾斜方差，至多c²/4。由零点处函数值及一阶导数均为0，得到 \(Ee^{sX}\le e^{s^2c^2/8}\)。对依次揭示独立U_j的Doob差分应用此条件界，得到
\[
 P\{|R_U(t)-ER_U(t)|>a\}\le2e^{-Ka^2/2}.
\]
取K²个等距位移及 \(a=4\sqrt{\log K/K}\)，并集失败率至多2K⁻⁶。

另一方面，BV平移不等式给
\[
 |R_U(t+h)-R_U(t)|\le\|P_U\|_\infty
         \|P_U(\cdot+h)-P_U\|_1\le2K|h|.
\]
\(\Phi_K\) 的Lipschitz常数至多K/2。从网格补到全部圆周只需额外3/K，证明(1)。这里没有对随机跳点作未经许可的并集。

## 4. 不依赖标签的全低频交叉控制

对1≤m<K，置 \(Z_m=\sum_je^{-2\pi im a_j}\)。根单位求和给EZ_m=0。独立性和 \(|e^{iz}-e^{iw}|\le|z-w|\) 给
\[
 E|Z_m|^2\le K(\pi m/K)^2=\pi^2m^2/K.
\]
非零Fourier系数满足 \(\widehat P_U(m)=Z_m/(2\pi im)\)，故 \(E|\widehat P_U(m)|\le1/(2\sqrt K)\)。

令 \(H_L(x)=P_U(x+L/2)-P_U(x-L/2)\)，\(B_L=\sum_{m\le M}c_m\sin(\pi mL)\cos(2\pi mx)\)。对全部L,J，
\[
 |\langle H_L,B_J\rangle+\langle B_L,H_J\rangle|
 \le4A\sum_{m=1}^{K-1}\frac{|\widehat P_U(m)|}{m^2}=:X_A,
 \quad EX_A\le\frac{A\pi^2}{3\sqrt K}.
\]
Markov给 \(P\{X_A>\Gamma_K\}\le1/2\)，其中 \(\Gamma_K=2A\pi^2/(3\sqrt K)\)。与(1)相交得到前述正概率。这个事件覆盖所有M<K与整个系数包络，不是只覆盖实验中计算的16个频率。A=0时交叉项恒为0，单独处理即可。

箱平滑只乘绝对值不超过1的Fourier因子。若两个来源另作不同平移，附加的Fourier相位也不改变此帽。

## 5. 平滑、完整核及原常数核对

令 \(P_{U,\delta}=P_U*\beta_\delta\)、\(R_{U,\delta}=R_U*(\beta_\delta*\beta_\delta)\)。三角概率核的 \(E|z|=\delta/3\)，而 \(\Phi_K\) 的Lipschitz常数为K/2。因此good事件上
\[
 \|R_{U,\delta}-\Phi_K\|_\infty
       \le\varepsilon_K+K\delta/6.
\]
准确展开四项有
\[
 \langle H_{L,\delta},H_{J,\delta}\rangle
       =2R_{U,\delta}((L-J)/2)-2R_{U,\delta}((L+J)/2).
\]
在I=[10/9,11/9]、K≥10时，后一位移的理想Φ为0。故对全部L,J∈I，完整核满足
\[
 \left|\langle H_{L,\delta}+B_L,H_{J,\delta}+B_J\rangle
 -\left[\tfrac13(1-K|L-J|/2)_+^3+\langle B_L,B_J\rangle\right]\right|
 \le4\varepsilon_K+\Gamma_K+\tfrac23K\delta.                \tag{2}
\]
这是全连续窗口的二点核，不需删除有限域缺陷lag。

现行原归档 `../../028/source_interface_baseline.tex` 第四、五步已经使用cubic核，而不是triangular核。本式与该局部形状相同。设未平滑极限中 \(r=Kv|L-J|/2\)，则高高核的v倍损失趋于r；归一化方差 \(\gamma L^2/3=9L^2/29\) 后，\(v(1-\kappa)\to29r/(9L^2)\)。采用 \(\sqrt\lambda B_r-\lambda r/2\) 的记号时，\(\lambda=58/(9L^2)\)，与原记号一致。固定w平滑仍须保留其局部卷积轮廓，不能直接把未平滑导数当成固定δ导数。

## 6. 可以迁移的预算与没有证明的下界

对原系数c_m，\(1+\sum_{m\ge1}b_m\sin^2(\pi mL)=\gamma L^2\) 在I上成立。取M≈√v、\(\delta=(wKv)^{-1}\)，则(2)的随机核误差乘v，在
\[
 v^2\log K/K\longrightarrow0
\]
时趋零，低频尾vT_M也趋零，固定w的平滑预算保持O(1/w)。因此独立随机抖动提供一个没有field存在性门槛的来源二点替代。

若 \(v\asymp n^\beta\)、\(K\asymp n^{(1-\beta)/2}\)，这个来源二点条件支持 \(0<\beta<1/5\)。这比原q>H(B)≈B⁴的充分构造所给β<1/9更宽。这里只移除了该来源核障碍；它没有自动证明完整倾斜、极限、pair求和和L1截断链，也没有证明任何 \(n^{1/2+\eta}\) 或正对数下界。单来源首阶 \(K\sqrt v\) 抵消仍然存在。

实验仅计算固定seed下K15、63、255、1023的FFT网格自相关和至多16个低频的实际Fourier和。FFT误差、有限网格与BV积分误差没有区间认证，不能据此宣称这几个实现已获全部good事件证书。上述存在性和连续量词来自解析证明，不是数值拟合。
