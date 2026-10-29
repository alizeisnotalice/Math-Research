# 第052轮：真实cos幅度混合探路与单产品逐点支配的反例

本目录先实际运行 `probe.py`，再记录以下结论。最强一般严格结论是显式的likelihood曲率/Taylor式，以及非退化prior通常不能被某一个单幅度产品逐点支配。没有得到weak比较的反例，也没有把数值中的弱泛函比较当成定理。方差平坦尺度族的严格共同排序比较由root的同轮 `../flat_budget/PROOF.md` 另行给出，本稿不重复该大引理。

## 1. 真正来源、同一查询和同一prior

来源取

\[
\rho_\alpha(x)=1+.9\alpha\cos(2\pi m x),\qquad
m=1\text{ 或 }3,
\quad\alpha\in\{0,.2,.35,.5,.65,.8,.95,1.1\}.
\]

所有密度非负且质量1，最小source密度至少.01。中心区间**全边长**L平均准确是

\[
Y_{\alpha,L}(x)=1+\alpha d_L(x),\qquad
d_L(x)=.9\frac{\sin(\pi mL)}{\pi mL}\cos(2\pi m x).
\]

取三个L=10/9、7/6、11/9，n=16/64/128，n维输出为

\[
T_{\alpha,L}(x)=\prod_{i=1}^n(1+\alpha d_L(x_i)).
\]

同一全局幅度priorν在全部尺度上固定；先取产品后作凸混合：`M_L=∫T_{α,L}dν`，没有每坐标重新抽幅度。四个prior分别为正幅度均匀、两端点等权、高幅度重权、常数幅度0与高幅度等权。

每配置在同一32768个Haar查询上，比较 `max_L M_L`、每个α的 `max_L T_{α,L}`。实际按16批、每批2048查询运行，避免保存大n坐标矩阵。经验weak泛函按排序值乘经验尾质量求sup；严格阈值在有限经验分布中用从下逼近排序值解释。结果仅是MC，不是向外概率证书、有限维迁移或一般上界。

## 2. 实际结果：没有观察到weak反例，但共同winner不能白拿

共24项配置。观察到

- `W_MC(mixture max)/prior平均W_MC(product max)` 最大约.9930435；
- `/best supported product` 最大约.8530246。

本轮没有找到反例；这些比值小于1不能证明总体弱泛函比较。

n128、正幅度均匀prior的代表数据：

| m | mixture weak_MC | prior平均product weak_MC | best supported weak_MC |
|---:|---:|---:|---:|
| 1 | .4997803 | .5242111 | .7479724 |
| 3 | .5643139 | .5842094 | .7921373 |

与此同时，在正幅度范围内改变α时，**同一查询**的winning scale发生改变的比例为7.8%–32.1%。m1、n128时约32.05%。所以不能无条件把尺度max移入prior积分，并指定一个同时适合所有α的winner。

精确二阶量的解析公式为

\[
q_L=n\int d_L^2=\frac{n(.9)^2}{2}
\left(\frac{\sin(\pi mL)}{\pi mL}\right)^2.
\]

m1、n128时三个值数值约.497684、.964743、1.452780；m3时约.354543、.428775、.293011。固定cos来源下这种尺度间差随n线性放大，**不**满足接受源族的 `q_L=v+O(1)` 渐近平坦条件。本实验展示的是共同排序的缺失，不是证明方差平坦是任何可能weak比较的必要条件；它更不是接受源的反例。

## 3. log-concavity给出的严格信息与它不给的信息

对任意实数据d_i、`|d_i|≤h`、`0≤α≤A`、Ah<1，

\[
\partial_\alpha^2\log T_\alpha
=-\sum_i\frac{d_i^2}{(1+\alpha d_i)^2}\le0.
\]

因此每个查询的likelihood关于幅度log-concave；其最大化幅度可以依赖查询，不能由此选择同一个全局幅度代替prior。

写 `S_L=Σd_L(x_i)`、`V_L=Σd_L(x_i)^2`，幂级数还给全体α,L一致的精确余项包络

\[
\left|\log T_{\alpha,L}-\left(\alpha S_L-\frac{\alpha^2V_L}{2}\right)\right|
\le\frac{nA^3h^3}{3(1-Ah)}.
\]

若另有 `|V_L-v|≤D`，才可以在该好集上把不同L的quadratic penalty换成共同α²v/2，另支付A²D/2。对于非负α，共同penalty使所有理想likelihood按同一 `S*=max_L S_L` 排序。这正是root平坦方差证明使用的附加结构；不是仅由α的log-concavity自动得到的。

Jensen对log只给lower方向，不给 `∫T_αdν≤T_{Eα}`。weak-L¹也不能把函数相加直接作常数1三角式。任何上界必须明确用共同排序、阈值split或坏集输出质量付款。

## 4. exact n2例：通常不存在单产品的逐点支配

令单坐标base在两个原子±1上均匀，`d=δ sign`、δ=1/10，n=2。取非退化prior在u=1/5和v=1上各1/2。所有产品非负且质量1。

记平均幅度a=(u+v)/2=3/5。对++原子，

\[
\tfrac12T_u+\tfrac12T_v-T_a
=\delta^2\operatorname{Var}(\alpha)=1/625>0.
\]

实际四原子的Fraction值在JSON：混合在++/−−处比均值幅度产品多1/625，在两个异号原子处少1/625。

更一般，假设某个α_*满足此混合≤T_{α_*}逐点。两者质量均1，必须四原子处全部相等。比较++与−−的差强迫α_*=Eα；再比较++就与上式矛盾。因此**没有任何单个α产品逐点支配这个非退化混合**。

这个例子不是weak范数比较的反例，也不排除带1+o(1)或固定常数的比较。δ→0时差为O(δ²)，所以尤其不能把它解释为渐近log增益机制。

## 5. 范围和复现

本目录已完成真实有限核的同查询幅度prior探路、exact逐点反例及likelihood余项接口。关于接受源的“不损失无界log因子”结论，需要root同轮条件定理中的方差平坦、好集的anchor输出质量预算及有限尺度族条件；本目录数值没有认证这些条件。

未证明单产品weak范数上界，未处理一般conditioned input，未给一般C_n上界或新下界。

```sh
/opt/homebrew/opt/python@3.15/bin/python3.15 probe.py
```

仅依赖已有NumPy和标准库；脚本/JSON/本NOTE只写本目录。
