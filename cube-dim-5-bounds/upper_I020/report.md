# I020：四子集完成项分组与分数 Hölder 的新完整 I 上界

本批取得新的完整辅助积分上界：

\[
 I\le 1.42147007315294853246451860546810\ldots\times10^{-5}
 <\boxed{\frac{2221047}{156250000000}}.
\]

该显式包围严格改善第006轮的 (1.57928752881910\ldots\times10^{-5})，约为旧包围的0.9000704731倍；仍未达到 (11/10^6)。这是完整Green、全部时间、整个正unit face的推导，不消费截时/空间采样预算。没有给新的 (g) 或 (C_5) 算子上界，(C_5\le32) 保持。新机制及其有限证书交主独审，不以本报告替代独立验收。

本目录由本代理独占；主clone和其他作者文件只读，没有提交或改写他人工作。输入及SHA见 `registration.json`。已读主round007排除：不使用 (alpha d^{-3/14}+A+Bd+Cd^2) 标量majorant族，也不把其抽象双点障碍当成真实形状分布。本批重新利用具体矩阵的Cauchy–Binet完成项及独立stable变量。

## 1. 固定接口与尺度抵消

仍是16条未单位化方向 (v=(1,\pm1,\pm1,\pm1,\pm1))，独立标准 (alpha=3/2) stable 时间，Lebesgue测度。记 (eta=3/4)，(S_i) 为独立正stable，(mathbb Ee^{-\lambda S_i}=e^{-\lambda^\beta})，

\[
 A=\sum_iS_iv_iv_i^T,\quad T=\sum_iS_i,\quad M=A/T,\quad d=\det M.
\]

形状概率 (R) 是以 (T^{-\beta}/Z) 倾斜再推送，(Z=1/[16\Gamma(7/4)])。已审019的轨道平均/Minkowski给

\[
 I\le c^{10/7}\frac{\pi^2}{60}J^{10/7},\quad
 c=\frac3{2^{15/2}\pi^{5/2}},\quad J=\mathbb E_Rd^{-b},\quad b=\frac3{20}.
\]

因为 (5b=\beta)，有**精确**的尺度抵消

\[
 J=Z^{-1}\mathbb E(\det A)^{-b}.
\]

不假定 (T) 与 (P=S/T) 独立。本批只改进右侧负矩，不改面角因子或一般输入传递。

## 2. 按共同四子集精确分组

令 (W=\det(\sum_i v_iv_i^T)=16^5)。对每个四子集 (J\subset[16])，取其5维整数cofactor normal (n_J)，使 (n_J\cdot v_j=0) 对 (j\in J) 成立。若rank小于4则 (n_J=0)，不计入以下正权组。

\[
 w_{J,i}=(n_J\cdot v_i)^2=\det(V_{J\cup\{i\}})^2,\quad
 L_J=\sum_iw_{J,i}S_i,
\]
\[
 C_J=\sum_iw_{J,i}=16\|n_J\|^2,\qquad
 \theta_J=\frac{C_J}{5W}.
\]

每个五方向基在其五个四子集中恰出现五次，所以

\[
 \det A=\frac15\sum_{|J|=4}\left(\prod_{j\in J}S_j\right)L_J,
 \qquad\sum_JC_J=5W,\quad\sum_J\theta_J=1.
\]

定义 (F_J=(\prod_{j\in J}S_j)(L_J/C_J)>0)。加权AM–GM给

\[
 \det A\ge W\prod_JF_J^{\theta_J},\quad
 (\det A)^{-b}\le W^{-b}\prod_JG_J^{\theta_J},\quad G_J=F_J^{-b}.
\]

在 (S_i=1) 时，所有 (F_J=1) 且左端 (det A=W)，此步取等；分组没有任意丢弃或重新标定CB质量。

关键是 (L_J) 不含选中的四个 (S_j)，所以它与这些四变量独立。令

\[
 p_{J,i}=w_{J,i}/C_J,\quad
 R_J=\sum_{i:w_{J,i}>0}p_{J,i}^{\beta},\quad
 m(s)=\mathbb ES^{-s}=\frac{\Gamma(1+s/\beta)}{\Gamma(1+s)}.
\]

(L_J/C_J) 的Laplace变换恰为 (e^{-R_J\lambda^\beta})，故对任意 (s>0)，

\[
 \mathbb E F_J^{-s}=m(s)^5R_J^{-s/\beta}<\infty.
\]

普通加权Hölder直接得到一个新完整界

\[
 J\le 2\Gamma(7/4)m(b)^5\prod_JR_J^{-\theta_J/5}
 =1.52555406301960771424\ldots,
\]

相应 (I\le1.44189314937231422369\ldots\times10^{-5})。这里确实让每组五个负矩抬到 (b)，但保留线性stable和的精确尺度后，全式仍严格优于旧证书；没有只比较其中一个有利因子。

## 3. 真实变量scope与分数 Hölder

组函数的真实依赖坐标是

\[
 A_J=J\cup\{i:w_{J,i}>0\}.
\]

不在该集合的 (S_i) 完全不进入 (G_J)。方向群在16个方向上传递，并同时置换 (J,C_J,A_J)，因此每个坐标的加权占用相同。本批精确有限几何给

\[
 \sum_{J:i\in A_J}\theta_J=\rho=\frac{53}{64}\quad(i=1,\ldots,16).
\]

取 (h_J=\theta_J/\rho)，则每个坐标恰有 (sum_{J:i\in A_J}h_J=1)。独立产品概率上，逐坐标使用广义Hölder可证

\[
 \mathbb E\prod_Jf_J^{h_J}\le\prod_J(\mathbb Ef_J)^{h_J}
\]

当 (f_J) 只依赖 (A_J) 时成立：对一个坐标，所有涉及它的指数和不超过1，先应用Hölder并积分掉该坐标；其他函数提出积分。对剩余坐标重复，指数和条件仍成立，最后成为各函数的完整期望。可先将非负函数截断，再单调取极限；不要求不同组函数独立。

令 (f_J=G_J^{\rho})，于是

\[
 \mathbb E\prod_JG_J^{\theta_J}
 \le\prod_J(\mathbb EG_J^{\rho})^{\theta_J/\rho}
 =m(b\rho)^{5/\rho}\prod_JR_J^{-\theta_J/5}.
\]

所有组负矩有限，所以该步骤同时证明完整期望的有限性，包括全部退化形状邻域。它不是对 (R) 强行使用独立性；独立性只在尺度抵消后原始 (S) 产品概率中使用。

最终新矩界为

\[
 \boxed{J\le J_{020}=2\Gamma(7/4)
 \left[\frac{\Gamma(373/320)}{\Gamma(1439/1280)}\right]^{320/53}D},
 \qquad D=\prod_JR_J^{-\theta_J/5}.
\]

其中 (b\rho=159/1280)、(b\rho/\beta=53/320)。这与普通Hölder的提升指数不同；不把负矩指数 (b\rho) 错留为 (b)。

## 4. 有限精确几何及闭式常数

预登记批次仅检查全部 (inom{16}{4}=1820) 个四子集的4阶cofactor：100个rank不足4，1720个rank4。没有再次枚举4368个五子集。非零primitive normal只有以下三个完整符号/置换类型：

| normal绝对值排序 | 四子集数 | 总 (	heta) | 非零完成方向数 | scope大小 | 非零 (p_{J,i}) 的多重集 |
|---|---:|---:|---:|---:|---|
| (1,1,0,0,0) | 1160 | 1/2 | 8 | 12 | 8个1/8 |
| (1,1,1,1,0) | 480 | 3/8 | 10 | 14 | 8个1/16、2个1/4 |
| (2,1,1,1,1) | 80 | 1/8 | 12 | 16 | 7个1/32、4个1/8、1个9/32 |

每个 (C_J)、占用坐标及分组概率都用整数/Fraction核验；总 (C_J=5,242,880=5\cdot16^5)，16个坐标的占用都精确等于53/64。另对一个明确四子集，用独立5阶Fraction行消元核对全部16个完成行列式平方与cofactor点积平方。程序保存类型代表、精确概率/权重及全部几何记录流的SHA，不归档全量重复枚举。

三个 (R) 的闭式分别为

\[
 R_1=8^{1/4},\quad R_2=1+1/\sqrt2,\quad
 R_3=\frac{7+8\sqrt2+3\sqrt3}{2^{15/4}},
\]
\[
 D=R_1^{-1/10}R_2^{-3/40}R_3^{-1/40}
 =0.89938058742941192127\ldots.
\]

因此最终公式只有固定Gamma和代数常数，不含空间/概率求积或未知的形状分布量。

## 5. 实际结果、证据边界与复现

本机Python实际完成一次固定384-bit Arb批，约0.070秒；退出码0。完整结果球为

\[
 J_{020}=[1.51039608935354144429209079702134053211\ldots\pm1.82\times10^{-85}],
\]
\[
 I_{020}=[1.42147007315294853246451860546810367934\ldots\times10^{-5}
 \pm3.45\times10^{-91}].
\]

实际比较确认 (I_{020}) 严格低于已锁旧球的**下端**，不是比较显示小数。向上有理上界 (2221047/156250000000=0.0000142147008) 经Arb与其上端严格比较后取得。完整球、runtime和源hash见 `results.json`。

新界比目标仍大，不声称达到 (I\le11/10^6)，也不能由这个上包围倒推真实 (I) 大于目标。它可以直接输入同模型一般 (L^1)/全尺度传递；但没有本批新 (g)，故没有 (C_5<32)。若进一步推进，四子集完成项只剩一个线性和可精确积分；再减少Hölder损失需要保留两个或更多完成方向的联合正结构，不能退回已排除的标量majorant族。

复现：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_I020/grouped_negative_moment.py
```

依赖Python标准库与python-flint 0.9.0；没有共享MATLAB、后台作业或图窗。必要白名单仅报告、登记、脚本、结果、manifest五项；不复制输入原件，不改变主clone。此批现锁定，后续不同k/两层分组另登记，不覆盖020。
