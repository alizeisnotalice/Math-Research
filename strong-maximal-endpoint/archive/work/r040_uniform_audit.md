# r040：正中偏差 uniform-η prefactor 审计

日期：2026-11-06。目标区间为
\(\eta\in[A\sqrt{\log t/t},\delta]\)。结论：新的率函数
\(J_{new}=I_*-\varepsilon D\) 的 Hamiltonian 负余量是正确的；参考半群归一化也确实
消掉了裸的 \(1/t\) prefactor 项。但固定 \(\eta\) Fourier 证明不能直接均匀化：
其中关于 \(v/f\)、\(Q/f\)、call 的尾积分必须以 two-scale Mills/saddlepoint 展开重做。
在该补充完成前，不能登记 \(O(\eta/t+\eta^4)\) 的严格 uniform transport lemma。

## 1. 率函数与形式 prefactor 机制

\[
J_{new}=I_*-\varepsilon D,qquad D\sim A_3\eta^3,quad D_y=O(\eta^2).
\]
Taylor Hamiltonian 计算给
\[
H_{J_{new}}=-\varepsilon D+
 \tfrac12\varepsilon^2\Lambda_*''(D_y)^2+\cdots<0
\]
于充分小固定 \(\varepsilon,\delta\)，故 r039 的正四阶差 \(F\) 在率层已被正确移除。

若把 transport 比值相对于**精确** \(P_a\) 半群归一化，则在 \(\eta=0,c=a\)
恒有 \(P_av_{t-1}/v_t=1\)。所以单独 saddlepoint 写法中来自 \(t^{-1/2}\) 的
generic \(1/t\) 项不应保留：光滑 action 窗在 \(c-a=O(\eta)\) 下只给
\(O(\eta/t)\)。\(\beta(\eta)\log t\) 的附加项可产生
\(O(\varepsilon\eta^2\log t/t)\)，但在 \(t\eta^2\ge A^2\log t\) 时可并入
\(O(\eta^4)\)。这支持所追求的误差尺度，**条件是**下述 tail expansion 也保持该
消去。

## 2. 固定倾斜论证失效的位置

固定 \(\eta>0\) 时，可将 \(v/f\) 的尾积分截到 \(z=O(\log t)\)，剩余由
\(e^{-\eta z}\) 任意阶压小。下端
\(\eta=A\sqrt{\log t/t}\) 不满足此点：
\[
 e^{-\eta C\log t}=e^{-CA(\log t)^{3/2}/\sqrt t}\longrightarrow1,
\]
而 \(Q/f\) 的主质量在 \(z\asymp1/\eta\gg\log t\)。因此不能把
`bounded-z Fourier expansion + outer exponentially small` 用于 h0/h1 或 call。

最简单的正态模型已经显示奇异余项。令
\[
 Q_G/f_G=\int_0^\infty e^{-\eta u-u^2/(2\nu t)}du.
\]
写 \(x=t\eta^2\)，其 Mills 展开为
\[
 {Q_G\over f_G}={1\over\eta}\left(1-{1\over\nu x}
 +{3\over\nu^2x^2}-{15\over\nu^3x^3}+\cdots\right)
 ={1\over\eta}-{1\over\nu t\eta^3}
 +{3\over\nu^2t^2\eta^5}-\cdots.\tag{1}
\]
故写成 \(h_0+h_1/t+O(t^{-2})\) 时，\(h_1\asymp\eta^{-3}\)，且余项不是
uniform \(O(t^{-2})\)。WolframScript 15 实际将该积分化为 Erfc 并给出 (1)。

## 3. 必需的替换与 matched 参数

需要在 \(x=t\eta^2\ge A^2\log t\) 下建立 uniform tilted Mills expansion：保留
高斯 \(e^{-u^2/(2\nu t)}\) 因子，或展开到有限但足够高的 \(x^{-1}\) 阶，并对
\(Q/f,R/f,\) call 及 auxiliary \(C_t\) 同时给出余项。之后再在 exact-reference
normalized ratio 中逐项核对这些 Mills 修正；不能预先声称它们不影响 \(\eta/t\)
取消。它们可能以 \(tD\cdot(t\eta^2)^{-1}=O(\eta)\) 进入 prefactor，虽未必破坏
结果，却需计算。

对所有 \(t\le N\) 也不能以 \(t|a_N-\mu|=o(1)\) 替换 matched 参数；在 \(t\asymp N\)
此量很大。应全程以 \(a_N,\ell_N\) 写精确基准，并使用
\(|a_N-\mu|\lesssim\sqrt{\log N/N}\le C\eta\)（在上述下界和 \(t\le N\)）
处理动作邻域。有限 \(N\) call 项在固定 \(\eta\) 为 \(O(1/N)\)，但同样应在该
uniform bookkeeping 内保留。
