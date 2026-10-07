# 高势来源一次截断：三轮有理组件守卫预登记

2026-10-07。登记先于执行。只新增本 prefix；不改主账、旧结果或 MC 数据。

## 解析对象与范围

原完整捕获分母 m、阈值 m>τR^n 和一份全局来源集合 B={V_hi>H} 给 fullness 子支的径向界
J_full≤W_B ℓ(I_full/(ηW_B))，其中 ℓ(r)=r (r≤1)，ℓ(r)=1+log r (r≥1)。在 T=e^(H/2) 处的凹函数切线为 ℓ(r)≤log T+r/T，因此 J_full≤(H/2)W_B+(ηT)^−1 I_full。fragment 尚未付。

本次只核这些标量构件和事先指定常数组合，不生成空间 input，不评估原 winner/FIRST/fullfuture/CPGP/LCA/history，不做 MC，不声称证明 sqrt(n) 阶数。

## 固定三轮

n=16,64,256；s=√n=4,8,16；T=2^(s+2)，H=2log T，η=1/4，δ₀=1/8，候选尚未证 fragment allowance δ₁=1/4，ν=1/4，K=s。log2 用 atanh(1/3) Fraction 级数上下界，项数分别 16,24,32；余项上界 2z^(2N+1)/((2N+1)(1−z²))。

每轮检查：

1. log2 区间嵌套与宽度、H>0 且 H<1+nlog2（不靠已知强包络截点变成平凡尾）。
2. 切线在 r=0,1/2,1,T/2,T,2T,T² 的方向；r=T 用代数恒等而非区间强求等号。
3. 正 full fee (ηT)^−1=4/T；δ₀+fee≤1/4；若 fragment allowance 真能支付，则 δ₀+fee+δ₁<1−ν−1/K。最后一项只核充分条件常数，不断言 fragment 合同成立。
4. b=c/m 的 low/full/fragment 分类与 fragment 的 m<ηW_B/δ₀，及 ηW_B/δ₀≤τa^n 时 fragment 不可能。仅用有理 scalar fixture，不冒充可实现的空间 gates。

只执行一次脚本，保存 JSON 结果与含脚本 SHA256 的收据；断言失败时保留失败，不自适应更换参数。没有新空间估计阶，所以不再追加模拟或旧数值重跑。
