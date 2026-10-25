# upper_detgap023：一个条件点的 logdet Jensen 正修正试点

只取最大方向 v=(1,1,1,1,1)、unit-time 最大值 z=t=1、x=(1,1,1,1,0)
一个条件案例。不重跑旧 order 的 retained-value 网格，不修改旧锁。

理论先核两个完整上界：全谱三次 log 多项式和 A≤Amax 下的 Hessian/Kronecker
正积分。有限试点采用384-bit Arb，重用 audit_order022 已核的 m(1),a(1),F(1)。
仅用CDF正Stieltjes箱包围M2/M3，不推高阶截断矩级数：区间[1/4,1]128等分，
[0,1/4]以0与(1/4)^jF(1/4)完整包围。CDF直接reciprocal-Gamma系数求和，
每残类64项及严格余项，不调用作者旧引擎。

Hessian正积分只在同一条件点用64个u箱及正权重∫(1−u)du包围。
输出两gap/安全乘子与局部条件kernel；不得把这个条件点乘子用于整个旧箱，
不得声称新全域g、G或C5界。一次批次，120秒checkpoint上限，不追加参数搜索。
