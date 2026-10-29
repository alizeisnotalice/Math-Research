# 单项定理核验记录

日期2026-10-29；仅核验本轮脚本采用的经验Bernstein公式，不进行一般文献综述。

- 官方COLT2009原论文PDF：https://www.cs.mcgill.ca/~colt2009/papers/012.pdf （9页）；读PDF第1页Theorem1假设、第2页Theorem4及其样本方差定义，第4页Theorem11定位。未声称全文精读。
- 题录交叉核对：https://arxiv.org/abs/0907.3740 ，作者Andreas Maurer / Massimiliano Pontil，2009年7月21日，会议备注COLT2009。题名与PDF一致。
- Theorem4一侧失败概率δ，log(2/δ)，系数7/[3(N−1)]；对C−Z另付一次方向费用。本文不采用Corollary5中不同形式的log因子。
- 本轮底层抽样与浮点阈值尚未向外认证，故仅给理想精确模型的条件适用解释。
- 不下载重复全文，不引入未核验DOI或页码范围。完整BibTeX见references.bib。
