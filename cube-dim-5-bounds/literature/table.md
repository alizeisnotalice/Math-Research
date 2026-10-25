# 五维公开文献与门槛

检索日期：2026-10-25。对象固定为 Lebesgue 测度上的中心轴平行方体、共同半径、一般 \(L^1\) 输入与全部正尺度。表中“公开最佳”仅指本轮已实际检索并核验的来源，不是全球穷尽性断言。

|来源|核验范围|五维可用结论|证据状态|
|---|---|---|---|
|J. M. Aldaz, *A remark on the centered n-dimensional Hardy--Littlewood maximal function*, CMJ 50 (2000), 103--112, Proposition 1.4, p.110|正是上述 cube/Lebesgue/弱 \((1,1)\) 对象；原 PDF 的命题与证明页已用 pypdf 读取|\(C_5\ge A_5\)，其中 \(1.7481246343484800401516719055882944<A_5<1.7481246343484800401516719055882945\)；完整有理端点见 `aldaz_a5_certificate.json`|原文已读相关页；本地准确根区间重算通过。PDF SHA-256 `05f80011467c1a734142e0d244e149a6abec75454b58b90e641bc75fb7d13286`|
|J. M. Aldaz, *Kissing numbers and the centered maximal operator*, arXiv:1811.10478v2 (2020), Theorem 1.2 and Corollary 3.4|其结论针对所有局部有限 Borel 背景测度，因此包含 Lebesgue；\(\ell_\infty\) balls 即轴平行 cubes|\(C_5\le2^5=32\)。该值对更大背景类的 uniform constant 是 sharp，但并不说明 Lebesgue 的 \(C_5=32\)|原文定理与 Cor.3.4 已读；HTML SHA-256 `4892ca4b28c99f3e3965ee9ea34ec8ebf81047a93906643601eabe7f9df92d61`|
|J. M. Aldaz, *The weak type (1,1) bounds for the maximal function associated to cubes grow to infinity with the dimension*, Ann. Math. 173 (2011), 1013--1023|同一 Lebesgue centered-cube 定义；p.1013 定理及维数单调说明已读|\(C_5\ge C_4\)，可把已有内部 R4 证书传递到五维；其渐近结论不提供新的显式五维数值界|PDF SHA-256 `5539d5f3421aa082df37ae91bf9e88b660bb15c356fb0f18db720c60ac08795d`|
|Y. Lin, *Two partial balayage principles for weak-type estimates*（作者网页，2026-10-25 访问）|网页列出的 cube 上界仅为二维 axis-parallel squares；球的高维公式不是 cube 结果|没有可直接代入五维 cube 的数值上下界；可作方法背景，不能用于五维门槛|实际阅读表格与 cube 小节；HTML SHA-256 `5b99ba8a8b6780de4b9fc15425e58b01c8419e1941c311ce4afe412caa4373d1`|

## 当前门槛

- 本轮已核验的公开数值比较为 \(A_5\le C_5\le32\)。
- 内部继承下界应使用 R4 的**完整准确区间**及其已审维数传递证明；不能用本表的 Aldaz 小数代替 R4 账本，也不能用公开式把 R4 重新命名为公开发表结果。
- 本轮检索还查阅了 2026 年 Lin/Palomar 记录、arXiv 和关键词检索；它们显示新的 cube 结论仅在二维，未找到可核验的五维专用改善。该负结果不是“不存在其他结果”的结论。
