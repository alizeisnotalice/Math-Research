# 只读源码盘点与真正新增范围

upper_theory_cont 本轮只读盘点，没有重跑生产或修改任何原件。完整α1 Fourier密度早已构造：body_fan_geometry.py的7室×3840=26880锥经004独审覆盖，body_fan_arb_anchor.py完整计算pt，不只是quotient。旧归一ray1762中280是真顶点、1482人工；005 star消去有准确证明。真顶点三轨道11000/16、11110/24、21111/36已准确分类，不能重启当新研究。

本次未在已核档案找到仅280真顶点的实际facet/边界simplex剖分；005_nonvertex_star_cancellation.md的pulling仅存在性。故新astra_vertex017只实现该280-vertex内核及准确几何覆盖，不重复完整密度存在性、原26880 fan或旧Fejer通过基线。

已有有限正时间下积分007：点(1,936/1009,673/1009,68/1009,917/1009)，eps1/128,T1，175 Fejer节点、degree24，H下端约.00014615177，仅单点。旧H1/24与H1/32全域.00011目标有严格反例，不能复跑同起点期望通过。1/64、1/128八点仅有限点过，没全域面覆盖。

原009 whole-box逐锥相消诊断只做interior16盒×225时间节点，79.34秒，盒包围远大于中心Arb值；near-S2/boundary/corner没运行。新I-alpha1 shape闭合不自动消除这种全盒相消，故新017第一阶段只几何/内核及两个新资格anchor，不开global mesh。旧303595776节点均匀空间计划从未启动，不当已失败网格或存活句柄；磁盘009二进制存在不意味着进程活着。

原件hash（真实只读盘点提供，主未重跑）：

- numerics/body_fan_geometry.py: 84dad4a7c74297fcb4e63377276d4a56e7dfa75c038063c2e3f319954aa1ce43
- numerics/body_fan_arb_anchor.py: 348f0d348e084c4fb5ce62ce7194a9a6dfa58b5e1244e6a883a0fe7ae3f1f9c3
- audit/005_nonvertex_star_cancellation.md: c611c15c7bce10b5691b52a371b46d770bc73ffcb9f1f7d880f6adfb06e8f01c
- numerics/007_final_verify.json: b8d07d5b8353c5da5b6979e1457ae4eddd428c2cb9c2329e2ab4c6ca8288b3a3
- round003/verify_cutoff_counterexamples.json: fee97fc3911700e2d4e5bde84c424649c750d79065e08ebbdcd681a7df2dae03
- numerics/009_diagnostic.json: de23ffcd360d3caf7ade03bb8d58138e24b34d6015cccd0924c71e239d4c6699

旧浮点body_fan_check.py行列错误与旧cosine漏2^-5判断均已撤回，不允许重用为有效来源。source005_vertex_rays_exact.py位于独立audit clone，保持只读。新内核Cauchy仍须用同ψ、(2π)^-5归一化，两anchor相合不能冒充全域g。
