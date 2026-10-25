# 第005轮：有序正 face 的单形顶点积分证书

令 \(\Delta=\{0\le y_5\le\cdots\le y_1\le1\}\)，并令 \(u_i=(1^i,0^{5-i})\)（\(0\le i\le5\)）。body32 对五个自由坐标的置换不变：置换仅重标同尺度 Cauchy 参数。因此，除有限个坐标相等的超平面（五维 Lebesgue 零测）外，正 face 是120个置换 \(\sigma\Delta\) 的不交并，故对 \(q=6/5\)
\[
 \int_{[0,1]^5}H_\epsilon(1,y)^q\,dy
 =120\int_\Delta H_\epsilon(1,y)^q\,dy. \tag{1}
\]

第004轮的严格空间接口是 \(F(y)=H_\epsilon(1,y)\exp(\lVert y\rVert_2^2/(32\epsilon^2))\) log-convex。令 \(S=\operatorname{conv}(v_0,\ldots,v_5)\) 为任意五维单形，且 \(H_\epsilon(1,v_i)\le U_i\) 是严格顶点上界。写 \(y=\sum_i\lambda_i v_i\)，则 log-convexity 再接 AM--GM 给
\[
 H_\epsilon(1,y)^q\le
 \exp\left({q\over32\epsilon^2}\left[\sum_i\lambda_i\lVert v_i\rVert_2^2-\lVert y\rVert_2^2\right]\right)
 \sum_i\lambda_i U_i^q. \tag{2}
\]
定义准确单形方差
\[
 D_S=\max_{\lambda_i\ge0,\ \sum\lambda_i=1}
 \left[\sum_i\lambda_i\lVert v_i\rVert_2^2-\left\lVert\sum_i\lambda_i v_i\right\rVert_2^2\right].
\]
重心坐标的积分均为 \(|S|/6\)，所以
\[
 \int_S H_\epsilon(1,y)^q\,dy
 \le {|S|\over6}\exp(qD_S/(32\epsilon^2))\sum_{i=0}^5U_i^q. \tag{3}
\]
这是真正的连续误差界：\(D_S\) 对有理单形可通过枚举全部非空 KKT active supports 得到准确有理值；它不是 sparse 插值。

初始 \(\Delta\) 上，\(y_j=\sum_{i\ge j}\lambda_i\)，故方差正好是 \(\sum_{j=1}^5y_j(1-y_j)\le5/4\)，并在 \(\lambda_0=\lambda_5=1/2\) 达到。单个粗单形在小 \(\epsilon\) 下仍有不可用的大指数损失，所以需要有明确覆盖账的细分。

实际脚本 simplex_cover_audit.py 从 \(\Delta\) 出发、用确定性最长边中点二分直到120个 leaves。每一步恰为两个内点不交单形且有准确体积相加；因此 leaves 精确覆盖 \(\Delta\)。120个置换产生14,400个物理单形，但置换对称意味着只需在基本域保存120个 leaf 的顶点区间并缓存共享节点。脚本以有理 KKT 枚举每 leaf 的63个非空 vertex supports。实际运行（4.4秒、零 kernel evaluation）验证：fundamental volume \(1/120\)，总正 face volume 1，初始 \(D_\Delta=5/4\)，120 leaves 的 \(D_S\) 范围为 \(25/192\) 至 \(65/256\)。下一阶段可在至多120个基本域 leaf 上评估顶点区间，以 (3) 加总并乘120；不得因数值小而删除未覆盖 leaf。

下界不能只使用单形的六个顶点。log-半凸仅给内部上包络；第004轮 ODE 比较沿坐标线需要每条截面的两个端点，而一般单形截面端点落在 facet 内部，并非六个顶点。因此 \(g\) 仍须用包住单形的 axis-parallel box（32个顶点）施行 ODE，或使用较粗的一阶平移下界。单形化直接节约的是 \(I\) 的顶点上界，不能冒充为同强度的 \(g\) 认证。

本文件尚未输入任何 kernel 顶点区间，因而未给出 face 全域积分、cube 全域最小值或新的 \(C_6\) 界。
