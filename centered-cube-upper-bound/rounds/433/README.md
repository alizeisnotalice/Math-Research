# 第433轮：真实共同父图与实际尺度的静态双来源预算

一般一致 n^(1/2+o(1)) 上界仍未闭合。本轮在完整父—双来源层建立两项一般静态矩：

- 由正Lebesgue中心测度定义真实共同父图及来源度U；以实际父质量归一的指数距离矩≤[1+2n log(b/a)]W。
- 在实际R_A处使用局部来源质量，保留全部尺度恢复项；以网格Jensen将局部增长熵支付为n log5·W，得到第二项静态矩。

原sourcewise S推前后，两项阈值支顺序支付，总新增2S(1+√n)W，无新增X_q，不再含终端Θ²。完整低指标交集仍待一般预算；没有把带权矩当作无权全部闭合，也不退款旧费用。

[研究归档](research.zip)含详细LaTeX源、两项证明、独立审计、三份Python脚本和本机MATLAB输出。依赖原静态来源帽见[第413轮](../413/README.md)，当前原余流与旧账见[第432轮](../432/README.md)。

解压后运行 `python3 r433_geometry_check.py`、`python3 r433_joint_check.py`、`python3 r433_parent_check.py`。MATLAB运行 `r433_parameters.m`，结果写入 `/tmp/`。脚本中的连续模型、标量证书和浮点诊断均说明范围，未认证完整原历史。

仅新增本README和ZIP、更新项目总览；保留LaTeX源，不编译PDF。后续按用户新分工集中完整剩余交集控制，不以继续增加局部分支作为主要成功标准。

归档SHA256：`fdbc2402bedc2dd16c9836599887392b4dc7cc5bbcb118e98ecfcb2ce3e5d67d`。
