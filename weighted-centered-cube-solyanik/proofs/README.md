# 平移证书与历史证明索引

当前正文以 `paper/homogeneous_extension.tex` 的统一齐次权定理替代此前三段重复证明。
仍有独立量化价值的匹配范数平移界与方体精确二阶系数保留于 `translation_certificates.tex`；
在本目录以 XeLaTeX 编译 `certificates.tex` 可独立复核。

旧轮报告按当时文件名记录。相应已推送原件可从以下历史提交读取，不改写历史：

- [第15轮方体首项与一维二阶](https://github.com/alizeisnotalice/Math-Research/blob/3627dc8/weighted-centered-cube-solyanik/paper/power_endpoint.tex)
- [第16轮统一二阶与精确平移系数](https://github.com/alizeisnotalice/Math-Research/blob/376e5aa/weighted-centered-cube-solyanik/paper/power_second_order.tex)
- [第17轮任意范数球推广](https://github.com/alizeisnotalice/Math-Research/blob/b108f38/weighted-centered-cube-solyanik/paper/norm_geometry.tex)

删除当前版本中的重复专章不会删除相应结论：它们均由更强的一般定理涵盖。
历史本地第17轮交付及第18轮清理前备份也保留原有可编辑源文件。

## 第21轮清理的临界旧专章

正文第11轮的一般无原子永久见证证明已覆盖有限质量、一维及乘积权三个特例。以下旧专章长期未被主入口引用，现从当前paper目录移出，结论与原件仍保留：

- [第8轮有限总质量证明](https://github.com/alizeisnotalice/Math-Research/blob/b375524/weighted-centered-cube-solyanik/paper/finite_mass_critical.tex)
- [第9轮一维无限质量证明](https://github.com/alizeisnotalice/Math-Research/blob/13d9fe4/weighted-centered-cube-solyanik/paper/critical_infinite_extension.tex)
- [第10轮乘积权证明](https://github.com/alizeisnotalice/Math-Research/blob/d0053e5/weighted-centered-cube-solyanik/paper/product_critical_extension.tex)

本机第21轮交付包的support/historical_proofs/round008-010/保留可编辑副本；当前一般证明为paper/general_critical.tex。相关历史计算脚本保持算法，仅更新证明定位注释。


## 第24轮统一角向稳定证明

原第23轮paper/global_l1_cubic.tex的三次二次常数现作为更强全p包络的推论保留，统一源文件为paper/angular_envelope.tex；共同Jensen步骤只证明一次。旧源可从[第23轮提交](https://github.com/alizeisnotalice/Math-Research/blob/e3576e03106fed235e5503dd7641c97e273a983f/weighted-centered-cube-solyanik/paper/global_l1_cubic.tex)读取，且本机第23轮独立交付LaTeX中保留原件。文件重命名和证明合并未删有效结论，旧轮报告按当时路径及公式编号理解。


## 第29–30轮当前端点依赖

- `homogeneous_extension.tex` H5：集合浓度H5a直接给所需角向矩；测度扰动与G1给Δ_v；H6环带Tonelli、H7–8共同核心、H9截断体积归一给匹配上下界。二阶H11仍要求H1，未随首项扩大假设。
- `homogeneous_applications.tex` A6–7：一维幂权浓度、方体切片和锥矩给有符号角幂的严格次临界区间。
- 同文件AC1–5：临界权在中心方体的精确局部包络与全平移方体粗浓度，自含切片推导。
- AC6：只重复H6的短类求和（不使用H5a严格范围的结论），长类以全w浓度保证实际半径有界，再在共同核心内用AC3；收缩孔的严格中心见证给匹配下界。
- AC10依赖既有共同中心不等指数精确传递；AC11只用AC3的显式误差和平衡短/长类，仅为可证余项，未给出最优阶或精确第二项。

有限数值验证不进入上述逻辑依赖，完整弱型A1分类也不是AC6的必要前提。原件及历史公式编号按各轮当时状态保存。
