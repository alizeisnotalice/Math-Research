# cb9：三维体对角Cauchy锥/Green接口预登记

2026-10-25。只新增本目录cb9_*，不改旧证据/主账本、不扩招提交推送；不是唯一工作者。已读主第9轮登记、input_snapshot并逐文件实际核SHA；发现ledger与snapshot不匹配并已报告主代理，公式接口文件匹配。未用ledger作新旧比较，及body_diagonal_frame_interface。外部design仅待证猜想。源hash见cb9_sources.json。

固定V={(1,±1,±1)}，L=−Σ|v·D|，标准Fourier inverse (2pi)^−3。首批仅固定这个frame：准确生成6方形面各2三角与8三角面共20 simplicial cones，核det、全facet描述、radial coverage、体积归一；不做域形状/核成本优化。

先交cb9_frame_green.py可调用green(x)、green_face(y,z)、green_array(X)。验证p_t锥公式并保留全部复项，再推G的真表示。generic log-divided-difference只在允许的全项积分条件下使用；zero/repeated forms按总和或真实正积分处理，不分别丢发散项。尝试将p_t合为单个全实有理式，并以正Cauchy卷积/Feynman参数作独立正常数和退化校准。body-diagonal线奇性不得改成有限值。

有限交叉：最多12固定泛型有理点（程序中预列）与6固定退化类别：原点、非零body diagonal、坐标轴、coordinate plane、一个zero form如(1,1,1/3)、repeated nonzero cone forms如(1,2/3,1/3)。泛型点核p_t两表示和G解析/真实正积分；退化点分别核相应finite limit或infinity。每点最多两精度，禁止随机补点/无限近线扫描。锥集合是准确有限计算，解析覆盖证明不由采样替代。

Python主动实际，venv=/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python。MATLAB仅交cb9前缀便携脚本给主串行，不clear/path/reset。候选/接口交原独审；数值Green不等于全cube g/J或新U3/U4。停止于20锥准确覆盖+可调用数值接口+预列18点以内交叉+完整公式/源hash/缺口；更大精度或新机制另登记。
