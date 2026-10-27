# 八维中心方体极大算子论文

研究已按用户要求暂停。本目录为已完成结果的中文论文及便携复核材料，没有启动新参数实验。

## 正式结论

2188811179426669/10^15 < C8 < 56627260639/2000000000。

对象固定为Lebesgue、中心轴平行方体、所有坐标共同半径、一般L1及全部正尺度。下界采用暂停前已经独立全标签认证的E098构造，完整源分数与准确封装数据见supplement/paper_claims.json。Green奇异层、概率源及比较原理单列；三个点值及更紧长时界不被当作g/I证书或新Cn上界。论文不声称Cn或整个见证族最优，不把内部成果称为公开发表记录。

## 文件

- main.tex、sections/、references.tex：可编辑论文源码。
- main.pdf：实际XeLaTeX编译的论文。
- supplement/paper_claims.json：原始准确分数、直接由原分数向外舍入的紧凑界、来源hash。
- supplement/verify_paper.py：只依赖Python标准库，重新核本稿唯一完整下界并及上界有理对数尾。
- supplement/verify_green.py：可选的精确全锥和/积分核验，需要python-flint与SymPy。
- lower_audit_*：暂停前已完成的原独审证据副本。原完整审计器含原项目来源追溯，便携数学复核请使用verify_paper.py。

## 复现

在本目录运行：

```sh
python3 -B supplement/verify_paper.py
python3 -B supplement/verify_paper.py --scalar-only
python3 -B supplement/verify_green.py  # 可选，需python-flint、SymPy
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

完整标准库复核已实际PASS：32768标签、210446缓存状态、34054包含检查，用时约6.09秒。此封装复核复用原独立固定轴扫掠算法，不称第三套独立实现，不做搜索。Green便携复核已实际PASS：1152锥、三个精确点值及全空间峰。编译与数学核验分别记录。

主定理下界高于已核公开Aldaz A8及开局内部七维完整证书，上界低于已核通用256；参考范围不是全球最优性声明。文献采用已核原始资料与正式题录。作者姓名未由用户提供，源码author留空。

## 暂停边界

R006 tail第二参数工作未纳入数值主定理。研究目标状态paused，后续只按用户新指令恢复。论文准备中的源核对、包装回归和编译不是新研究轮次。
