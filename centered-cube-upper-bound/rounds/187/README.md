# 第187轮：整个软参数尾的全局核变差预算

一般 n^(1/2+o(1)) 一致弱型上界尚未闭合。本轮独审通过：真实soft/scale两参数核的L1上确界预算，完整保留内部混合导数与移动硬面；先接原完整删除行、再合并唯一j，支付整个原互反正J-share软参数尾，无需good重叠下界。

取h_n=min(1/2,c sqrt(log(n+1)/n))，整个s≥1−h_n片的来源费用为O_c(sqrt(n)log(n+1))W_A，全部常数及原外系数见归档。结合第186轮，主缺口集中到较硬参数区的极小good重叠交通及其他原分支。辅助score演化保留一般阈值BV和删除费用，但阈值forcing的较小来源价仍未证，不能算第二个付款。

本机检查：6组40/60位矩积分、30组product上界公式、15组selector公式（8有效范围）及3组有理事件模型全部通过；归档隔离重放两个完整JSON逐字节一致。有限计算非区间证书，也不认证原森林/FIRST/history或一般上界。

[research.zip](research.zip)：53,981字节，15成员；SHA256 `15130a0090ea988e25fd17a9a9fcce36cdec90bac109ad5ad9a6ac82c21dc7fa`。含完整证明、详细LaTeX文本、独审和必要数值材料；不编译PDF。工作副本完整备份仓库外并逐字节核对。

解压后 `python3 reproduce.py --verify-only`；在已有mpmath 1.3.0环境运行 `python3 reproduce.py --output-dir <不存在的新目录>` 可重放。实际本机PYTHONPATH和证据范围见包内numerical_review.md。

本轮仅新增README与紧凑归档，并更新项目最新索引。1 Astra high、2 Sol high、1 Sol medium按独立任务交付后停止；无可用扣额明细。
