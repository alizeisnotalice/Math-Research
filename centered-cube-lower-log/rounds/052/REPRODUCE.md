# 第52轮复现

实际环境：本机Python3.15.0、NumPy2.5.4、SciPy1.18.1。其余只用标准库；本轮不需要mpmath/MATLAB。没有保存原始查询数组或临时缓存。

```sh
PYTHONDONTWRITEBYTECODE=1 python3 verify_saved.py --output /tmp/cube052-verification.json
```

验证器首先检查payload及第45/50/51轮证明引用哈希；在临时副本中真正重跑五个脚本，再比较全部五份JSON摘要。原归档不被改写。耗时字段跳过，其他整数/有理字符串精确比对，浮点以声明容差比对。相同MC seed实现计算复现，不代表数学上的概率误差认证。

各脚本可独立运行：flat_budget/probe.py输出12项精确有理预算（D有限值仅为假设）；check_order.py输出64项Fraction共同排序及无共同排序反例；amplitude/probe.py运行24项真实连续核prior探路；tailgate/probe.py读取FROZEN_PLAN后运行独立训练/测试计数门；verify_exact.py执行1824项单侧尾方差有理检验。

必要结果摘要保留有限阈值profile以供人工核查训练/测试选择及粗网格误差。成本pilot只保留在本机工作目录，不是复现依赖，未归档。最小FROZEN_PLAN的cost_pilot字段仅记录其历史发生时间关系，不要求读取该文件。

严格结论见flat_budget/PROOF.md和tailgate/PROOF.md。连续周期核实验和抽象有限原子反例分别标注；不得将后者误当真实立方体来源。
