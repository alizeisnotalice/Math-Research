# 第4轮：证明与最小复现入口

[本轮报告](../rounds/004.md) / [精确区间](strict-interval.json)。当前上下界针对一般L¹、固定Lebesgue背景、全部正尺度。没有精确C₂或全局达到性结论。

- [核证明](proofs/kernel.md) / [独审](reviews/kernel.md) / [便携终审](reviews/kernel-portable.md)。
- [改权及校准jet证明](proofs/weights.md) / [独审](reviews/weights.md) / [便携终核](reviews/weights-portable.md) / [重复日志去重补签](reviews/weights-portable-cleanup.md)。
- [核包](kernel/README.md)：需python-flint 0.9.0；默认重核192位C1/C2/积分及历史新全C4覆盖，另给完整C4重跑。
- [有限输入](weights/README.md)：两步的实际有理输入、原双核日志和有限化标量。
- [有符号校准分数](jet/README.md)：生成的R3 records不重复归档；有理一阶并集、直接副本计数与有限步核。

从项目根，可复算最终有限下界：

```sh
python3 -I checks/lower/exact_check.py round004/weights/transfer_34_27.json --Dr 1099511627776
python3 checks/lower/lower_check.py round004/weights/transfer_34_27.json 1099511627776
```

第二面积实现需NumPy；它不依赖第一实现的输出。输入A的复算只需替换同目录candidate_05.json。周期输入的有限L¹传递必须读证明，不能只由两个程序的PASS推断。默认jet读取R3生成目录，建议按jet/README指定独立可写输出目录。stats.declipped是未使用占位，不代表无裁边；数学计算直接对整个端点jet裁切。

冻结正文和独审按原字节保存，其旧路径或“待审”措辞属于历史来源；当前裁决见配对最终审文。便携代码的变化限于路径、参数与输出位置，源映射和省略资料见source-bindings.json。原始用户提示、未验证内部材料、环境、缓存和大规模探索结果不公开。
