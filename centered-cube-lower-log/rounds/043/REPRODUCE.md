# 第43轮复现

依赖Python、NumPy、SciPy。倾斜数值原运行Python3.12.14/NumPy2.3.5；独审与高维pilot原运行Python3.15.0/NumPy2.5.4/SciPy1.18.1。使用已有本机环境，不安装新依赖。root的捆绑Python未装SciPy，本轮统一包核验实际使用已有的Python3.15环境。

```sh
python3 verify_saved.py --output /tmp/cube-round043-check-new.json
```

默认核验：全部payload SHA；在临时副本重算9000项倾斜摘要、exact来源、先验mismatch探针、谱价、独立Fraction/恒等式检查、高维pilot的质量及小DP检验、最终探索性CI汇总。不会重新训练或重放全部heldout。数值按容差比较；精确分数与整数逐字比较。原文件只读，输出目标必须不存在。

```sh
python3 verify_saved.py --replay-test --output /tmp/cube-round043-heldout-new.json
```

此模式额外重放全部六组冻结测试的七阈值guard计数，约25秒，无新独立样本、不重训。它仍只是冻结浮点分类器的复现，不会补出高维区间/连续query格证书。

完整训练复现时，复制 `round043_binary` 到一个全新输出目录，移除该副本中原SUMMARY/BALANCED_SUMMARY、PLAN和FROZEN/RESULT文件，依次运行 `probe.py`、`balanced_probe.py`、`replay_guards.py`、`combine.py`。两个计划和所有阈值/seed在源码中；balanced组是看到free组坍缩后追加的独立pilot。已有文件的原归档不可覆盖。不同NumPy版本若改变随机流或搜索并列，必须另列新结果，不能宣称复现原冻结选择器。

`round043_tilt/probe.py`只读引用第42轮 `round042_channel/probe.py`，先尝试同级本地研究目录，再使用仓库中 `rounds/042/...`。核验器把同SHA helper映射到临时视图；没有在第43轮重复归档旧脚本。只用它的确定性源/核函数，不执行第42轮主扫描。

六份binary独立RESULT与两份SUMMARY中的对应row完全重复，所以Git只保存SUMMARY中的原件。核验器在临时视图逐row恢复这六份JSON，并用 `shared_binary_result_hashes.json` 核对其字节SHA，再调用原verify_saved；这种映射不是新增模型执行。原本地RESULT、日志和旧MANIFEST留在工作输出目录，不上传Git。

9000项倾斜扫描和二类别pilot都是数值探索；严格理论来自自含证明及有理反例。高维pilot没有全流程向外FP或continuous-cell认证，其CP区间不可当作Lebesgue弱型证书。辅助trace的指数增长不能当作C_n增长。本轮没有新C_n阶数结论。
