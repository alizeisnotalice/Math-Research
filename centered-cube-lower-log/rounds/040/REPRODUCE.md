# 复现

Python3和NumPy即可。默认入口：

```sh
python3 verify_saved.py --output /tmp/cube-round040-new-audit.json
```

它在临时副本中重算504项CI、ESS、source-only帽、RN、成对统计及84个same-h比较，重放一维source验证、相关性和混合倾斜诊断。原件只读，不抽高维新样本、不建立大缓存。结果文件必须尚不存在。约需数秒至十余秒，依本机环境而变。

完整pilot约需0.6GB内存，原主运行243.4秒：

```sh
python3 replay_pilot.py --output /tmp/cube-round040-new-pilot
```

该输出目录必须不存在。脚本复制冻结source/builder/driver/helper到新目录，程序先写自己的PLAN再抽counts。R=4/8/16/32，全部233尺度，N256和seeds均不变；重放属于同协议复算，不是独立确认。

原主运行Python3.12.14/NumPy2.3.5；部分生产验证与独审Python3.15.0/NumPy2.5.4。不同NumPy的multinomial流不承诺逐位一致。对新运行比较样本/统计hash时，要同时保留版本及完整调用顺序；运行时间与RSS不要求一致。

`R*_paired.npz`保留逐样本模型贡献、logD、全component和的log、proposal标签，允许独立核验保存统计。raw counts及每个b/j的核统计没有归档，不能从小payload独立复原全部RN分量或全部原查询；需要完整协议重放。npz均以allow_pickle=False读取。

source验证和一维诊断需要38轮冻结source，默认核验器从仓库`../038/round038_highdim/source_tables.py`映射到临时旧路径，不重复归档。主pilot本身只依赖40轮四份冻结代码。

source的解析中点误差与经验FP guard分开；完整BLAS、CDF和随机提议未向外认证。固定阈值CI及有限尺度包络不等于完整weak排序。源码和证据修正只涉及一份独审元数据断言，实际pilot没有重跑或修改。
