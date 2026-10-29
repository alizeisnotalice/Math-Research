# 复现与共享证据

需要Python和NumPy。默认入口：

```sh
python3 verify_saved.py --output /tmp/cube-round041-audit-new.json
```

它检查归档hash，在临时副本中重算variance/MGF小例、固定目标方法统一CI例、全部逐尺度上包、最终1008个CI及配对统计。不会重建高维source或抽取新尾样本。原件只读，输出文件须尚不存在。

四个约3.7MiB的`R*_upper_inputs.npz`只存完整离线复核upper必需的逐尺度log均值/variance及两个log-product预算，没有n维查询点、counts或R×J×Q特征。保留它们可以在NumPy随机流变更后仍验证同一批上包，不能用将来另抽样替代这些原观测。

实际230秒重放得到的paired NPZ与40文件字节完全相同；本轮保存四份真实重放RESULTS元数据，其中的paired hash直接指向40文件。因此默认核验将40的同hash字节映射到临时base_replay视图，不重复归档。这个映射只核对共享证据，不是再独立执行模型，也不重新建立原运行身份；原运行身份来自已保存receipt、实际重放metadata及当时两路审查。

完整同协议重建可另运行：

```sh
python3 replay_full.py --output /tmp/cube-round041-full-new
```

新目录必须不存在；脚本在那里重新生成counts/source统计，要求三种流hash及全部旧payload数组与40逐位相同，然后重新生成upper及最终分析。约需4分钟、0.6GB RAM。原主运行Python3.12.14/NumPy2.3.5；跨版本随机流不承诺一致，hash不符会拒绝接受。它仍是旧数据同协议重放，不是新独立确认。

方差旧稿、写作时尚未应用的uniform候选和独审原字节均保留；最终采用状态以`round041_post/FINAL_PLAN.json`及总报告为准。phase1重复RESULTS没有再归档，因为FINAL包含全部相同lower数据；保留phase1脚本和登记，原全量阶段记录仍在本地工作目录。

完整FP未认证。指数上移floor只处理一个下溢方向；高精度CI反演不修复source/BLAS/erfc或包络夹逼的数值证据缺口。有限h/J的结论不扩成全尺度weak或C_n新阶。

复算 FINAL_RESULTS 时先逐字段核对数值（运行时间与运行环境字段除外），随后在临时视图恢复已归档原件，供后续独审核对其冻结 SHA；不把含新运行时间的 JSON 宣称为逐字节相同。
