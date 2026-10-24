# 第20轮：Tomar v4复现与第一轮缩放优化

公开源锁定bca281a706bbf6375b04d601ad19ba4ad0bf28b6。upstream-v4仅保留复现必要的原封不动核验器、公开证书、格式及Apache许可证；不是完整上游仓库，不包含Lean工程。SOURCE_MANIFEST.json保存来源与文件哈希。

要求Python3.12、python-flint==0.9.0。实际运行摘要见RUN_SUMMARY.json；四次完整核验均在本机完成，第一套两份完整短日志保留，第二套逐叶报告/全日志保留本地且存哈希。second_public.json是实际报告的最小C2摘录，非全报告。

从本目录进入upstream-v4后，可独立核验新证书：

```sh
python3 -I checker/cert_run.py ../round01/scaled_certificate.json --Q 8 --maxdepth 7 --Kdepth 6 --procs 6 --only MKA
python3 -I checker/independent/hl_check.py ../round01/scaled_certificate.json --claim 3.59091 --procs 6
```

将证书路径换为cert/cert_h32_a4740_R3320_3p59123.json、第二命令claim换成3.59123，即复现公开基准。第一程序PASS不自动比较成本阈值；其完整质量球必须与准确阈值比较。exact_cost.py从本目录的固定证书独立用整数40次根与Fraction隔离公开准确成本；不调用两核验器或Arb。新成本为其100000/100009倍，准确端点和严格增量见round01/scaling_arithmetic.json。

round01/make_scaled.py记录唯一候选生成规则；现存证书时故意拒绝覆盖。需重新生成时在独立副本中移走旧候选再执行，不能将此生成操作当作新搜索。两核验器共享Arb及数学公式，故只是实现层独立；整数成本隔离不共享Arb。没有本机Lean完整编译。

此次是三轮计划第1轮：一个预定有理缩放因子、两次新完整核验各600秒上限；没有失败后补参数。后两轮的具体预算将在运行前依据反馈确定。旧第20轮下界草稿因用户改目标而停止，未运行、不归档。
