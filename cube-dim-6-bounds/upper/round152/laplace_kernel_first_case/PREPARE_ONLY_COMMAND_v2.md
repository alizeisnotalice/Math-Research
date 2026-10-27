# Terra唯一prepare-only命令（Sol未执行）

此版本保留原五个文件，追加明确root与group-reduced贡献口径。
原prepare_for_terra.py保留供来源审计，当前准备入口仅用v2。

下述root是147已存Terra模式资料所在的独立brascamp108 clone内
新152目录；命令会按白名单复制source packet并冻结terminal源。
目录中若已有不同源、准备输出或真实任务，程序停止而不覆盖。
该命令不会导入Arb/SymPy、执行single_case、density生成或科学
预检，也不会创建parent_release/started/certificate。

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round152/laplace_kernel_first_case/prepare_for_terra_v2.py --root /Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round152/laplace_kernel_first_case --primary-root /Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds --prototype-root /Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward --audit-root /Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-pushforward147-audit/cube-dim-6-bounds/audit
```

只有主与独审确认新冻结hash后才另登记科学launch。首例数值仍
由保留原single_case.py执行，新的frozen_sources会绑定追加口径：
它仅是group-reduced21fan中cone0单pair正贡献，不能解释为原始
固定pair在全R6上的事件积分，也不是完整beta2上端。

Sol此轮只进行文档/源文件的行政hash封存；未执行本命令、prepare、
科学预检或任何数值。
