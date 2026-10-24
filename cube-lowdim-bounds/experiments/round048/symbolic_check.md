# 已实际执行的三项符号辅助记录

本文件从本轮已有工具记录整理，**没有重跑**。仅核三项解析恒等式/极限，不读取或认证58原子完整构造，不计算矩形面积，不形成新下界证书。

## 原执行命令

运行时使用本机数学venv：

```bash
/Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python - <<'PY'
from pathlib import Path
from hashlib import sha256
import sympy as s
p=Path('<用户提供的 cube2_lower_construction_20261024(1).md 的绝对路径>')
a,rho,e,eta=s.symbols('a rho e eta', positive=True)
N,k,U,M=s.symbols('N k U M', positive=True)
res=s.expand(4*a*a-4*(rho+e)**2-4*(a-rho-e)*(a+rho+e))
limit=s.limit((N-2*k)**2*U/(N*N*M),N,s.oo)
gap=s.expand((2-2*eta)-(1+2*eta))
assert res==0 and limit==U/M and gap==1-4*eta
print({'document_sha256':sha256(p.read_bytes()).hexdigest(),'candidate_sibling_exists':(p.parent/'candidate_final.json').exists(),'sympy_version':s.__version__,'radius_factor_residual':str(res),'finite_window_limit':str(limit),'other_cluster_exclusion_gap':str(gap)})
PY
```

上面仅将原命令中用户微信临时目录的绝对路径作占位脱敏；数学语句、解释器、元数据校验及输出字段均为原执行内容。原路径及说明文件来源由本轮主侧保留，文档SHA绑定如下；本摘录不冒称脱敏后的占位命令已再次执行。

## 原实际输出

执行正常退出0，命令耗时约0.04065秒；SymPy版本1.14.0。原输出为：

```text
{'document_sha256': '1ed2840d8673da3596d998e1e190c040cf8f6e3225e6bb7b809fbde2b7e83a42', 'candidate_sibling_exists': False, 'sympy_version': '1.14.0', 'radius_factor_residual': '0', 'finite_window_limit': 'U/M', 'other_cluster_exclusion_gap': '1 - 4*eta'}
```

## 证据边界

1. `radius_factor_residual=0`仅核4a²−4(ρ+ε)²的差因式分解。严格半径余量a−ρ≥1/D来自isqrt/floor的整数论解析论证，不是本次SymPy独立认证。
2. `finite_window_limit=U/M`核固定边界宽度k时的截断比例极限；实际见证原子偏移是否有相应统一有限k仍需完整构造/源码证明，本次缺件。
3. `other_cluster_exclusion_gap=1−4η`核光滑背景顶点团排除距离差；η∈(0,1/4)、包含关系、质量及严格超水平集极限由解析推导承担。

`candidate_sibling_exists=False`仅表示执行时说明文档同目录没有candidate_final.json，不证明该文件在整台电脑上不存在。文档SHA只绑定说明文本，不替代候选证书SHA或准确核验结果。以上三项均不能认证文档宣称的1.7031197176，也没有搜索、LP、核/生成元点评价、面积重算或新文献检索。
