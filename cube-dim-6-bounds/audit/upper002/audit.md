# 独立审计：第 002 轮六维 body fan 的精确 section 覆盖

审计日期：2026-10-25。对象是
`Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round002/fan/` 的
`exact_cells_results.json`，读取 SHA-256 为
`18d8aa578f84283d00e1c9c30522fb14cc49812092afd5716f2cc9ae68ba2798`。
生产器与其同作者 checker 本身不构成独立证书；本审计不导入生产器，也不从
certificate 读取 crossing forms、rays 或 cell 分类作为构造输入。

## 结论

**接受有限 fan 的 section 分解、cell/ray 绑定、体积和 peak 积分。**

独立脚本从 32 个方向
\((1,\varepsilon_1,\ldots,\varepsilon_5)\) 重建所有 32 个仿射 body form，
再以 six section vertices 的正负值筛出 12 个 crossing forms。它将这 12 个 form
与 6 个 chamber facets 合并，精确求解全部
\(\binom{18}{5}=8568\) 个交点系统：3270 个奇异系统、18 个落在闭 chamber 的
不同顶点。

脚本对全部 \(2^{12}=4096\) 个 word 重新建立兼容顶点集合。一个 word 只有在
该集合重心对六条 chamber 不等式和十二条指定 strict sign 都严格成立时才记为
满维。得到 21 个 cell；每个 cell 的**全部**兼容顶点恰为六个。审计逐一确认：

* cell 的 stored rays 与这个由 word 重新得出的顶点集合逐项相同；
* 六射线的齐次 determinant 严格非零，因此是六个仿射独立 section 顶点；
* 所有六条 chamber 不等式、十二条 crossing 的 word 弱不等式都在每个顶点成立；
* 每个原始 32 body form 在整组六顶点取同一弱符号。

这些检查修补了“单独读取 stored rays，或单独读取 certificate crossforms”无法
证明覆盖的缺口。certificate 只在以上重建完成后才用于逐字段相等比较。

## 精确数值

单纯形 section 体积的准确和为
\[
 \sum_{21\ \mathrm{cells}} |\det R|/5! = 1/120,
\]
正是排序 chamber 的体积。对每个射线的 \(\psi\) 值从 32 个 body form 重新计算，
从而得到
\[
 \sum |\det R|/\prod_j\psi(r_j)
 ={520066141\over2247297614217216000},
\]
并以 signed-permutation 群大小 \(2^6 6!\) 得全空间指数积分
\[
 {520066141\over48769479475200}.
\]
因此本 fan 所导出的 \(p_1(0)\) 表达为
\[
 {1\over(2\pi)^6}{520066141\over48769479475200}.
\]

这只是 finite Fourier-fan 的 peak 数据，不是 \(G\) 的 cube 下界、face 积分
上界，也不是 \(C_6\) 新上界。

## 实际运行与范围

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/audit/upper002/independent_fan_audit.py
```

实际耗时约 0.5 秒，使用 Python Fraction；没有展开 967680 个全空间锥，未调用
MATLAB。输出 PASS，审计代码 SHA-256 为
`6a9300f3ea5a667ce3a0969cf841eadcca836b43ab83322d6de3e61c2ebbed66`，结果 SHA-256
见同目录结果文件。生产器的当前源码、同作者 checker 和当前结果分别为
`947664af7de2e6d3d8bbd56a0d7ab36bedecea68fccd89572a48af5d1efce318`、
`234d7bdcd8f458e045c2c2b3e6a499659da61d94cde63e46011fc05637faf06d`、
`18d8aa578f84283d00e1c9c30522fb14cc49812092afd5716f2cc9ae68ba2798`。

主clone集成后将输入路径改为项目相对位置并实际复跑PASS（退出码0，约0.36秒）；可移植脚本SHA及输入SHA以当前independent_fan_audit_results.json为准。上列原SHA保留为独立代理原始证据。
