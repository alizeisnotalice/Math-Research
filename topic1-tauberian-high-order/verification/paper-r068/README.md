# 第68轮正式正文：源、实际编译与视觉验收

交付七源 identity `42eead43fe9e13ca26d0ea30aa6165c1e845c3dfac006d5a64c7107f8ad63647`。PDF为73页、703251B，SHA `5f84aaa0a3d18fad3e9a86eea5a275bf9c5222a13d8171f2d10a7e47d6ff98a1`。本包不把编译通过当数学证明；[独立数学重建](../../reviews/068-paper-math-audit.md)和[当前语义节点](../../evidence/068-dependencies.json)另给实际推导。

## 可逆来源与复现

`python3 verify_source_chain.py`实际逐字逆核三阶段：final-trim→layout→math-freeze→旧67页源。两个新增文件的旧状态是不存在（manifest的null），并非旧版空文件。所有检查在内存完成，不改论文。末尾机械b7表删除、层数/前缀符号修正与三处排版修补均有明确旧新文本。

在`../../paper/`运行`latexmk -norc -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex`。本机实际使用TeX Live 2026/XeLaTeX；需要Songti SC、Times New Roman及标准mathrsfs包。`compile-delivery-receipt.json`给原实际命令、exit0、源及PDF哈希；`compile-delivery.log`是原最终main.log。完整工具结构化返回保留在本地，SHA载于receipt，不重复归档它内含的长日志。未使用早前6126或中间7a46版的编译记录代最后5f84版。

## 实际渲染与检查

最终73页均用pdftoppm在15个实际worker渲染至100dpi；`render/render-result.json`保每页尺寸、worker PID、时间。与此前已全页看过的布局版本比较，仅36/47/72/73页改变，另外69页PNG逐字一致；四页另150dpi查看通过。主代理亦实际看了所有73页联系表、多个150dpi关键页及最终改变的四页，未见裁切、碰撞或异常空页。最终编译零告警、零hbox/缺字，225标签/256引用全部解析。

如需独立重渲染，将`render/render_all.py`和`../../paper/tauberian.pdf`复制到一个临时目录，PDF命名`input.pdf`后运行该Python脚本（需Pillow、pdfinfo、pdftoppm）。它检查固定输入SHA和73页，实际批次脚本原样保留。整页PNG及此前比较所需旧PNG保存在仓库外；`comparison.json`是实际比较记录，比较脚本因依赖旧PNG未列为独立核心。

`o9_algebra.py/stdout`是数学独审实际Fraction窄核，覆盖指定n的移动峰、N9层价和真实分母；它不证明全输入几何。原始编译/审报保留执行时机器路径，不作为运行依赖。`manifest.json`绑定本包、数学审报与交付源/PDF。
