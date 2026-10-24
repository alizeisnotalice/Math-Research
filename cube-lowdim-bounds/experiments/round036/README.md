# 第36轮：固定B₁的修正版积分路径补验

本轮不搜索新核，固定第20轮scaled_certificate.json（SHA21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172）。只修正公开hlcert.py两处复积分回调的analytic标志传播，原始源与历史输出保持不变。

复现采用Python3.12/python-flint0.9.0。在新临时根目录建round01、round17、upstream-v4；第20轮核放round01，第20轮公开原源码保持upstream-v4/checker层级。补丁生成器及其精确diff、哈希、预算和实际终态见本目录。不得把已有first.run.json/first.log/first_report.json放入新运行目录，包装器拒绝覆盖。

完整运行固定MKA、Q8、A深度7、K深度6、6工人、原程序96bit，硬预算180秒。仅一次；没有在失败后加深、提精度或延长。已保存的独立闭式核验及整数根成本证书承接第20轮，本轮不重新计算它们。

即使此补验通过，也只补齐既有基线的第二条证据链，不是新上界或3.4目标完成。两路径仍共享Arb与解析公式；通用比较和原点外分布接口沿用已审外部定理。

实际唯一运行59.046120416秒正常exit0、PASS；K2628格失败0，A219732格失败0。质量输出不是独立成本证明。first.run.normalized.json仅将命令中的本机绝对路径替换成相对路径；原本地run的SHA单独保留，result.json仍绑定该原始run哈希。首尾源码/核哈希和全部运行字段未变，不把归一化文件的SHA冒充原run。

复现时只复制make_patch.py、run_once.py、approved_inputs.json和summarize_once.py到新的round17；先执行 `python -B make_patch.py`，按生成manifest核对哈希，再执行 `python -B run_once.py` 和 `python -B summarize_once.py`。不要复制已归档patch_manifest/固定副本/结果后再运行生成器；它拒绝覆盖。修正版源由原件和精确补丁重建，避免重复归档整份上游代码。
