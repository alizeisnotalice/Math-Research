# 第37轮复现、部分证据与范围

本轮外区常数模型准入未闭合。实际385盒/770组件/6160准确行，311接受叶，首深度4失败后按计划终止，231待处理。不能把部分叶当整父格覆盖或新弱型上界。

Python3.12、python-flint0.9.0。新临时目录中：第20轮scaled_certificate.json放round01；第32轮historical_scaled_second.json.xz解压到根目录scaled_second.json；第20轮公开independent目录按upstream-v4/checker/independent放置。本轮export_cover.py/run_once.py/verify_saved.py放round18。输入均按实际源SHA锁定，无私人附件。

空输出目录运行 `python -B run_once.py`（唯一300秒、1500盒、深度4、192bit），终态后 `python -B verify_saved.py`（准确Fraction保存值复核，内部90秒检查）。不得把已保存结果复制到运行目录再启动；程序拒绝重复。固定常数下界接口顺序(u0,v0,u1,v1)，不调用任何Taylor/T/LP。

若只核验已有结果，将saved_models.tar.xz解压到round18，与result.json、export.run.json、export_cover.py、verify_saved.py配齐，删除新目录中拟输出的saved_verification.json后只运行verify_saved.py；无需重新求模型。它复核全部385模型及6160行、树中点分割和哈希；因为complete=false，不声称做过全父域slab覆盖。核验仍需round01准确核；不是独立重建超越函数区间。

压缩包保存原leaf_cover及385个完整最终模型，便于后续直接复用311接受叶、诊断失败链、复核6160行；不保存重复partial检查点。leaf_cover中的partial文件名/hash仅为历史运行阶段记录，最终复核不读取它们。原始partial与全部工作材料仍在本地。压缩内容已逐文件与本地原件核对。

failure_breakdown.json只从已存模型重组失败node615的最小行、未扣样条余项值、余项及区间宽度，不增加点或模型评价。展示小数不是符号判据。
