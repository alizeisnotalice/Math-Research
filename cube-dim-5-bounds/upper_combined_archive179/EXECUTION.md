# 179A 实际无损归档

唯一实际编码及独立冷恢复均成功；这里仅证明字节传输无损，不新增数学资格。批准的固定计划 SHA 为 `c037e917cf87506e0a7e3e0b0969d9854704ced4e47e36c2de201368ee5e5504`。原计划和所有旧冻结件保持原样；本文件是执行附录。

- controller PID 38003，编码 PID 38005、exit 0、7.4812300839985255 秒；独立恢复 PID 38040、exit 0、7.237960667000152 秒；总计 15.624832083994988 秒。只运行各一次，全部进程已退出，无错误或备用回退。
- 302 literal + 32 原存字段投影，共 334 个 blob。3 包合计 51,863,424 字节；各包 9,834,475 / 29,472,573 / 12,556,376 字节。未压 tar 各为 33,546,240 / 33,546,240 / 28,395,520 字节，均不超过 32MiB；每个物理文件小于 50,000,000 字节。
- 1,525 个逻辑路径、其中 1,173 个原数学 JSON，全部逐原长度和 SHA 精确恢复。95 个已 tracked 的仓库依赖只核字节 SHA，不再次携带。全部编码源前后 SHA 一致。
- 恢复器不打开原 producer scratch；H 二进制解码与数学算法调用均为 0。恢复通过复制已存字段、按 key 查表和原顺序序列化完成，不计算 factor、mean、root、slope、gain、floor 或 DFS。

## 最小导入

`delivery_whitelist.json` 为逐文件源路径、目标相对路径、长度和 SHA。只复制其列出的文件，再带同目录的 `delivery_whitelist.json`、`execution_manifest.json`。三包和 control 均放在 `packet/`，证据放在 `packet/evidence/`。不导入 staging、明文 tar、恢复沙盒、blob cache、进度或操作 stdout/stderr。162/164/166 的原源码、证明、授权、注册、master 和各15退出码及相应独审证据，已作为所选逻辑文件保存在包内。

不重复携带 167 原 H 或 opaque 169 原 tar，也不携带 180 大输入或 journal。已有依赖路径及各原长度/SHA 见 `packet/control/repository_refs.json`；解码需要含这些原字节的研究仓库。单纯 raw H 文件 SHA 核验不构成 H 真值区间的新证明。

## 便携字节恢复入口

在有 Python 3 标准库的机器，将下列三个绝对路径占位换成实际位置。输出目录必须尚不存在；保留其输出作为验证沙盒，无须提交该目录。此命令只供读者复核，本次实际运行已完成，不应为了补记录再运行。

```sh
python3 -B /ABS/179/packet/control/decode179.py \
  --archive /ABS/179/packet \
  --manifest-sha 8b8a27676cad711ffa6ae49792a3abd7718aea8ddd8e4fdbffca06bfb01c83f0 \
  --repository /ABS/Math-Research-cube-dim-5 \
  --output /ABS/NEW_EXCLUSIVE_RESTORATION_DIRECTORY
```

解码会核 control、三包、334 blob、所有原逻辑文件和95仓库引用。实际逐文件全恢复证据是 `packet/evidence/restoration_checks.json`，SHA `ab14d41d7317d9a559dd68f74c788f045202fda09d81262ae2d4af6aca8fada9`；实际解码 receipt SHA `ce3138c8601d4c16d5420b7ff766078dc890740591b7e1c78be34356d527637a`。

原 prepared gzip 缓存的历史字节不在恢复承诺内。完整数值字段由共享169表、小投影及 tracked167 索引保留；有新 SHA 的输入视图不得冒充旧缓存 SHA。此边界与所有1,173个原数学JSON必须原字节恢复的承诺不同，后者已经全部通过。未修改旧冻结 checker 的 SHA 门，不声称恢复后直接绕过其历史缓存接口。

未提交或推送。上界数学接纳和论文结论仍由主负责人已有的严格审计链承担。
