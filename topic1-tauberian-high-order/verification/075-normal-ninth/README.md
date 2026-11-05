# R075 实际计算与逆补丁

本包只复现固定竞争常数下的正常第三反馈有限账本与字节级逆补丁；解析几何和独审见 `../../proofs/075-normal-ninth.md`、`../../reviews/075-normal-ninth.md`。局部源/退款估计为 O(v^10 S)=o(v^9 S)，不识别真实 C_n 九阶系数，也不主张全输入九阶渐近。

## 文件与运行

- `third_feedback.wl` / `third_feedback.log`：实际 Wolfram 15 ARM PID 86683；核有限递推幂次，不生成额外 JSON。
- `independent_checks.py` / `independent_checks.json`：标准库 Fraction 独立有限核；无输入文件，原输出已保存。625 个有理输入通过；包含 normal loss=1/2、full clipped loss=0 的同 eroder 反例。
- `shared_eroder.diff` / `shared_eroder_check.json`：v1→v2 两处替换及原作者逆恢复记录。旧长源不重复归档。
- `verify_inverse.py` / `verify_inverse.json`：本包新增的可携带逆校验。它读取本轮 public proof、diff 与 check JSON，在内存逆应用两个 unified-diff hunk，并比对原源 SHA/字节数；不写 source。执行：

```sh
python3 verify_inverse.py > verify_inverse.json
```

该校验器无用户绝对路径。原作者 `shared_eroder_patch.py` 含本机绝对路径，未纳入公开包；本轮提供相对路径的只读逆校验，未运行该原始写文件脚本。

## 范围

计算只核代数账本和历史补丁字节身份，不证明有效源几何、CDF 或 Tauberian 结论。文件逐项哈希见包内 `manifest.json`（路径相对项目根目录）。
