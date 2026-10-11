# A04 本地来源解析

索引中的路径相对于交付包根目录。`EVIDENCE_ROOT` 指包含 `workspace_revised/`、`audit_current/` 与 `evidence/` 三个目录的包根，不是 `evidence/` 子目录。不要从已安装 Skill 目录拼接 `../../evidence`。

在解压后的完整交付包中，可显式指定包根并解析 P-56 的索引条目 `E1952`：

```sh
export EVIDENCE_ROOT="/path/to/unpacked/数学工具64项_完整验收提升_20261007"
python3 "$EVIDENCE_ROOT/audit_current/delivery/resolve_evidence.py" \
  --package-root "$EVIDENCE_ROOT" \
  --skill-dir "$EVIDENCE_ROOT/workspace_revised/math-a04-shared-nonanticipative-mip" \
  --evidence-id E1952
```

已安装的 Skill 可通过用户本地 root registry 查找同一包。registry 应符合 `math-skill-evidence-root-registry-v1`，并将 `math64-20261007` 指向解压后的包根；模板见 `audit_current/delivery/source-root-registry.example.json`。以下调用假设该包仍在 `MATH64_SOURCE_PACKAGE_ROOT`，该变量用于找到随包提供的校验器；为使 `--root-id` 从 registry 解析根目录，请不要同时设置 `EVIDENCE_ROOT`：

```sh
unset EVIDENCE_ROOT
export MATH64_SOURCE_PACKAGE_ROOT="/path/to/unpacked/数学工具64项_完整验收提升_20261007"
python3 "$MATH64_SOURCE_PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py" \
  --root-id math64-20261007 \
  --registry "$HOME/.codex/math-skill-evidence-roots.json" \
  --skill-dir "$HOME/.codex/skills/math-a04-shared-nonanticipative-mip" \
  --evidence-id E1952
```

`E1952` 指向 P-56 PDF，SHA-256 为 `56e47d0d5a099c44df71d7e1246fbdaaa784b1afba8261d9795af85c657dbe5e`。校验器仅在路径安全且文件 SHA 与索引一致时返回本地路径；通过解析不表示来源已阅读、结论已证明或 Skill 已安装。用户本地 registry 及 Skills 安装仅在单独获授权后配置。
