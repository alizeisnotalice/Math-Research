# -*- coding: utf-8 -*-
"""R47: 修订 P-2f81500e5dfa70c6 卡中两处非原文条件表述（append-only，保留原文）。"""
import json, hashlib, shutil, os

paths = [
 'workspace_revised/math-i02-filtration-projections-compensators/references/papers/HK-P-2f81500e5dfa70c6.json',
 'workspace_revised/math-j01-inhomogeneous-jump-generator/references/papers/HK-P-2f81500e5dfa70c6.json',
 'workspace_revised/math-j05-frozen-coordinate-iterated-compensation-interface/references/papers/HK-P-2f81500e5dfa70c6.json',
]

old_assump = "A is Borel and 0∉A, hence is bounded away from zero."
new_assump = ("A is Borel and 0∉A.【修订 2026-10-06 R47（卡↔PDF 比对）】PDF p4 原文仅假设 0∉A；"
 "卡原句\"hence is bounded away from zero\"系错误推断（0∉A 推不出，反例 A={1/n}）。"
 "原文证明中 \"lim ΔX_{t_n} ∈ Ā ⇒ ΔX_{t*}≠0 (since 0∉A)\" 一步实际需要 0∉closure(A)（即 A 与 0 有正距离），"
 "该条件为闭合原文证明 gap 所需的附加假设，使用本命题时应显式补充；卡另一定义卡已用正确的 0∉closure(A) 表述。")

old_def = "for cadlag X and A bounded away from zero, only finitely many A-jumps occur on each compact time interval."
new_def = ("for cadlag X and A bounded away from zero, only finitely many A-jumps occur on each compact time interval."
 "【修订 2026-10-06 R47】原文 Lemma 2（PDF p5）仅假设 A∈B*（0∉A），未写 bounded away from zero；"
 "数学上 0∉A 不足以保证结论（A={1/n} 可有无穷多 A-jumps），本卡条件是必要的数学修补而非原文转录，"
 "原文 p5 \"A⊆R\\(−ε,ε) for all ε>0\" 的量词系原文疏漏（应为\"存在 ε\"）。")

for p in paths:
    c = json.load(open(p))
    changed = []
    for d in c.get('definitions', []):
        if d['name'] == 'jump measure' and old_def in d['statement']:
            d.setdefault('statement_original_r47', d['statement'])
            d['statement'] = d['statement'].replace(old_def, new_def)
            changed.append('definition:jump measure')
    for t in c.get('theorem_cards', []):
        if 'Proposition 1' in t.get('label', ''):
            new_list = []
            hit = False
            for a in t['assumptions']:
                if a == old_assump:
                    t.setdefault('assumptions_original_r47', list(t['assumptions']))
                    new_list.append(new_assump)
                    hit = True
                else:
                    new_list.append(a)
            if hit:
                t['assumptions'] = new_list
                changed.append('theorem_card:Proposition 1')
    if changed:
        json.dump(c, open(p, 'w'), ensure_ascii=False, indent=1)
        h = hashlib.sha256(open(p, 'rb').read()).hexdigest()[:16]
        print(os.path.relpath(p), '->', changed, 'sha256[:16]=', h)
    else:
        print(os.path.relpath(p), '-> NO MATCH (需人工检查)')
