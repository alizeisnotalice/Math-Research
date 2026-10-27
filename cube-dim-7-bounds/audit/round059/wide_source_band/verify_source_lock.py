#!/usr/bin/env python3
"""Definition/provenance-only source lock; no mathematical computation rerun."""
import hashlib
import json
import re
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round059/wide_source_band')
LOCAL=Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round056/one_coordinate/LEMMA.md')
OLD='06a3ea9ac8e754ca0f3e63664aea99de20cd49a9c9080c49c8ceb3312e35f8ff'
NEW='fe95f369634bbc753f5e44a32d6a2b74b288748ebd4bfe1f6c3cac43b086bd8f'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
expected={ROOT/'PROOF.md':NEW,
    ROOT/'REVISION.patch':'8b2764ab61fa102df12ee17a6b860f7c8059c094befbc37275b8d36d53cde3f7',
    ROOT/'DEPENDENCY_ADAPTATION.md':'8e97ecbdf6c8ca33eedf06116617b1c6fb51edc5821433a150387ff856061429',
    ROOT/'certificate.json':'291fc41ee7fa40294992fe4e77ef011ebfa1f6657d8e0846278a8f3286a41421',
    ROOT/'certify.py':'332cd3ae1cc7a2094a5ffa61d9cbeebd7cd0a9ea319c7410bad8532fb93fea28',
    LOCAL:'97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec'}
assert all(sha(path)==digest for path,digest in expected.items())

newlines=(ROOT/'PROOF.md').read_text().splitlines(keepends=True)
patch=(ROOT/'REVISION.patch').read_text().splitlines(keepends=True)
oldlines=[];cursor=0;pos=0;hunks=0
while pos<len(patch):
    line=patch[pos]
    if not line.startswith('@@ '):
        pos+=1;continue
    match=re.match(r'@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@',line)
    assert match
    newstart=int(match.group(3))-1
    assert newstart>=cursor
    oldlines.extend(newlines[cursor:newstart]);cursor=newstart
    pos+=1;hunks+=1
    while pos<len(patch) and not patch[pos].startswith('@@ '):
        piece=patch[pos]
        if piece.startswith(' '):
            assert newlines[cursor]==piece[1:]
            oldlines.append(piece[1:]);cursor+=1
        elif piece.startswith('+'):
            assert newlines[cursor]==piece[1:]
            cursor+=1
        elif piece.startswith('-'):
            oldlines.append(piece[1:])
        elif piece.startswith('\\'):
            pass
        else:
            raise AssertionError('unexpected diff line')
        pos+=1
oldlines.extend(newlines[cursor:])
reverse_hash=hashlib.sha256(''.join(oldlines).encode()).hexdigest()
assert reverse_hash==OLD and hunks==4
old_receipt=json.loads((HERE/'new_margins_receipt.json').read_text())
assert old_receipt['sha256'][str(ROOT/'PROOF.md')]==OLD

files=list(expected)+[HERE/'new_margins_receipt.json',HERE/'check_new_margins.py',
    HERE/'review.md',Path(__file__).resolve()]
out={'status':'FINAL_SOURCE_LOCK_PASS','old_proof_sha256':OLD,'final_proof_sha256':NEW,
    'reverse_patch_reconstructs_initial_proof':True,'hunk_count':hunks,
    'mathematical_estimates_or_parameters_changed':False,'new_math_run':False,
    'original_new_margin_receipt_preserved':True,
    'sha256':{str(path):sha(path) for path in files},
    'scope':'Incremental source/definition/provenance check only. Full analytic acceptance is in review.md. No old arithmetic/grid/baseline repeated and no explicit n threshold created.'}
(HERE/'source_lock.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'final_proof_sha256':NEW,'reverse_patch_reconstructs_initial_proof':True,'hunk_count':hunks},indent=2))
