"""Read-only source identity, references and baseline-delta audit; no gate reruns."""
import difflib
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

HERE=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/cube-cover-manuscript/cube-dim-7-bounds')
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7')
NOTES=AUTHOR/'upper/round077/manuscript_notes'
digest=lambda data: hashlib.sha256(data).hexdigest()
source=(AUTHOR/'paper/main.tex').read_bytes()
assert digest(source)=='989d28ff83bd51f2c2dfe1822d9ecb59dd71502e5deb361db845cd87359edc9d'
text=source.decode()
oldline='这里相对尺度 $k$ 中的 $r_{\\rm hi}$ 不可遗漏。'
newline='这里 $k$ 的归一化包含 $r_{\\rm hi}$。'
assert text.count(newline)==1 and oldline not in text
recovered=text.replace(newline,oldline).encode()
assert digest(recovered)=='7f37e6594779f2223ab2f560bbf67de7d4c243fc83b116f478dbd2ab8faf1105'

baseline=subprocess.check_output(['git','cat-file','blob','d4a09760b285ebc9fc39cbfd87d90aff2a7aed18'],cwd=MAIN)
assert digest(baseline)=='8135adfb51470228ed18371c6c5f1a22caac0b35523b0cd49c50b13cfb9b386b'
old=baseline.decode()
labels=re.findall(r'\\label\{([^}]+)\}',text)
references=re.findall(r'\\(?:eqref|ref|pageref|autoref)\{([^}]+)\}',text)
keys=re.findall(r'\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}',text)
cites=[v.strip() for group in re.findall(r'\\cite(?:\[[^\]]*\])*\{([^}]+)\}',text) for v in group.split(',')]
assert len(labels)==len(set(labels)) and set(references)<=set(labels)
assert len(keys)==len(set(keys)) and set(cites)<=set(keys)
assert not re.search(r'round0|SHA256|[0-9a-f]{64}',text)
stack=[]
for match in re.finditer(r'\\(begin|end)\{([^}]+)\}',text):
    op,env=match.groups()
    if op=='begin': stack.append(env)
    else:
        assert stack and stack.pop()==env, (op,env)
assert not stack

def labelled_displays(s):
    ans={}
    for match in re.finditer(r'\\begin\{(equation\*?|align\*?|gather\*?)\}.*?\\end\{\1\}',s,re.S):
        for label in re.findall(r'\\label\{([^}]+)\}',match.group()): ans[label]=match.group()
    return ans
oldmath,newmath=labelled_displays(old),labelled_displays(text)
retained={label: oldmath[label]==newmath.get(label) for label in oldmath}
assert all(retained.values()), [k for k,v in retained.items() if not v]
wording=json.loads((NOTES/'wording_changes.json').read_text())
for entry in wording:
    assert entry['new'] in text and entry['old'] not in text, entry

sources={}
for group in ['sources.json','supplement_hashes.json']:
    entries=json.loads((NOTES/group).read_text())
    for key,entry in entries.items():
        path=Path(entry['source'])
        value=digest(path.read_bytes())
        assert value==entry['sha256'], (key,value,entry['sha256'])
        sources[key]={'path':str(path),'sha256':value}
primary=AUTHOR/'revision_v3/literature/originals'
for name in ['stein-stromberg1983.pdf','Aldaz2017_1605.05596v2.pdf','aldaz2020v2.pdf']:
    p=primary/name
    sources['primary/'+name]={'path':str(p),'sha256':digest(p.read_bytes())}
preflight=digest((HERE/'adaptation_identity_receipt.json').read_bytes())
delta=[]
for op,i,j,k,l in difflib.SequenceMatcher(None,old.splitlines(),text.splitlines(),autojunk=False).get_opcodes():
    if op!='equal': delta.append({'op':op,'old_start':i+1,'old_lines':j-i,'new_start':k+1,'new_lines':l-k})
out={
 'status':'PASS','pid':os.getpid(),'source_sha256':digest(source),'lines':len(text.splitlines()),
 'recovered_initial_source_sha256':digest(recovered),
 'R068_git_blob':'d4a09760b285ebc9fc39cbfd87d90aff2a7aed18','R068_sha256':digest(baseline),
 'labels':len(labels),'reference_occurrences':len(references),'bibliography_keys':keys,'citation_occurrences':len(cites),
 'duplicate_labels':[],'undefined_references':[],'undefined_citations':[],
 'balanced_environments':True,'unchanged_old_labelled_display_blocks':len(retained),
 'wording_replacements_present':len(wording),'baseline_delta_ranges':delta,
 'sources':sources,'preflight_receipt_sha256':preflight,
 'script_sha256':digest(Path(__file__).read_bytes()),
 'scope':'Static final-source and evidence identity run only. No author14, author53, old13, old local/math verifier, search, numerical optimization, MATLAB, or compilation executed.'
}
(HERE/'frozen_source_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['status','pid','source_sha256','labels','reference_occurrences','bibliography_keys','unchanged_old_labelled_display_blocks','wording_replacements_present']},indent=2))
