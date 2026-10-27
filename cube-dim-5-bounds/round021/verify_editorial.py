"""Read-only editorial checks; not a rerun of the numerical certificates."""
from pathlib import Path
from collections import Counter
import re,json,hashlib,subprocess
r=Path(__file__).resolve().parents[1]
old=subprocess.check_output(['git','show','e5b1c4187038104ef17c4ff4e5ee73d14c68a20b:cube-dim-5-bounds/paper/main.tex'],text=True)
new=(r/'paper/main.tex').read_text()
locks=json.loads((r/'round021/source_locks.json').read_text())
assert all(hashlib.sha256((r/p).read_bytes()).hexdigest()==h for p,h in locks.items())
def canon(s):
 s=re.sub(r'\\tag\{[^}]+\}','',s)
 return re.sub(r'\s+','',s.replace('U_{175}','U_*').replace('s_{164}','s_0'))
def displays(s): return re.findall(r'\\\[(.*?)\\\]',s,re.S)
# Reconstruct J from the new line-broken exact decimal-integer representation.
parts=re.search(r'N_J=\{\}&(\d+)\\cdot10\^\{55\}\\\\\s*&\+(\d+),\\\\\s*D_J=\{\}&(\d+)\\cdot10\^\{55\}\\\\\s*&\+(\d+)',new)
assert parts
na,nb,da,db=map(int,parts.groups()); n=na*10**55+nb; d=da*10**55+db
oldJ=re.search(r'J_\*=\s*\\frac\{(\d+)\}\s*\{(\d+)\}',old)
assert (n,d)==tuple(map(int,oldJ.groups()))
a=Counter(map(canon,displays(old))); b=Counter(map(canon,displays(new)))
# The old J display is replaced only for readability; insert its exact old form.
jdisplay=next(x for x in a if 'J_*=\\frac' in x)
b[jdisplay]+=1
missing=list((a-b).elements())
assert len(missing)==1 and 'C_5' in missing[0] and '1.9943' in missing[0] # compressed abstract only
assert all(x in new for x in ['\\section{引言}','\\section{预备知识：半群与有限源比较}','\\section{主要结果的证明}','\\section{应用与进一步问题}'])
assert all(x not in new for x in ['研究现已暂停','双侧突破','PASS','\\mathrm{self}','173候选'])
G=11231606208775453; I=4472478819809024; K=28920678
assert K**4*G**5>3125*I**4*10**44
print(json.dumps({'status':'PASS','scope':'editorial preservation and one scalar display inequality only','old_display_count':sum(a.values()),'old_displays_preserved_modulo_names_tags':46,'abstract_compressed':1,'J_integer_reconstruction':True,'source_locks':locks,'research_resumed':False},ensure_ascii=False,indent=2))
