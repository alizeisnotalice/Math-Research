#!/usr/bin/env python3
"""Build a source-locked portable handoff, excluding caches/private files."""
import argparse, hashlib, json, subprocess, zipfile
from pathlib import Path
BASE='77be2e97b708bbe58f365b291983a78a32d7ea1d'
VERIFY='''#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,sys
root=Path(__file__).resolve().parent
m=json.loads((root/'PACKAGE_MANIFEST.json').read_text(encoding='utf-8'))
errors=[]
for e in m['files']:
 p=root/e['path']
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=e['sha256']:
  errors.append(e['path'])
if errors:
 print('FAIL:',*errors,sep='\\n');sys.exit(1)
print('PASS:',len(m['files']),'files; source commit',m['source_commit'])
print('Integrity check only: not a proof check or new mathematical result.')
'''
def main():
 a=argparse.ArgumentParser();a.add_argument('output',type=Path);args=a.parse_args()
 root=Path(__file__).resolve().parents[2]
 names=subprocess.check_output(['git','ls-tree','-r','--name-only',BASE,'--','cube-finite-dim-bounds','AGENTS.md'],cwd=root,text=True).splitlines()
 files={name:subprocess.check_output(['git','show',f'{BASE}:{name}'],cwd=root) for name in names}
 for name in ('README_资料包说明.md','新对话启动提示词.md','build_package.py'):
  f=Path(__file__).parent/name;files['cube-finite-dim-bounds/handoff/'+name]=f.read_bytes()
 files['verify_package.py']=VERIFY.encode()
 files['START_HERE.md']='''# 高维研究交接包
先运行 `python verify_package.py`，再阅读 `cube-finite-dim-bounds/handoff/README_资料包说明.md`。
把同目录 `新对话启动提示词.md` 中 `<目标维度>` 替换为5、6、8等单个整数，连同本ZIP交给新对话。
原项目保持暂停；新维度使用独立项目。文件完整性检查不等于数学认证。
'''.encode()
 manifest={'date':'2026-10-25','source_commit':BASE,'scope':'Complete tracked cube-finite-dim-bounds snapshot plus handoff instructions; no other research projects or private originals','files':[{'path':k,'bytes':len(v),'sha256':hashlib.sha256(v).hexdigest()} for k,v in sorted(files.items())]}
 files['PACKAGE_MANIFEST.json']=(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n').encode()
 args.output.parent.mkdir(parents=True,exist_ok=True)
 with zipfile.ZipFile(args.output,'w',zipfile.ZIP_DEFLATED,compresslevel=9) as z:
  for name,data in sorted(files.items()):
   info=zipfile.ZipInfo(name,date_time=(2026,10,25,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,data)
 print(json.dumps({'path':str(args.output),'files':len(files),'bytes':args.output.stat().st_size,'sha256':hashlib.sha256(args.output.read_bytes()).hexdigest()},ensure_ascii=False))
if __name__=='__main__':main()
