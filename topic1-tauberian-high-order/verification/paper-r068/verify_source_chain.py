"""Check the delivery bytes and reverse all three editorial stages in memory."""
import hashlib,json,sys
from pathlib import Path
package=Path(__file__).resolve().parent
paper=Path(sys.argv[1]) if len(sys.argv)>1 else package.parents[1]/'paper'
evidence=Path(sys.argv[2]) if len(sys.argv)>2 else package
sha=lambda b:hashlib.sha256(b).hexdigest()
freeze=json.loads((evidence/'delivery-freeze.json').read_text())
state={row['path']:(paper/row['path']).read_bytes() for row in freeze['files']}
for row in freeze['files']:
    assert sha(state[row['path']])==row['sha256'],row['path']
assert sha((paper/'tauberian.pdf').read_bytes())==freeze['pdf_sha256']
report={'delivery_identity':freeze['identity_sha256'],'pdf_sha256':freeze['pdf_sha256'],'stages':[]}
for name,key,new_key,old_key in [
 ('final-trim-manifest.json','files','new_sha256','old_sha256'),
 ('layout-manifest.json','files','layout_sha256','math_freeze_sha256'),
 ('source-edit-manifest.json','source_files','new_sha256','old_sha256')]:
    manifest=json.loads((evidence/name).read_text());rows=[]
    for row in manifest[key]:
        path=row['path'];assert sha(state[path])==row[new_key],(name,path,'new')
        lines=state[path].decode('utf-8').splitlines(keepends=True)
        for op in reversed(row['operations']):
            start,end=op['new_start_line0'],op['new_end_line0']
            assert ''.join(lines[start:end])==op['new_text'],(name,path,start)
            lines[start:end]=op['old_text'].splitlines(keepends=True)
        old=''.join(lines).encode('utf-8')
        if row[old_key] is None:
            assert old==b'',(name,path,'absent old file')
            state[path]=None;rows.append({'path':path,'old_sha256':None,'bytes':None})
        else:
            assert sha(old)==row[old_key],(name,path,'old')
            state[path]=old;rows.append({'path':path,'old_sha256':sha(old),'bytes':len(old)})
    report['stages'].append({'manifest':name,'files':rows,'exact_reverse':True})
print(json.dumps(report,ensure_ascii=False,indent=2))
