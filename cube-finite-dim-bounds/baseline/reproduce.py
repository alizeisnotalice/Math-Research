from pathlib import Path
import json,time,sys
from algebraic_bounds import calculate
p=Path(__file__).resolve().parents[1]
old=json.loads((p/'baseline/original_certificates.json').read_text())
out=[]
for expected in old:
 t=time.monotonic(); actual=calculate(expected['dimension'])
 assert json.loads(json.dumps(actual))==expected, 'certificate mismatch'
 out.append({'dimension':expected['dimension'],'exact_match':True,'seconds':time.monotonic()-t,'witnesses':len(actual['witness_counts'])})
print(json.dumps(out,indent=2))
if '--save' in sys.argv: (p/'results/baseline_run.json').write_text(json.dumps(out,indent=2)+'\n')
