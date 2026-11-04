import json
from pathlib import Path
p=Path(__file__).with_name('coefficients.json')
a=json.loads(p.read_text())
print('cases:',len(a))
print('all series match:',all(x['series_matches'] for x in a))
print('zero dimensions:',[x['n'] for x in a if x['zero']])
print('distinct kernels:',len({x['kernel_id'] for x in a}))
print('distinct pids:',len({x['pid'] for x in a}))
print('n2 expected:',a[0]['conditional_coefficient']=='-16015360/243')
