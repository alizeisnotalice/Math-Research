#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
import satellite_screen as S
HERE=Path(__file__).resolve().parent
eps=1/1024
S.A=S.A*(1+eps)**.2
S.HALF=S.A/2
C=json.loads((HERE.parent/'lower_round003/a2_equals_a3_results.json').read_text())['display_only'][0]
rows=[]
for first in range(2):
 for k in range(5):
  if first==k==0:continue
  row=S.run(first,k,eps)
  row['difference_float']=row['ratio_float']-C
  row['finite_difference_response']=row['difference_float']/eps
  rows.append(row)
  print(json.dumps(row),flush=True)
  (HERE/'codilation_results.json').write_text(json.dumps({'status':'floating_diagnostic_only','rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
