"""Recheck the four archived box volumes, without rerunning the search."""
import json
from pathlib import Path
import numpy as np
import search_positions as s

data=json.loads(Path(__file__).with_name('result.json').read_text())
checks=[]
for row in data['dimensions']:
    d=row['dimension']
    baseline={'spacings':s.base_spacings(d).tolist(),
              'delta_in_spacing_units':np.zeros((4,d-2)).tolist(),
              'score':row['baseline_reproduced']}
    for label,item in [('baseline',baseline),('best_noncontrol',row['best_noncontrol'])]:
        check=s.verify_item(item)
        assert check['status']=='checked' and check['agrees']
        checks.append({'dimension':d,'label':label,**check})
print(json.dumps(checks,indent=2))
