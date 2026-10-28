"""Numerical planning only; q=9/29 is a reference, not measured certification."""
from pathlib import Path
import math,json
rows=[]
for K in [1023,4095,16383]:
 a=.5/K;q=9/29;n=9e8*(K/1023)**2
 for _ in range(100):n=16*math.log(n)/(a*a*q)
 rows.append(dict(K=K,planning_n=n,planning_v=16*math.log(n),planning_n_a2=n*a*a,planning_M=math.ceil(math.sqrt(16*math.log(n)))))
Path(__file__).with_name('planning.json').write_text(json.dumps(rows,indent=2))
