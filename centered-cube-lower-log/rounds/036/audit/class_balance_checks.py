"""Exact box-crossing/category-mass audit of saved compensation jitters."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
src=ROOT/'compensation/source_jitters.json'
rows=[]
for source in json.loads(src.read_text()):
    K=source['K'];U=[F(x)for x in source['U']];R=3
    assert K%R==0
    for v in [8,32,128]:
        width=F(1,4*v)
        mass=[F(1,R)for _ in range(R)]
        left=[];right=[]
        for j,u in enumerate(U):
            if u<width/2:
                f=(width/2-u)/width/K;mass[j%R]-=f;mass[(j-1)%R]+=f;left.append(j)
            if 1-u<width/2:
                f=(width/2-(1-u))/width/K;mass[j%R]-=f;mass[(j+1)%R]+=f;right.append(j)
        assert sum(mass,F(0))==1
        err=max(abs(x-F(1,R))for x in mass)
        rows.append({'K':K,'v':v,'R':R,'left_crossings':left,'right_crossings':right,
                     'mu_category_masses_exact':list(map(str,mass)),
                     'max_mu_class_error_exact':str(err),'max_mu_class_error_display':float(err)})
out={'scope':'Exact dyadic-U box transfer arithmetic only; no simulation.',
     'source_jitters_SHA256':hashlib.sha256(src.read_bytes()).hexdigest(),
     'rows':rows,'conclusion':'Low-frequency zero category integral preserves an already balanced base, but cannot fix these nonzero base imbalances.'}
(HERE/'class_balance_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print('Nine exact box/category checks saved; nonzero imbalances in',sum(x['max_mu_class_error_exact']!='0'for x in rows),'cases')
