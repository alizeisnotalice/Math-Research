"""Independent exact integer inclusion-exclusion and prior actual recursion checks."""
from pathlib import Path
from itertools import combinations
import random,json,time,numpy as np
from round007_union_fast import volume
from round007_three_orbit import evaluate,A,face_keys
from round005_exposure import geometry
P=Path(__file__).resolve().parent

def main():
    rng=random.Random(20261025710)
    for _ in range(32):
        b=[]
        for k in range(6):
            z=[]
            for j in range(4):
                lo=rng.randrange(-3,4);z.extend([lo,lo+rng.randrange(1,5)])
            b.append(z)
        exact=0
        for n in range(1,7):
            for ss in combinations(b,n):
                v=1
                for j in range(4):v*=max(0,min(x[2*j+1] for x in ss)-max(x[2*j] for x in ss))
                exact+=(-1)**(n+1)*v
        assert volume(b)[0]==exact
    obs=json.loads((P/'round007_recursive_observations.json').read_text());models=json.loads((P/'round007_model_results.json').read_text())['rows'];patches,_,_=geometry();rows=[];first=models[0];pool,_=face_keys([np.array(first['phi1']),np.array(first['phi2'])],patches)
    for old in obs['observations']:
        model=next(x for x in models if x['name']=='A_complement' and abs(x['share1']-old['share1'])<1e-10);s=old['sigma'];row,keys,selected=evaluate(A*(1+s*np.array(model['direction'])),np.array(model['phi1']),np.array(model['phi2']),s,old['share1'],pool=pool);gap=row['ratio']-old['ratio'];assert abs(gap)<1e-12
        rows.append(dict(sigma=s,share1=old['share1'],old_ratio=old['ratio'],sweep_ratio=row['ratio'],difference=gap,sweep_seconds=row['seconds'],old_seconds=old['seconds']));print(rows[-1],flush=True)
    (P/'round007_sweep_validation.json').write_text(json.dumps(dict(exact_integer_six_box_fixtures=32,all_passed=True,actual_recursive_crosschecks=rows),indent=2)+'\n')
if __name__=='__main__':main()
