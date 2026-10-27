#!/usr/bin/env python3
import json,math
from pathlib import Path
p=Path(__file__).resolve().parent;r=json.loads((p/'results.json').read_text());out=[]
for i in range(3):
 a=[x for x in r['results'] if x['handle'].startswith(f'p{i}')];n=sum(x['n'] for x in a); sums=[sum(x['sum'][j] for x in a) for j in range(7)]; square=[sum(x['sum_squares'][j] for x in a) for j in range(7)];seeds=[[x['sum'][j]/x['n'] for j in range(7)] for x in a];tail=[[sum(x['kernel_bin_sums'][j][k] for x in a for k in (7,8))/sums[j] for j in range(4)]];u=[sum(x['u_bin_sums'][k][j] for x in a for k in range(3))/sums[j] for j in range(7)]
 out.append(dict(point_id=i,n=n,mean=[v/n for v in sums],sample_variance=[max(0,(square[j]-sums[j]**2/n)/(n-1)) for j in range(7)],seed_mean_min=[min(x[j] for x in seeds) for j in range(7)],seed_mean_max=[max(x[j] for x in seeds) for j in range(7)],seed_means=seeds,shape_term_ge_0_01_contribution_fraction=tail[0],U_lt_0_001_contribution_fraction=u,high_precision_recomputations=sum(len(x['fallbacks']) for x in a)))
(p/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
