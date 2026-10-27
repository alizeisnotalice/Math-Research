#!/usr/bin/env python3
import json,math
from pathlib import Path
p=Path(__file__).resolve().parent;d=json.loads((p/'results.json').read_text());a=d['results'];n=sum(r['completed_shapes'] for r in a);out=[]
for k,name in enumerate(d['labels']):
 s=sum(r['sum'][k] for r in a);s2=sum(r['sum_squares'][k] for r in a);mean=s/n;var=(s2-s*s/n)/(n-1);assert var>=0
 sm=[r['sum'][k]/r['completed_shapes'] for r in a];out.append(dict(label=name,mean=mean,sample_variance=var,sample_sd=math.sqrt(var),descriptive_sem=math.sqrt(var/n),minimum=min(r['minimum'][k] for r in a),maximum=max(r['maximum'][k] for r in a),stream_means=sm,stream_mean_min=min(sm),stream_mean_max=max(sm)))
z=dict(status=d['status'],sample_size=n,stats=out,raw_violations=[sum(r['raw_violations'][j] for r in a) for j in range(8)],tolerance_violations=[sum(r['tolerance_violations'][j] for r in a) for j in range(8)],fallback_count=sum(len(r['fallbacks']) for r in a),interpretation='finite pseudo-random sample diagnostics; SEM and block variability are descriptive, not deterministic error certificates')
(p/'summary.json').write_text(json.dumps(z,indent=2)+'\n')
