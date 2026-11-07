#!/usr/bin/env python3
"""Bounded local source scan; continuous-curve dual ceilings are floating only."""
import argparse,importlib.util,json,math,sys,time
from pathlib import Path
import numpy as np
from scipy.special import erf
ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--output',type=Path,default=Path('/tmp/probe61.json'));args=ap.parse_args()
sp=importlib.util.spec_from_file_location('compute60_probe',args.project/'code/compute60.py');m=importlib.util.module_from_spec(sp);sys.modules[sp.name]=m;sp.loader.exec_module(m)
old=json.loads((args.project/'rounds/r060/compute60.json').read_text());w=np.array(old['results'][0]['fixed_weight_refinements']['corrected_weights']);rows=[];start=time.time()
for g in [.15,.10,.05,.025]:
 r=m.run_gamma(g,w);w=np.asarray(r['weights']);v=math.exp(-2*g*16**2/64)+math.sqrt(math.pi*64/(2*g))*erf(16*math.sqrt(2*g/64));d1=r['dual_scaled']-(m.Q-1)*float(m.COEFF@w)
 r.update({'continuous_cone_volume':v,'dual_scaled_target1':d1,'continuous_score_ceiling_float':math.log(v)+m.TAU-math.log(d1) if d1>0 else None})
 rows.append(r);args.output.write_text(json.dumps({'status':'partial','seconds':time.time()-start,'rows':rows,'scope':'Finite-window floating dual only. No interval certificate or continuous feasibility.'},indent=2,default=m.jsonable)+'\n');print(g,r['continuous_score_ceiling_float'],flush=True)
# Preserve the best ceiling candidate and inspect additional intermediate receivers
# under the SAME fixed source, rather than assuming nodal feasibility interpolates.
args.output.write_text(json.dumps({'status':'complete','seconds':time.time()-start,'rows':rows,'scope':'Finite-window floating dual only. No interval certificate or continuous feasibility.'},indent=2,default=m.jsonable)+'\n')
