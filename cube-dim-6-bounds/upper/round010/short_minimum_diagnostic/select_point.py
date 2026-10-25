#!/usr/bin/env python3
"""Rational projected linear forecasts; no spatial kernel evaluation."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,time
import flint
from flint import arb,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
oldpath=PROJECT/'upper/round007/g_time_moments/results.json';old=json.loads(oldpath.read_text())
gr=json.loads((HERE/'gradient_results.json').read_text());pre=json.loads((PROJECT/'upper/round009/short_potential/postprocess_results.json').read_text())
assert gr['status']=='COMPLETE_TWO_NEW_GRADIENT_PASSES' and gr['new_full_passes']==2
c=[F(v,gr['denominator'])for v in gr['center_numerators'][1:]];D=10**12

def candidate(r,g):
 norm=max(map(abs,g));ys=[max(F(1,10**8),min(F(1)-F(1,10**8),x-r*a/norm))for x,a in zip(c,g)]
 ys=sorted(ys,reverse=True);ns=[min(D-10000,max(10000,round(y*D)))for y in ys]
 for i in range(1,5):
  if ns[i]>=ns[i-1]:ns[i]=ns[i-1]-1
 assert all(0<ns[i]<D for i in range(5)) and all(ns[i]>ns[i+1]for i in range(4))
 return [D]+ns

def forecast(ns,g,k):return k+sum(float(a)*float(F(n,D)-x)for a,n,x in zip(g,ns[1:],c))
rows=[]
for row in gr['gradient_passes']:
 T=F(row['tau']);late=[arb(z)for z in row['gradient_x2_to_x6']];vecs=[]
 for a in old['gradients'][:3]:vecs.append([arb(z)-g for z,g in zip(a['gradient_x2_to_x6'],late)])
 vs=[[float(z)for z in v]for v in vecs];base=vs[0];n=math.sqrt(sum(x*x for x in base));norm=max(map(abs,base))
 cos=[sum(a*b for a,b in zip(base,v))/(n*math.sqrt(sum(x*x for x in v)))for v in vs]
 dif=[max(abs(a-b)for a,b in zip(base,v))/norm for v in vs]
 strong=[i for i,x in enumerate(base)if abs(x)>.02*norm]
 signs_ok=all(all(base[i]*v[i]>0 for i in strong)for v in vs)
 stable=signs_ok and min(cos)>.995 and max(dif)<.05
 g=[F(round(x*10**12),10**12)for x in base]
 k=float(arb(next(x for x in pre['rows']if F(x['T'])==T)['K_T_center']['ball']))
 target=.98*float(T/4096);maxr=F(3,20);r=None
 if stable and forecast(candidate(maxr,g),g,k)<=target:
  lo,hi=0,150000
  while hi-lo>1:
   mid=(lo+hi)//2
   if forecast(candidate(F(mid,10**6),g),g,k)<=target:hi=mid
   else:lo=mid
  r=F(hi,10**6)
 ns=None if r is None else candidate(r,g)
 rows.append({'T':str(T),'band_gradients_arb384':[[z.str(80,more=True)for z in v]for v in vecs],
 'cosines_float_diagnostic':cos,'relative_max_differences_float_diagnostic':dif,'stable':stable,
 'rounded_rational_direction_gradient':list(map(str,g)),'line_step':None if r is None else str(r),
 'linear_target_98pct_threshold':target,'forecast_at_max_step':forecast(candidate(maxr,g),g,k),
 'candidate_numerators':ns,'denominator':D,'linear_K_forecast':None if ns is None else forecast(ns,g,k),
 'center_K_float':k,'strong_component_signs_agree':signs_ok})
valid=[r for r in rows if r['line_step'] is not None]
if valid:
 best=min(valid,key=lambda r:F(r['line_step']));selected={k:best[k]for k in ['T','line_step','rounded_rational_direction_gradient','candidate_numerators','denominator','linear_K_forecast','center_K_float']}
 # Exact target duplicate check across accessible registrations/results.
 token=str(selected['candidate_numerators'][2]);matches=[]
 for root in [PROJECT/'upper',Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper'),Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper')]:
  for p in root.rglob('*.json'):
   if HERE in p.parents:continue
   if token in p.read_text():matches.append({'path':str(p),'sha256':sha(p)})
 assert not matches,'Possible duplicate target'
 selected.update({'point_index':1,'duplicate_scan_token':token,'duplicate_matches':matches,'registered_before_kernel':True})
 (HERE/'point1.json').write_text(json.dumps(selected,indent=2)+'\n')
else:selected=None
out={'status':'POINT1_REGISTERED'if selected else'STOP_NO_RELIABLE_FORECAST','rows':rows,'selected':selected,
 'model_is_not_a_K_certificate':True,'new_kernel_passes':0,'inputs':[{'path':str(p),'sha256':sha(p)}for p in [oldpath,HERE/'gradient_results.json',PROJECT/'upper/round009/short_potential/postprocess_results.json',HERE/'registration.md',Path(__file__)]]}
(HERE/'selection_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'selected':selected,'diagnostics':[{'T':r['T'],'cos':r['cosines_float_diagnostic'],'rel_diff':r['relative_max_differences_float_diagnostic'],'r':r['line_step'],'forecast':r['linear_K_forecast']}for r in rows]},indent=2))
