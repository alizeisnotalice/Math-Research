from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time
import mpmath as mp
p=Path(__file__).resolve().parent;mp.mp.dps=120;H=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();s=json.loads((p/'frontier4_summary.json').read_text());best={}
for z in s['receipts']:
 j=json.loads((p/z['file']).read_text())
 for x in j['top32']:
  k=x['il']
  if k not in best or x['cost_float']>best[k]['cost_float']:best[k]=x
pool=s['top64_float_proposals']+list(best.values());unique={tuple(x[k] for k in ('il','iq','isigma','ichi','hden')):x for x in pool}
def v(x):
 x=F(x);return mp.mpf(x.numerator)/x.denominator
def S64(x):
 acc=t=mp.mpf(1)
 for j in range(1,65):t*=x/j;acc+=t
 return acc
rows=[]
for row in unique.values():
 L=v(row['L']);q=v(row['qmax']);sigma=v(row['sigma']);chi=v(row['chi']);D=mp.mpf(row['hden']);eJ=mp.power(2,-20);eB=mp.power(2,-19);eN=eB;eS=eJ;eT=eB;eD=eJ;ea=mp.mpf(1)/65536;eR=eJ;trim=mp.power(2,-19);rlo=1+trim;rhi=mp.mpf(3)/2*q-trim;qmin=mp.mpf(2)/3;x=(1-sigma)*L;eta=(1-q)/(2*q);gap=chi*x-mp.mpf(183)/1000-mp.mpf(2211)/50000*x/(1-eJ);yS=((1-mp.exp(-L))*eta-chi)*x;yT=(sigma/2-chi*(1-sigma))*L;cap=mp.log(q/qmin)*rhi*D*L/S64(yS);tau=rhi*mp.sqrt(L)*D*L/S64(yT);retain=gap**2/(mp.mpf(531)/100000*x+gap**2);P0=1-eB-6/D;ul=P0*(1-eN)/(1+eJ);uh=(1+eN)/(P0*(1-eJ));gamma=(1-eD)*(1-ea)*mp.log(rhi/rlo)*(2*retain*(1-eR)*ul-eS-eT)/(2*(1+eR)*uh);price=mp.log(1-sigma);remainder=mp.exp(-L)+mp.exp(-sigma*L)+2*mp.sqrt(L)*mp.exp(-sigma*L/2);cost=gamma+price-remainder;ok=all([gap>0,yS>0,yT>0,2/(D*L)<trim,cap<eS,tau<eT,cost>0]);rows.append({'parameters':{k:row[k] for k in ('L','qmax','sigma','chi','hden','il','iq','isigma','ichi')},'cost':mp.nstr(cost,80),'gamma':mp.nstr(gamma,80),'retention':mp.nstr(retain,80),'cap_prefix_upper':mp.nstr(cap,80),'tau_prefix_upper':mp.nstr(tau,80),'gap':mp.nstr(gap,80),'checked_selected_numeric_gates':ok})
rows.sort(key=lambda x:mp.mpf(x['cost']),reverse=True);out={'scope':'post-run 120-decimal candidate diagnostics only, no interval proof or Cn theorem; unspecified full finite guards unverified','source_sha256':H(Path(__file__)),'screen_summary_sha256':H(p/'frontier4_summary.json'),'precision_decimal':120,'prefix':'S64 positive partial exponential polynomial, reciprocal gives upper diagnostic','selection':'union global top64 and each L index best, deduplicated','rows':rows,'executed_unix':time.time()};f=p/'frontier4_recheck.json';f.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'selected':len(rows),'checked':sum(x['checked_selected_numeric_gates'] for x in rows),'best':rows[0],'receipt_sha256':H(f),'source_sha256':out['source_sha256']}))
