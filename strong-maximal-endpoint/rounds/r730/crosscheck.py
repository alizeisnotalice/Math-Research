from pathlib import Path
import json,re
import mpmath as mp
mp.mp.dps=100
here=Path(__file__).resolve().parent
p=json.loads((here/'rank_tv_bounds.json').read_text());w=json.loads((here/'rank_tv_wolfram.json').read_text())
def number(s):return mp.mpf(re.sub(r'`{1,2}[0-9.+-]*','',s).replace('*^','e'))
rows=[];stability=[]
for r in w['results']:
 pr=next(x for x in p['results'] if x['n']==r['n'] and x['quad_tolerance']==1e-11)
 pairs=[('conditional_survival','conditional_survival'),('unconditional_survival','unconditional_survival'),('central_probability','central_probability_cdf'),('tv_lower_expression','tv_lower_expression')]
 rows.append(dict(n=r['n'],digits=r['precision_digits'],differences={k:mp.nstr(abs(number(r[k])-mp.mpf(repr(pr[j]))),30) for k,j in pairs}))
for n in sorted(set(x['n'] for x in w['results'])):
 a=next(x for x in w['results'] if x['n']==n and x['precision_digits']==50);b=next(x for x in w['results'] if x['n']==n and x['precision_digits']==80)
 stability.append(dict(n=n,difference=mp.nstr(abs(number(a['tv_lower_expression'])-number(b['tv_lower_expression'])),30)))
out=dict(status='complete',matched_cases=len(rows),scope='Floating/high-precision consistency only, not directed interval arithmetic',python_wolfram=rows,wolfram_precision_doubling=stability,max_python_wolfram_difference=mp.nstr(max(mp.mpf(v) for r in rows for v in r['differences'].values()),30),max_cdf_positive_integral_difference=max(abs(r['cdf_integral_difference']) for r in p['results']))
(here/'crosscheck.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
