from pathlib import Path
import json,re
import mpmath as mp
mp.mp.dps=100
here=Path(__file__).resolve().parent
p=json.loads((here/'cube_halo.json').read_text());w=json.loads((here/'cube_halo_wolfram.json').read_text())
def wl(s):return mp.mpf(re.sub(r'`{1,2}[0-9.+-]*','',s).replace('*^','e'))
rows=[]
for z in w['results']:
 r=next(t for t in p['results'] if t['n']==z['n'] and t['target_digits']==z['target_digits'])
 rows.append(dict(n=z['n'],digits=z['target_digits'],saved_decimal_log_difference=mp.nstr(abs(wl(z['log_halo'])-mp.mpf(r['log_halo'])),30)))
st=[]
for n in p['dimensions']:
 a=next(t for t in p['results'] if t['n']==n and t['target_digits']==40);b=next(t for t in p['results'] if t['n']==n and t['target_digits']==70)
 st.append(dict(n=n,saved_decimal_log_difference=mp.nstr(abs(mp.mpf(a['log_halo'])-mp.mpf(b['log_halo'])),30)))
out=dict(status='complete',scope='Cube only, floating/saved-decimal comparison not interval certification',matched_rows=len(rows),max_same_target_difference={str(k):mp.nstr(max(mp.mpf(t['saved_decimal_log_difference']) for t in rows if t['digits']==k),30) for k in [40,70]},rows=rows,python_precision_doubling=st,maximum_float_tail_relative=mp.nstr(max(mp.mpf(t['relative_analytic_tail_bound_float']) for t in p['results']),30),maximum_scipy_log_difference=max(abs(t.get('scipy_minus_mp_log',0)) for t in p['results']))
(here/'cube_crosscheck.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','matched_rows','max_same_target_difference','maximum_float_tail_relative','maximum_scipy_log_difference']}))
