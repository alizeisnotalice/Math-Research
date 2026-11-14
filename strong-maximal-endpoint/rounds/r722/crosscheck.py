"""Reproduce high-precision cross-software and within-software summaries."""
from pathlib import Path
import json,re
import mpmath as mp
mp.mp.dps=120
HERE=Path(__file__).resolve().parent
def wl(s):return mp.mpf(re.sub(r'`[^* ]*','',s).replace('*^','e'))
def pyvals(z):return [mp.mpf(z['baseline']['c']),mp.mpf(z['baseline']['r']),mp.mpf(z['C_fixed']),mp.mpf(z['C_singleaxis']),mp.mpf(z['C_all']),mp.mpf(z['H_eigenvalues']['transverse_multiplicity_n_minus_1']),mp.mpf(z['H_eigenvalues']['longitudinal_multiplicity_1'])]
def wlvals(z):return [wl(z[k]) for k in ['baseline_c','baseline_r','fixed_curvature','single_axis_curvature','all_axis_curvature','H_transverse','H_parallel']]
if __name__=='__main__':
 p=json.loads((HERE/'curvature.json').read_text());w=json.loads((HERE/'wolfram.json').read_text());rows=[];cross=[];pdiff=[];wdiff=[];names=['c','r','C_fixed','C_singleaxis','C_all','H_transverse','H_parallel']
 for z in w['results']:
  if z.get('status')!='ok':raise ValueError(z)
  key=(z['n'],z['pattern']);hi=next(r for r in p['results'] if (r['n'],r['family'])==key and r['amplitude']=='1' and r['precision_dps']==70)
  delta=[abs(a-b) for a,b in zip(wlvals(z),pyvals(hi))];rows.append({'n':key[0],'pattern':key[1],'WL_precision':z['precision'],'WL_minus_Python70_abs':dict(zip(names,map(lambda t:mp.nstr(t,30),delta)))})
  if z['precision']==80:cross.extend(delta)
 for z in p['results']:
  if z['precision_dps']!=70:continue
  lo=next(r for r in p['results'] if (r['n'],r['family'],r['amplitude'],r['precision_dps'])==(z['n'],z['family'],z['amplitude'],40));pdiff.extend(abs(a-b) for a,b in zip(pyvals(z),pyvals(lo)))
 for z in w['results']:
  if z['precision']!=80:continue
  lo=next(r for r in w['results'] if (r['n'],r['pattern'],r['precision'])==(z['n'],z['pattern'],50));wdiff.extend(abs(a-b) for a,b in zip(wlvals(z),wlvals(lo)))
 out={'status':'matched','matched_WL_rows':len(rows),'independent_high_precision_structures':len(cross)//len(names),'max_WL80_Python70_abs_difference':mp.nstr(max(cross),40),'max_Python70_40_abs_difference_all72cases':mp.nstr(max(pdiff),40),'max_WL80_50_abs_difference':mp.nstr(max(wdiff),40),'scope':'Numeric comparison, not interval certification','comparisons':rows}
 (HERE/'crosscheck_summary.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='comparisons'})
