"""Reduce full local r094 batches to the minimal evidence used in the report."""
import argparse,json,hashlib,math
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--laplace',default='/tmp/laplace94');ap.add_argument('--normal',default='/tmp/normal_probe94');ap.add_argument('--wolfram',default='/tmp/wolfram94');ap.add_argument('--tilted',default='/tmp/r094-tilted-replay');a=ap.parse_args();out=Path(a.output);lp=Path(a.laplace);nr=Path(a.normal);wf=Path(a.wolfram);tf=Path(a.tilted)
 full=load(wf/'verification94.json');assert full['pass_all'];selected=[r for r in full['curves'] if r['j']%64==0];assert len(selected)==67
 wmin=min(full['curves'],key=lambda r:float(r['rows'][1]['curvature_ratio'].split('`')[0]));wsmall={k:v for k,v in full.items() if k not in ['curves']};wsmall.update(curves=selected,full_curve_count=len(full['curves']),all_full_curves_pass=all(r['pass'] for r in full['curves']),minimum_curvature_ratio_record=wmin,raw_full_output_sha256=sha(wf/'verification94.json'),selection='j divisible by 64; 67 of all 4225 curve records; all 24 tilt records retained');dump(out/'wolfram/verification94.json',wsmall)
 lrows=[json.loads(s) for s in (lp/'rows.jsonl').open()];assert len(lrows)==4225;ls=load(lp/'summary.json');ls.update(raw_rows_sha256=sha(lp/'rows.jsonl'),selected_rows=[r for r in lrows if r['dyadic_exponent_numerator']%64==0]);dump(out/'python/curvature94.json',ls)
 lo=hi=None;neg=0;ci_positive=0;by_n={};seeds=set();count=0
 for line in (nr/'rows.jsonl').open():
  x=json.loads(line);count+=1;seeds.add(x['seed']);by_n[x['N']]=by_n.get(x['N'],0)+1;neg+=x['gain_vs_reference']<0;ci_positive+=x['DKW_gain_interval'][0]>0
  if lo is None or x['gain_vs_reference']<lo['gain_vs_reference']:lo=x
  if hi is None or x['gain_vs_reference']>hi['gain_vs_reference']:hi=x
 assert count==30848 and len(seeds)==count
 ns=load(nr/'summary.json');ns.update(raw_rows_sha256=sha(nr/'rows.jsonl'),references=load(nr/'references.json'),cases_by_N=by_n,negative_observed_gain_count=neg,positive_DKW_lower_bound_count=ci_positive,minimum_gain_record=lo,maximum_gain_record=hi,total_paired_vectors=count*2000000,unique_seeds=len(seeds),method='PCG64; SHA256(repr(spec))[:8]; Gamma/common-thinned-Binomial CRN; empirical quantile with two-CDF DKW and 2/S interpolation allowance; no certified optimizer or physical halo assertion');dump(out/'python/normal94.json',ns)
 rows=load(tf/'rows.json');ts=load(tf/'summary.json');selected=[r for r in rows if r['a'] in (1,8) and r['xi'] in (-4,0,4)];assert len(rows)==315
 ts.update(selected_rows=selected,selection='a in {1,8}, xi in {-4,0,4}; all dimensions',raw_full_rows_sha256=sha(tf/'rows.json'));dump(out/'python/tilted94.json',ts)
 resources={}
 for name,path in [('curvature',Path('/tmp/laplace94_resources.json')),('normal',Path('/tmp/normal_probe94_resources.json')),('tilted',Path('/tmp/r094-tilted-replay-resources.json')),('wolfram',Path('/tmp/verify94_resources.json'))]:
  z=load(path);assert z['exit_code']==0 and z['stop_reason'] is None;resources[name]=z
 dump(out/'resources94.json',resources)
 print(json.dumps({'output':str(out),'laplace':4225,'Wolfram_curves':len(full['curves']),'tilted':315,'normal':count,'negative_gains':neg,'positive_CI_lower':ci_positive,'source':'minimal evidence; full raw batches retained outside Git'}))
if __name__=='__main__':main()
