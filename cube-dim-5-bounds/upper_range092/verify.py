#!/usr/bin/env python3
"""Independent rescaled integer square roots, rational ranges and fourth powers."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,ceil
import gzip,json,hashlib,multiprocessing as mp,os,time,sys
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_range092');N=1<<160;ALT=None;TS=None

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def sqrt_scaled(q,bits):
 D=1<<bits;n=isqrt(q.numerator*D*D//q.denominator);lo=F(n,D);hi=F(n+1,D);assert lo*lo<=q<hi*hi;return lo,hi

def qualify_catalog():
 cat=json.loads((H/'ratio_catalog.json').read_text());alt={};proof=[]
 for xs,r in cat.items():
  x=F(xs);a,b,c,d=map(F,[r['sqrt1_lower'],r['sqrt1_upper'],r['sqrt2_lower'],r['sqrt2_upper']]);assert a*a<=1+x*x/16<b*b and c*c<=1+x*x/256<d*d;assert all((v*(1<<192)).denominator==1 for v in (a,b,c,d));assert b-a==d-c==F(1,1<<192);pL=(a+x/4)**2*(c+x/16)**8;pU=(b+x/4)**2*(d+x/16)**8;assert F(r['lower'])<=pL<pU<F(r['upper'])
  # Different normalization: roots of16+x² and256+x², then divide by2^36.
  al,au=sqrt_scaled(16+x*x,224);bl,bu=sqrt_scaled(256+x*x,224);lower=(al+x)**2*(bl+x)**8/(1<<36);upper=(au+x)**2*(bu+x)**8/(1<<36);U=F(ceil(upper*(1<<224)),1<<224);assert F(r['lower'])<=lower<=upper<=U<=F(r['upper']);alt[xs]=U;proof.append({'x':xs,'independent_upper':str(U),'source_upper':r['upper'],'sqrt_bits':224})
 (H/'ratio_verification.json').write_text(json.dumps({'status':'pass','count':len(cat),'all_primary192_squares_and_powers_outward':True,'independent_rescaled224_roots':True,'old_exp_asinh_H_calls':0,'cases':proof},indent=2)+'\n');return alt

def fourth_root_interval(q,bits=224):
 D=1<<bits;scaled=q.numerator*D**4//q.denominator;n=isqrt(isqrt(scaled));lo=F(n,D);assert lo**4<=q;hi=lo if lo**4==q else F(n+1,D);assert q<=hi**4;return lo,hi

def init(alt,ts):
 global ALT,TS
 ALT=alt;TS=ts

def one(row):
 rank=row['rank'];p=S/'nodes'/f'{rank//1024:02d}'/f'{rank:05d}.json.gz';r=json.loads(gzip.decompress(p.read_bytes()));assert r['status']=='complete' and r['key']==row['key'];assert sha(row['packet_path'])==row['packet_sha256'];packet=json.loads(gzip.decompress(Path(row['packet_path']).read_bytes()));assert packet['original_raw_sha256']==row['original_raw_sha256'];point=list(map(F,packet['point']));L=row['level'];offsets=[(F(i,2**L)-x,F(i+1,2**L)-x) for i,x in zip(row['indices'],point[1:])];h=max(abs(x) for span in offsets for x in span);mu=max(abs((a+b)/2) for a,b in offsets);assert h==F(row['h_infinity']) and mu==F(row['mu_infinity']);m2=F(row['M2']);ends=packet['endpoints'];assert [e['t'] for e in ends]==TS;low_grid=0;upper=F(0);mean=F(0);terms=[];D=1<<288
 for i,e in enumerate(ends):
  a=F(e['t']);nxt=ends[i+1]['H'] if i+1<73 else {'lower':'0','upper':'0'};ml=F(e['H']['lower'])-F(nxt['upper']);mh=F(e['H']['upper'])-F(nxt['lower']);assert 0<ml<=mh;R=ALT[str(h/a)];term=ml/R;low_grid+=term.numerator*D//term.denominator;terms.append(term);upper+=mh*R;mean+=mh*min(R,1+mu/a+F(25,48)*m2*R/a**2)
 l,u,m=F(r['rangeL_i'],N),F(r['rangeU_i'],N),F(r['mean_i'],N);lf=max(F(row['old_rangeL_i'],N),F(low_grid,D));fallback=False
 if l>lf:lf=max(F(row['old_rangeL_i'],N),sum(terms,F(0)));fallback=True
 assert l<=lf and l>=F(row['old_rangeL_i'],N);assert min(F(row['old_rangeU_i'],N),upper,3)<=u<=F(row['old_rangeU_i'],N);assert min(F(row['old_M_i'],N),mean,u)<=m<=F(row['old_M_i'],N);assert 0<l<=m<=u;vol=F(1,2**(4*L));il=F(r['new_power_lower_i'],N);assert (il/vol)**4<=l**5
 plL,plU=fourth_root_interval(l**5);puL,puU=fourth_root_interval(u**5);chordU=vol*plU if l==u else vol*((u-m)*plU+(m-l)*puU)/(u-l);assert chordU<=F(r['new_chord_upper_i'],N);assert r['new_self_L']==max(row['old_self_L'],row['weight']*r['new_power_lower_i']);assert r['new_self_U']==min(row['old_self_U'],row['weight']*r['new_chord_upper_i']);assert r['new_self_L']<=r['new_self_U'];assert r['new_H_calls']==r['old_exp_asinh_calls']==0
 return {'rank':rank,'key':row['key'],'record_sha256':hashlib.sha256(serial(r).encode()).hexdigest(),'pid':os.getpid(),'exact_lower_sum_fallback':fallback}

def hierarchy(proofs,reg):
 state=json.loads(gzip.decompress((S/'initial_state.json.gz').read_bytes()));nodes=state['nodes'];front=state['frontier'];proof={p['rank']:p for p in proofs};previous=None;prevL=F(reg['initial_LH']);prevU=F(reg['initial_UH']);count=0
 for bid in range(14):
  tag=f'{bid:04d}';brp=H/('block'+tag+'_registration.json');br=json.loads(brp.read_text());rp=H/('block'+tag+'_receipt.json');rc=json.loads(rp.read_text());cp=H/('block'+tag+'_cells.json.gz');assert sha(cp)==rc['certificates_sha256'] and sha(brp)==rc['registration_sha256'];assert br['parent_receipt_sha256']==previous;rows=json.loads(gzip.decompress(cp.read_bytes()));assert [r['rank'] for r in rows]==list(range(br['start'],br['end']));assert len(rows)==br['count']<=4096
  for r in rows:
   assert hashlib.sha256(serial(r).encode()).hexdigest()==proof[r['rank']]['record_sha256'];n=nodes[r['key']];assert n['selfL']==r['old_self_L'] and n['selfU']==r['old_self_U'];n['selfL'],n['selfU']=r['new_self_L'],r['new_self_U']
  def visit(k):
   n=nodes[k]
   if n['children']:
    for c in n['children']:visit(c)
    n['effL']=max(n['selfL'],sum(nodes[c]['effL'] for c in n['children']));n['effU']=min(n['selfU'],sum(nodes[c]['effU'] for c in n['children']))
   else:n['effL'],n['effU']=n['selfL'],n['selfU']
   assert n['effL']<=n['effU']
  for k in state['roots']:visit(k)
  lo=F(sum(nodes[k]['effL'] for k in state['roots']),N);up=F(sum(nodes[k]['effU'] for k in state['roots']),N);assert lo>=prevL and up<=prevU;assert lo==F(rc['after']['LH_truncated_integral']) and up==F(rc['after']['UH_truncated_integral']);bu=F(rc['after']['B_upper']);E=F(reg['E076_upper']);iu=F(rc['after']['I_upper']);cu=F(rc['after']['cost_upper']);g=F(reg['valid_g']);assert bu>E and (bu-E)**5>=up**4;assert iu**4>=bu**5 and cu**4>=(5*bu/g)**5;assert (5**5*bu**5<32**4*g**5)==rc['after']['cost_strictly_below32'];assert (iu<F(reg['old_complete_I_cap']))==rc['after']['I_strictly_below_old_complete_cap'];previous=sha(rp);prevL,prevU=lo,up;count+=len(rows)
 actual=json.loads(gzip.decompress((S/'state.json.gz').read_bytes()));assert actual['nodes']==nodes and actual['frontier']==front;ml=max(nodes[k]['level'] for k in front);assert sum(nodes[k]['weight']*2**(4*(ml-nodes[k]['level'])) for k in front)==2**(4*ml);return lo,up

def main():
 began=time.monotonic();reg=json.loads((H/'registration.json').read_text());assert sha(H/'engine.py')==reg['engine_sha256'] and sha(H/'ratio091.py')==reg['ratio091_sha256'];assert sha(H/'ratio_catalog.json')==reg['catalog_sha256'];rows=json.loads(gzip.decompress((H/'nodes.json.gz').read_bytes()));assert len(rows)==56139;alt=qualify_catalog();pool=mp.get_context('spawn').Pool(12,initializer=init,initargs=(alt,reg['time_endpoints']));proofs=[]
 try:
  for p in pool.imap_unordered(one,rows):
   proofs.append(p)
   if len(proofs)%8192==0:print('VERIFIED',len(proofs),round(time.monotonic()-began,3),flush=True)
  pool.close();pool.join()
 except BaseException:pool.terminate();pool.join();raise
 lo,up=hierarchy(proofs,reg);result={'status':'pass','nodes':len(proofs),'catalog_cases':len(alt),'all_primary192_root_squares_and_powers_checked':True,'independent_rescaled224_ratio_roots':True,'independent224_quarter_roots_all_chords':True,'all_range_mean_and_aggregate_intersections_checked':True,'all_new_lower_powers_exactly_checked':True,'global_B_I_cost_root_powers_exactly_checked':True,'exact_lower_sum_fallbacks':sum(p['exact_lower_sum_fallback'] for p in proofs),'all14_hierarchy_prefixes_valid':True,'LH':str(lo),'UH':str(up),'physical_frontier':55819,'weighted_volume':'1','new_H_calls':0,'old_exp_asinh_calls':0,'all12_workers_joined':True,'worker_pids':sorted({p['pid'] for p in proofs}),'pid':os.getpid(),'seconds':time.monotonic()-began,'source_sha256':sha(__file__)};(H/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
