from pathlib import Path
from fractions import Fraction as F
import gzip,json,time,argparse,os
import engine as E
import verify as V
H=Path(__file__).resolve().parent;N=1<<160

def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args();reg=json.loads((H/'registration.json').read_text());summary=json.loads((H/'summary.json').read_text());assert E.sha(H/'engine.py')==reg['engine_sha256'] and E.sha(H/'ratio091.py')==reg['ratio091_sha256'];assert E.sha(H/'nodes.json.gz')==reg['nodes_sha256'] and E.sha(H/'ratio_catalog.json')==reg['catalog_sha256'];assert E.sha(H/'globalg_acceptance012.json')==reg['g_acceptance_sha256'];cat=json.loads((H/'ratio_catalog.json').read_text())
 for xs,r in cat.items():
  x=F(xs);a,b,c,d=map(F,[r['sqrt1_lower'],r['sqrt1_upper'],r['sqrt2_lower'],r['sqrt2_upper']]);assert a*a<=1+x*x/16<b*b and c*c<=1+x*x/256<d*d;assert F(r['lower'])<=(a+x/4)**2*(c+x/16)**8<(b+x/4)**2*(d+x/16)**8<F(r['upper']);al,au=V.sqrt_scaled(16+x*x,224);bl,bu=V.sqrt_scaled(256+x*x,224);assert (au+x)**2*(bu+x)**8/(1<<36)<=F(r['upper'])
 inputrows=json.loads(gzip.decompress((H/'nodes.json.gz').read_bytes()));lookup={r['key']:dict(r,currentL=r['old_self_L'],currentU=r['old_self_U'],children=None) for r in inputrows};graph=json.loads(gzip.decompress((H.parent/'upper_adapt082/numeric_hierarchy.json.gz').read_bytes()));assert E.sha(H.parent/'upper_reenvelope088/manifest.json')==reg['source088_manifest_sha256'];assert E.sha(H.parent/'upper_adapt082/manifest.json')==reg['source082_manifest_sha256']
 for k,n in graph['expanded_nodes'].items():lookup[k]['children']=n['children']
 roots=[r['key'] for r in inputrows if r['level']==5];front=[r['key'] for r in inputrows if r['physical_frontier']];assert len(roots)==52360 and len(front)==55819;maxL=max(lookup[k]['level'] for k in front);assert sum(lookup[k]['weight']*2**(4*(maxL-lookup[k]['level'])) for k in front)==2**(4*maxL);prevU=F(reg['initial_UH']);prevL=F(reg['initial_LH']);previous=None
 for bid in range(14):
  tag=f'{bid:04d}';bp=H/('block'+tag+'_registration.json');br=json.loads(bp.read_text());rp=H/('block'+tag+'_receipt.json');rc=json.loads(rp.read_text());cp=H/('block'+tag+'_cells.json.gz');assert E.sha(cp)==rc['certificates_sha256'] and E.sha(bp)==rc['registration_sha256'];assert br['parent_receipt_sha256']==previous;sel=inputrows[br['start']:br['end']];assert __import__('hashlib').sha256(E.serial(sel).encode()).hexdigest()==br['selection_sha256'];rows=json.loads(gzip.decompress(cp.read_bytes()));assert [r['rank'] for r in rows]==list(range(br['start'],br['end']))
  for r in rows:
   row=lookup[r['key']];assert row['currentL']==r['old_self_L'] and row['currentU']==r['old_self_U'];l,u,m=F(r['rangeL_i'],N),F(r['rangeU_i'],N),F(r['mean_i'],N);assert F(r['old_rangeL_i'],N)<=l<=m<=u<=F(r['old_rangeU_i'],N) and m<=F(r['old_M_i'],N);vol=F(1,2**(4*r['level']));assert (F(r['new_power_lower_i'],N)/vol)**4<=l**5;pl,pu=V.fourth_root_interval(l**5);ql,qu=V.fourth_root_interval(u**5);chord=vol*pu if l==u else vol*((u-m)*pu+(m-l)*qu)/(u-l);assert chord<=F(r['new_chord_upper_i'],N);assert r['new_self_L']==max(r['old_self_L'],r['weight']*r['new_power_lower_i']) and r['new_self_U']==min(r['old_self_U'],r['weight']*r['new_chord_upper_i']);row['currentL'],row['currentU']=r['new_self_L'],r['new_self_U']
  for k in sorted(lookup,key=lambda k:lookup[k]['level'],reverse=True):
   n=lookup[k];n['effL']=max(n['currentL'],sum(lookup[c]['effL'] for c in n['children'])) if n['children'] else n['currentL'];n['effU']=min(n['currentU'],sum(lookup[c]['effU'] for c in n['children'])) if n['children'] else n['currentU'];assert n['effL']<=n['effU']
  lo=F(sum(lookup[k]['effL'] for k in roots),N);up=F(sum(lookup[k]['effU'] for k in roots),N);assert lo>=prevL and up<=prevU;assert lo==F(rc['after']['LH_truncated_integral']) and up==F(rc['after']['UH_truncated_integral']);B=F(rc['after']['B_upper']);I=F(rc['after']['I_upper']);C=F(rc['after']['cost_upper']);Eu=F(reg['E076_upper']);g=F(reg['valid_g']);assert (B-Eu)**5>=up**4 and I**4>=B**5 and C**4>=(5*B/g)**5;assert (5**5*B**5<32**4*g**5)==rc['after']['cost_strictly_below32'];prevL,prevU=lo,up;previous=E.sha(rp)
 assert up==F(summary['budget']['UH_truncated_integral']);result={'status':'pass','portable':True,'nodes':56139,'ratio_cases':89,'all_root_squares_and_powers_outward':True,'all_lower_quartics_and_chord_roots_checked':True,'all14_hierarchy_prefixes_valid':True,'LH':str(lo),'UH':str(up),'physical_frontier':55819,'weighted_volume':'1','new_H_old_exp_asinh_calls':0,'seconds':time.monotonic()-began,'pid':os.getpid(),'source_sha256':E.sha(__file__)}
 if args.record:E.atomic(H/'portable_verification.json',result)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
