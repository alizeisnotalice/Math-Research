from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,json,time,os
from flint import arb,ctx
import engine as E
H=Path(__file__).resolve().parent;N=1<<160

def main():
 t=time.monotonic();ctx.prec=384;ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args();reg=json.loads((H/'registration.json').read_text());summary=json.loads((H/'summary.json').read_text());assert E.sha(H/'engine.py')==reg['engine_sha256'] and E.sha(H/'moment085.py')==reg['moment085_sha256'];assert E.sha(H/'nodes.json.gz')==reg['nodes_sha256'];assert E.sha(H/'globalg_acceptance012.json')==reg['g_acceptance_sha256'];
 for dirname,field in [('upper_global080','source080_manifest_sha256'),('upper_adapt082','source082_manifest_sha256'),('upper_integral079','source079_manifest_sha256')]:
  d=H.parent/dirname;assert E.sha(d/'manifest.json')==reg[field];mf=json.loads((d/'manifest.json').read_text());assert all(E.sha(d/k)==v['sha256'] for k,v in mf['files'].items())
 nodes0=json.loads(gzip.decompress((H/'nodes.json.gz').read_bytes()));assert len(nodes0)==56139;lookup={r['key']:dict(r,children=None,currentU=r['old_self_U']) for r in nodes0};graph=json.loads(gzip.decompress((H.parent/'upper_adapt082/numeric_hierarchy.json.gz').read_bytes()))
 for k,n in graph['expanded_nodes'].items():lookup[k]['children']=n['children']
 rootkeys=[r['key'] for r in nodes0 if r['level']==5];front=[r['key'] for r in nodes0 if r['physical_frontier']];assert len(rootkeys)==52360 and len(front)==55819;maxlevel=max(lookup[k]['level'] for k in front);assert sum(lookup[k]['weight']*2**(4*(maxlevel-lookup[k]['level'])) for k in front)==2**(4*maxlevel);previousU=F(reg['initial_UH']);previous_receipt=None
 for b in range(14):
  tag=f'{b:04d}';bp=H/('block'+tag+'_registration.json');br=json.loads(bp.read_text());rp=H/('block'+tag+'_receipt.json');rc=json.loads(rp.read_text());cp=H/('block'+tag+'_cells.json.gz');assert E.sha(cp)==rc['certificates_sha256'] and E.sha(bp)==rc['registration_sha256'];assert br['parent_receipt_sha256']==previous_receipt;selected=nodes0[br['start']:br['end']];assert __import__('hashlib').sha256(E.serial(selected).encode()).hexdigest()==br['selection_sha256'];rows=json.loads(gzip.decompress(cp.read_bytes()));assert [r['rank'] for r in rows]==list(range(br['start'],br['end']))
  for r in rows:
   row=lookup[r['key']];assert row['currentU']==r['old_self_U'] and row['old_self_L']==r['old_self_L'];l=F(r['rangeL_i'],N);u=F(r['rangeU_i'],N);m=F(r['new_mean_i'],N);assert l<=m<=F(r['old_mean_i'],N)<=u;v=arb(1)/2**(4*r['level']);q=arb(5)/4;lp=arb(E.rq(l))**q;up=arb(E.rq(u))**q;chord=v*(lp*arb(E.rq((u-m)/(u-l)))+up*arb(E.rq((m-l)/(u-l))));assert E.ep(chord,True)<=F(r['new_chord_integral_upper_i'],N);assert r['new_self_U']==min(r['old_self_U'],r['weight']*r['new_chord_integral_upper_i']);row['currentU']=r['new_self_U'];assert r['new_H_calls']==0
  for k in sorted(lookup,key=lambda k:lookup[k]['level'],reverse=True):
   n=lookup[k];n['effectiveU']=min(n['currentU'],sum(lookup[c]['effectiveU'] for c in n['children'])) if n['children'] else n['currentU'];n['effectiveL']=max(n['old_self_L'],sum(lookup[c]['effectiveL'] for c in n['children'])) if n['children'] else n['old_self_L'];assert n['effectiveL']<=n['effectiveU']
  up=F(sum(lookup[k]['effectiveU'] for k in rootkeys),N);lo=F(sum(lookup[k]['effectiveL'] for k in rootkeys),N);assert up<=previousU and up==F(rc['after']['UH_truncated_integral']);assert lo==F(reg['initial_LH']);B=E.ep(arb(E.rq(up))**(arb(4)/5)+arb(E.rq(reg['E076_upper'])),True);assert B<=F(rc['after']['B_upper']);g=F(1,10000);assert (5**5*F(rc['after']['B_upper'])**5<32**4*g**5)==rc['after']['cost_strictly_below32'];previousU=up;previous_receipt=E.sha(rp)
 assert up==F(summary['budget']['UH_truncated_integral']);result={'status':'pass','portable':True,'nodes':56139,'all_quantized_chords_checked':True,'all14_hierarchy_prefixes_checked':True,'physical_frontier':55819,'weighted_volume':'1','LH_unchanged':str(lo),'UH':str(up),'qualified_g':'1/10000','new_H_calls':0,'seconds':time.monotonic()-t,'pid':os.getpid(),'source_sha256':E.sha(__file__)}
 if args.record:E.atomic(H/'portable_verification.json',result)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
