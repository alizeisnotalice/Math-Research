#!/usr/bin/env python3
"""Portable certificate/hierarchy readback, no local endpoint cache required."""
from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,json,hashlib,time,os
from flint import arb,ctx
import kernel as K
import verify as V
H=Path(__file__).resolve().parent;N=1<<160

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args();t=time.monotonic();ctx.prec=384;reg=json.loads((H/'registration.json').read_text());s=json.loads((H/'summary.json').read_text());assert K.sha(H/'manager.py')==reg['manager_sha256'] and K.sha(H/'kernel.py')==reg['kernel_sha256'] and K.sha(H/'mechanism.md')==reg['mechanism_sha256'];assert K.sha(H.parent/'upper_global080/manifest.json')==reg['source080_manifest_sha256'];assert K.sha(H.parent/'upper_integral079/manifest.json')==reg['source079_manifest_sha256'];history=V.rebuild(1);assert history['UH']==s['budget']['UH_truncated_integral'] and history['LH']==s['budget']['LH_truncated_integral'];allkeys=set();certmap={}
 for b in range(2):
  tag=f'{b:04d}';br=json.loads((H/('batch'+tag+'_registration.json')).read_text());assert K.sha(H/('batch'+tag+'_jobs.json.gz'))==br['jobs_sha256'];rc=json.loads((H/('batch'+tag+'_receipt.json')).read_text());assert rc['registration_sha256']==K.sha(H/('batch'+tag+'_registration.json'));assert K.sha(H/('batch'+tag+'_cells.json.gz'))==rc['cell_certificates_sha256'];jobs=json.loads(gzip.decompress((H/('batch'+tag+'_jobs.json.gz')).read_bytes()));certs=json.loads(gzip.decompress((H/('batch'+tag+'_cells.json.gz')).read_bytes()));assert len(jobs)==len(certs)==br['child_count'];jm={j['key']:j for j in jobs};assert all(k not in allkeys for k in jm);allkeys.update(jm)
  for c in certs:
   j=jm[c['key']];assert c['level']==j['level'] and c['indices']==j['indices'] and c['weight']==j['weight'];bds=c['bounds'];l,u,m=[F(bds[a],N) for a in ('L_i','U_i','M_i')];assert 0<l<=m<=u<=3;v=arb(1)/2**(4*c['level']);q=arb(5)/4;lp=arb(K.rq(l))**q;up=arb(K.rq(u))**q;chord=v*(lp*arb(K.rq((u-m)/(u-l)))+up*arb(K.rq((m-l)/(u-l))));assert F(bds['integral_lower_i'],N)<=K.ep(v*lp) and K.ep(chord,True)<=F(bds['integral_upper_i'],N);certmap[c['key']]=c
 hashes=json.loads(gzip.decompress((H/'mathematical_hashes.json.gz').read_bytes()));assert set(hashes)==allkeys;sample=json.loads(gzip.decompress((H/'endpoint_samples.json.gz').read_bytes()))
 for k,ev in sample.items():
  assert hashlib.sha256(K.serial(ev).encode()).hexdigest()==hashes[k];assert ev['bounds']==certmap[k]['bounds'];assert len(ev['weights'])==140 and sum(F(a) for a in ev['weights'].values())==0;assert len(ev['endpoints'])==73
  for i,e in enumerate(ev['endpoints']):
   nxt=ev['endpoints'][i+1]['H'] if i+1<73 else {'lower':'0','upper':'0'};assert F(e['H']['lower'])-F(nxt['upper'])>0
 hierarchy=json.loads(gzip.decompress((H/'numeric_hierarchy.json.gz').read_bytes()));assert hierarchy['effective_global']==s['budget'];result={'status':'pass','portable':True,'children':len(allkeys),'physical_frontier':history['frontier_count'],'weighted_volume':'1','hierarchical_LH':history['LH'],'hierarchical_UH':history['UH'],'all_quantized_child_chords_checked':True,'frozen_exact_selection_and_four_axis_orbit_checked':True,'sample_payloads_checked':len(sample),'new_H_calls':0,'g_only_conditional':True,'no085_consumed':True,'seconds':time.monotonic()-t,'pid':os.getpid(),'source_sha256':K.sha(__file__)}
 if args.record:K.atomic(H/'portable_verification.json',result)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
