"""Read-only: no executable launches, compilation, full-group sums, or writes."""
import sys
sys.dont_write_bytecode=True
import hashlib,json
from fractions import Fraction as Q
from pathlib import Path
from check_anchors import P,bounds,contained
from validate_full import check,HA,HT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 manifest=json.loads((P/'manifest.json').read_text())
 for name,h in manifest['whitelist_sha256'].items():assert sha(P/name)==h,('archive hash',name)
 for name,h in manifest['external_readonly_sha256'].items():assert sha(Path(name))==h,('external hash',name)
 reg=json.loads((P/'registration.json').read_text());n=json.loads((P/'anchors_native.json').read_text());v=json.loads((P/'anchors_validation.json').read_text())
 assert n['parser_selftests_passed'] and len(n['anchors'])==len(v['cases'])==len(reg['anchors'])==6
 assert v['native_sha256']==sha(P/'anchors_native.json')
 for a,r in zip(n['anchors'],v['cases']):
  assert contained(a['value'],r['band_DD5_reference'])
  assert all(contained(x,y) for x,y in zip(a['partials'],r['partials_DD6_reference']))
  for label,q in r['density_exact_rationals'].items():
   l,u=bounds(a[label]);assert l<=Q(q)<=u
 assert all(bounds(v)==(Q(0),Q(0)) for v in n['anchors'][0]['partials'])
 build=json.loads((P/'build_record.json').read_text());assert sha(P/'band_engine.cpp')==build['source_sha256']
 rec=json.loads((P/'execution_record.json').read_text());assert rec['attempt']==1 and rec['exit_code']==0 and not rec['stderr']
 assert rec['source_sha256']==build['source_sha256'] and rec['binary_sha256']==build['binary_sha256']
 ha=json.loads(HA.read_text());ht=json.loads(HT.read_text())
 a=next(r for r in ha['time_nodes'] if r['tau']=='1/250');t=next(r for r in ht['rows'] if r['label']=='A007_reference')
 actual=check(json.loads((P/'native_result.json').read_text()),a,t)
 saved=json.loads((P/'full_validation.json').read_text())
 for k,val in actual.items():assert saved[k]==val,('full validation',k)
 assert saved['full_gradient_crosschecked'] is False
 print('PASS saved-only: whitelist/external hashes,6 scalar cases,full band containment,positive endpoint densities,Fraction Euler;0 kernels/0 writes')
if __name__=='__main__':main()
