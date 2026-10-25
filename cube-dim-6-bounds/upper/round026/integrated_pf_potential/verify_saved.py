"""Pure read verifier: no PF sum, positive-time sample, native call, or writes."""
import sys
sys.dont_write_bytecode=True
import json
from fractions import Fraction as Q
from pathlib import Path
from pf_integral import P,PF,F,sha,bounds,B,arb,ctx,endpoints,table
def contains(d,r):
 l,u=bounds(d);a,b=bounds(r);return l<=a<=b<=u
def main():
 m=json.loads((P/'manifest.json').read_text())
 for f,h in m['whitelist_sha256'].items():assert sha(P/f)==h,('local hash',f)
 for f,h in m['external_readonly_sha256'].items():assert sha(Path(f))==h,('external hash',f)
 reg=json.loads((P/'registration.json').read_text());r=json.loads((P/'execution_record.json').read_text())
 source22=next(Path(f) for f in reg['source_sha256'] if f.endswith('round022/reduced_point_preflight/results.json'))
 pre=json.loads(source22.read_text())
 assert pre['point_numerators']==reg['point_numerators']==[50,50,50,45,27,0] and pre['denominator']==reg['denominator']==50
 assert pre['zero_node_logical_terms']==88320 and pre['repeated_node_logical_terms']==201408
 assert r['native_executions']==r['positive_T_values']==1 and r['new_G_full']==0 and r['native_exit_code']==0 and not r['native_stderr']
 assert r['command'][-9:]==['1','25','50','50','50','50','45','27','0']
 s=json.loads((P/'scalar_checks.json').read_text());assert s['formal_scalar_cases']==len(s['cases'])==8 and s['positive_time_values']==0
 for row in s['cases']:assert row.get('residual',row.get('limit'))=='0'
 structural=table();assert len(structural)==417 # exact roots/zero channels only, no time substitution
 br=json.loads((PF/'results.json').read_text());fr=json.loads((F/'results.json').read_text());assert br['Gamma_G']==fr['Gamma_G']
 pf=json.loads((P/'pf_result.json').read_text());n=json.loads((P/'native_H.json').read_text());c=json.loads((P/'comparison.json').read_text())
 assert n['full_group_terms']==967680 and n['computed_DD_evaluations']==483840 and n['zero_node_terms']==88320 and n['confluent_terms']==201408
 assert c['status']=='PASS' and c['PF640_strictly_inside_G_minus_H512']
 # Independent exact endpoint Minkowski subtraction, without Arb radius rewrapping.
 gl,gu=bounds(pf['G_from_frozen_Gamma']);hl,hu=bounds(n['dyadic_interval']);kl,ku=bounds(pf['K_T'])
 assert gl-hu<kl<=ku<gu-hl
 assert kl>0 and bounds(pf['p_T'])[0]>0
 ctx.prec=1024
 G=B(Q(br['Gamma_G']))/(128*arb.pi()**5)
 density=B(Q(pf['density_raw_rational']))/(64*arb.pi()**6)
 assert contains(pf['G_from_frozen_Gamma'],endpoints(G))
 assert contains(pf['p_T'],endpoints(density))
 assert pf['T']=='1/25' and pf['precision_bits']==640
 print('PASS pure read: hashes,8 formal receipts,359 rational square roots,exact saved G-H endpoint inclusion,p-density conversion;0 PF sums/0 kernels/0 writes')
if __name__=='__main__':main()
