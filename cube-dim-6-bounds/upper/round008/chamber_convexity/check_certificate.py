#!/usr/bin/env python3
"""Independent exact geometry and rational Jensen-sign reconstruction."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
from flint import arb,ctx,fmpq,fmpz
from certify_counterexample import exact_sign_dd
HERE=Path(__file__).resolve().parent;ctx.prec=384
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def triangular(nodes):
 z=sorted(nodes);b=[fmpq(1 if x>0 else-1)for x in z]
 for k in range(1,6):b=[(b[j-k+1]-b[j-k])/(z[j]-z[j-k])for j in range(k,6)]
 return b[0]
def main():
 data=json.loads((HERE/'exact_certificate.json').read_text());registry=json.loads((HERE/'probe_registry.json').read_text())
 triple=json.loads((HERE/'negative_triple.json').read_text());assert data['status']=='PASS_STRICT_CHAMBER_CONVEXITY_COUNTEREXAMPLE'
 assert registry['new_points_computed']==5 and len(registry['tests'])==2
 anchors=[]
 for raw in [[-3,-2,-1,1,2,4],[-4,-3,1,2,3,5],[1,2,3,4,5,6]]:
  z=list(map(fmpq,raw));a=exact_sign_dd(z);b=triangular(z);assert a==b
  anchors.append({'nodes':raw,'exact_DD5':str(a),'barycentric_equals_triangular':True})
 pts=[];rs=[]
 for row,probe in zip(data['points'],triple['points']):
  assert row['point_numerators']==probe['numerators']and row['denominator']==probe['denominator']
  x=[Q(n,row['denominator'])for n in row['point_numerators']]
  assert x[0]==1 and all(x[j]>x[j+1]for j in range(5))and x[-1]>0
  assert row['complete_group_terms']==967680 and row['all_nodes_nonzero_and_distinct']
  r=fmpq(fmpz(row['G_times_pi5_numerator']),fmpz(row['G_times_pi5_denominator']));assert r>0
  value=arb(r.numerator)/arb(r.denominator)/arb.pi()**5
  assert(value-arb(row['G_arb384'])).contains(0)
  pts.append(x);rs.append(r)
 assert all(pts[0][j]+pts[2][j]==2*pts[1][j]for j in range(6))
 direction=[Q(0),Q(0),Q(0),Q(1),Q(-1),Q(0)]
 assert all(pts[2][j]-pts[1][j]==Q(1,100)*direction[j]for j in range(6))
 gap=(rs[0]+rs[2])/2-rs[1]
 assert gap==fmpq(fmpz(data['exact_Jensen_gap_times_pi5_numerator']),fmpz(data['exact_Jensen_gap_times_pi5_denominator']))and gap<0
 interval=arb(gap.numerator)/arb(gap.denominator)/arb.pi()**5
 assert interval<-arb(13)/4000000000
 assert interval>-arb(163)/50000000000
 out={'status':'PASS_EXACT_CHAMBER_COUNTEREXAMPLE_RECONSTRUCTION','full_kernel_replayed_by_checker':False,
  'independent_kernel_arithmetic':'3 actual full-group fmpq barycentric evaluations in exact_certificate.json',
  'new_G_points':5,'interior_centers':1,'directions_tested':2,
  'strict_chamber_verified':True,'whole_segment_in_strict_convex_chamber':True,'exact_midpoint_verified':True,
  'h':'1/100','direction':[0,0,0,1,-1,0],'Jensen_gap_exact_rational_factor_negative':True,
  'Jensen_gap_arb384':interval.str(100,more=True),
  'Jensen_gap_strict_rational_interval':['-163/50000000000','-13/4000000000'],
  'convexity_on_ordered_chamber':False,'new_C6_bound':False,
  'scalar_anchors':anchors,'certificate_sha256':sha(HERE/'exact_certificate.json'),
  'registry_sha256':sha(HERE/'probe_registry.json'),'checker_source_sha256':sha(Path(__file__))}
 (HERE/'checked_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
