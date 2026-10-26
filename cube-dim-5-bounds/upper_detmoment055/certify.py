#!/usr/bin/env python3
"""One unweighted beta-half two-node Gauss lower and pure Minkowski gate."""
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx
import hashlib,json,os,time
HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def root_interval(q,k,digits=80):
 q=F(q);assert q>0;den=10**digits;n=q.numerator*den**k//q.denominator
 lo=0;hi=1<<((n.bit_length()+k-1)//k)
 while lo+1<hi:
  mid=(lo+hi)//2
  if mid**k<=n:lo=mid
  else:hi=mid
 low,high=F(lo,den),F(lo+1,den);assert low**k<=q<high**k
 return low,high
def ar(q):
 q=F(q);return arb(q.numerator)/q.denominator
def serial(pair):return [str(q) for q in pair]
def main():
 started=time.monotonic();ctx.prec=384;reg=json.loads((HERE/'registration.json').read_text())
 for path,record in reg['sources'].items():assert sha(Path(path))==record['sha256']
 raw=json.loads((MAIN/'astra_upper016/hermite_moments.json').read_text())
 assert raw['backend_source_sha256']=='bfdf502cac067d2a1d9465c96ac15076c2e897cee89f8121bfed40977fc3fbf0'
 assert raw['canonical_moments']==55942 and raw['completed_W_power']==6
 m=[F(1)]+[F(raw['moments'][str(k)][0]) for k in range(1,4)]
 variance=m[0]*m[2]-m[1]**2;assert variance>0
 A=(m[0]*m[3]-m[1]*m[2])/variance;B=(m[1]*m[3]-m[2]**2)/variance
 orth=[m[2]-A*m[1]+B*m[0],m[3]-A*m[2]+B*m[1]];assert orth==[0,0]
 disc=A*A-4*B;assert disc>0;sl,sh=root_interval(disc,2)
 nodes=[((A-sh)/2,(A-sl)/2),((A+sl)/2,(A+sh)/2)]
 (l1,h1),(l2,h2)=nodes
 p=lambda x:x*x-A*x+B
 assert 0<l1<h1<l2<h2<1
 assert p(l1)>0>p(h1) and p(l2)<0<p(h2)
 assert h1<m[1]<l2
 denlow,denhigh=l2-h1,h2-l1
 weights=[((l2-m[1])/denhigh,(h2-m[1])/denlow),((m[1]-h1)/denhigh,(m[1]-l1)/denlow)]
 assert all(0<lo<hi<1 for lo,hi in weights)
 f_bounds=[(root_interval(1/hi,10)[0],root_interval(1/lo,10)[1]) for lo,hi in nodes]
 jlo=sum(w[0]*f[0] for w,f in zip(weights,f_bounds));jhi=sum(w[1]*f[1] for w,f in zip(weights,f_bounds))
 # Independent full Arb balls from exact rational A,B, not decimal nodes.
 x1=(ar(A)-ar(disc).sqrt())/2;x2=(ar(A)+ar(disc).sqrt())/2
 ws=[(x2-ar(m[1]))/(x2-x1),(ar(m[1])-x1)/(x2-x1)]
 assert 0<x1<x2<1 and all(w>0 for w in ws)
 for k in range(4):assert sum((w*x**k for w,x in zip(ws,[x1,x2])),arb(0)).contains(ar(m[k]))
 J=sum((w*x**(-arb(1)/10) for w,x in zip(ws,[x1,x2])),arb(0))
 assert ar(jlo)<J<ar(jhi)
 budget=json.loads((MAIN/'round010/orbit_budget.json').read_text())
 pl,ph=map(F,budget['pi_machin80_interval']);assert 3<pl<ph<F(22,7)
 g0=F(1,10000)
 coefficient=F(60**4*32**9,3125)
 gate=(g0*root_interval(coefficient*pl**7,5)[0],g0*root_interval(coefficient*ph**7,5)[1])
 gate_arb=ar(g0)*(ar(coefficient)*arb.pi()**7)**(arb(1)/5)
 assert ar(gate[0])<gate_arb<ar(gate[1])
 separated=jlo>gate[1]
 opposite=jhi<gate[0]
 # Lower bound on the pure Minkowski expression's cost at g0, not actual cost.
 cost4lo=F(3125)*jlo**5/(60**4*32**5*ph**7*g0**5)
 costlo=root_interval(cost4lo,4)[0]
 if separated:assert costlo>32
 assert time.monotonic()-started<120
 out={'status':'PASS_pure_Minkowski_budget_obstruction' if separated else 'LOCAL_NO_CONCLUSION_two_node_lower_not_separating',
  'pid':os.getpid(),'unweighted_moments_0_to_3':[str(q) for q in m],'moments_were_existing_exact_results':True,'new_graph_states':0,'graph_backend_calls':0,
  'source_cache_read_only_hash_verified':reg['cache']['sha256'],'cache_hit_claimed':False,'cache_expanded':False,
  'variance':str(variance),'p2_coefficients':['1',str(-A),str(B)],'p2_orthogonality_residuals':serial(orth),'discriminant':str(disc),
  'exact_sqrt_discriminant_interval':serial((sl,sh)),'nodes_rational_intervals':[serial(q) for q in nodes],'weights_rational_intervals':[serial(q) for q in weights],
  'node_function_tenth_root_intervals':[serial(q) for q in f_bounds],'Gauss_value_rational_interval':serial((jlo,jhi)),
  'nodes_Arb':[str(x1),str(x2)],'weights_Arb':[str(w) for w in ws],'Gauss_value_Arb':str(J),
  'f_fourth_derivative':'7161/10000 * d^(-41/10) > 0; global two-double-node cubic Hermite minorant on (0,1]',
  'J_lower_meaning':'E_R d^(-1/10) >= exact Gauss value >= rational interval lower endpoint',
  'g0':str(g0),'cost_target':'32','necessary_pure_Minkowski_J_gate_formula':'g0*(60^4*32^9*pi^7/3125)^(1/5)',
  'J_gate_rational_interval':serial(gate),'J_gate_Arb':str(gate_arb),'J_lower_minus_gate_upper_exact':str(jlo-gate[1]),
  'lower_exceeds_gate_upper':separated,'Gauss_value_upper_below_gate_lower':opposite,
  'pure_Minkowski_cost_expression_lower':str(costlo),'pure_Minkowski_cost_lower_display_only':float(costlo),
  'context_Jflag_full_interval':budget['J_flag_full_source_interval'],'context_current_I_full_interval':budget['new_I_budget_expression_interval'],
  'source_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),'source_inputs':reg['sources'],
  'runtime':{'seconds':time.monotonic()-started,'arb_bits':ctx.prec,'root_digits':80,'new_moment_production':False},
  'scope':'Only pure Minkowski J budget at g0/cost32. No lower bound for actual I or C5, and no exclusion of all-orbit improved bounds.'}
 (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k] for k in ['status','pid','unweighted_moments_0_to_3','nodes_Arb','weights_Arb','Gauss_value_Arb','J_gate_Arb','lower_exceeds_gate_upper','pure_Minkowski_cost_lower_display_only','runtime']},indent=2))
if __name__=='__main__':main()
