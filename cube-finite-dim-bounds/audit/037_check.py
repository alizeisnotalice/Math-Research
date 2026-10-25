"""Narrow fixed n4 coefficient sharpening; no author executable imports.
Sparse monomial expansion uses precomputed powers of linear forms, not compose.
035 exact whole-cone anchor and integral certificate are inherited, not rerun.
"""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import argparse,hashlib,json,time
from flint import fmpq_mpoly_ctx,fmpq,arb,ctx
HERE=Path(__file__).resolve().parent

def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def aq(v):
 v=F(v);return arb(v.numerator)/v.denominator

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--main',type=Path,default=HERE.parent);ap.add_argument('--save',type=Path,default=HERE/'037_results.json');args=ap.parse_args();root=args.main;start=time.perf_counter()
 source=load(HERE/'037_source_snapshot.json')
 for p,h in source.items():assert sha(root/p)==h,p
 C=fmpq_mpoly_ctx.get(('r','s','t','u'),ordering='deglex');v=C.gens();zero=C.constant(0);one=C.constant(1)
 cones=[[(1,0,0,0),(1,1,0,0),(2,1,1,0),(3,1,1,1)],[(2,1,1,0),(1,1,0,0),(1,1,1,0),(1,1,1,1)],[(2,1,1,0),(1,1,0,0),(1,1,1,1),(3,1,1,1)]]
 expected=[F(1,24),F(17,640),F(482191,13088544)];target=F(17,640);rows=[]
 for i,(row,rays) in enumerate(zip(load(root/'upper/u11_minimum_global_results.json')['physical_cases'],cones)):
  forms=[sum((rays[j][k]*v[j] for j in range(4)),zero) for k in range(4)]
  powers=[[one] for k in range(4)]
  for k in range(4):
   for degree in range(1,21):powers[k].append(powers[k][-1]*forms[k])
  def expand(d):
   ans=zero
   for key,value in d.items():
    e=tuple(map(int,key.split(',')));term=C.constant(fmpq(value))
    for k in range(4):term*=powers[k][e[k]]
    ans+=term
   return ans
  A=expand(row['numerator_dict'])*powers[0][3];D=expand(row['denominator_dict'])
  ad={tuple(map(int,e)):F(str(c)) for e,c in A.to_dict().items()};dd={tuple(map(int,e)):F(str(c)) for e,c in D.to_dict().items()};support=set(ad)|set(dd)
  assert dd and all(c>0 for c in dd.values()) and all(c>0 for c in ad.values())
  assert not(set(dd)-set(ad))
  best=min(ad[e]/c for e,c in dd.items());assert best==expected[i]
  extra=set(ad)-set(dd);assert all(ad[e]>0 for e in extra)
  diff={e:ad.get(e,F(0))-target*dd.get(e,F(0)) for e in support};assert min(diff.values())>=0
  assert all(ad.get(e,F(0))-best*dd.get(e,F(0))>=0 for e in support)
  rows.append({'case':row['case'],'A_terms':len(ad),'D_terms':len(dd),'A_only_terms':len(extra),'A_only_minimum':str(min((ad[e] for e in extra),default=F(0))),'minimum_ratio':str(best),'all_attaining_exponents':[list(e) for e in sorted(dd) if ad[e]/dd[e]==best],'uniform_difference_all_coefficients_nonnegative':True})
 assert [r['A_terms'] for r in rows]==[1739,1702,1702] and [r['D_terms'] for r in rows]==[1739,1683,1683]
 # The full1152cone rational residue calculation was independently done in035.
 representation=load(root/'audit/035_representation_results.json')
 anchors=representation['degenerate_fullsum_checks']
 hit=[r for r in anchors if r['x']==[1,1,1,0]];assert len(hit)==1
 assert hit[0]['zero_poles_cancel_exactly'] and F(hit[0]['G_times_pi3'])==target
 integral=load(root/'audit/035_integral_results.json');I0=F(92072351,500000000000)
 assert integral['cells']==2440 and F(integral['I4_upper_exact'])<=I0==F(integral['I_outward_rational'])
 # Branch-fixed Machin identity, exact finite alternating bounds.
 a=F(1,5);b=F(1,239);t2=2*a/(1-a*a);t4=2*t2/(1-t2*t2)
 assert (t4-b)/(1+t4*b)==1
 def alt(x,n):return sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n)),F(0))
 pm=16*alt(a,13)-4*alt(b,4);pup=F(355,113);assert pm<pup
 U0=F(14371,1000);power=256*I0**3*(F(640,17)*pup**3)**4;gap=U0**3-power;assert gap>0 and U0<F(14753,1000)
 ctx.prec=512
 true_cost=(256*aq(I0)**3*(aq(F(640,17))*arb.pi()**3)**4).root(3)
 result={'status':'PASS fixed n4 exact model minimum 17/(640*pi^3); C4 < 14371/1000 = 14.371',
 'source_sha256':source,'cases':rows,'global_minimum_coefficient':str(target),
 'model_minimum':'17/(640*pi^3)','anchor_inherited_035':hit[0],
 'I035_exact_upper':integral['I4_upper_exact'],'I_outward':str(I0),'pi_rational_upper':str(pup),'pi_machin_upper':str(pm),
 'pi_machin_rational_gap':str(pup-pm),'cost_cube_rational_upper':str(power),'strict_cube_gap':str(gap),'accepted_U4':str(U0),'cost_interval_512bit':str(true_cost),
 'new_Green_nodes':0,'new_integral_cells':0,'author_executable_calls':0,
 'shared_dependencies':'flint exact polynomial arithmetic and Arb; same locked global dictionaries already independently validated035; own035 anchor and full CI4B integral inherited; no author/Astra program executed',
 'wall_seconds':time.perf_counter()-start}
 args.save.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','cases','cost_interval_512bit','strict_cube_gap','wall_seconds']},indent=2))
if __name__=='__main__':main()
