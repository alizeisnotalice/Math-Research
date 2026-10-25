"""Two fixed m necessary ceilings; no full-space conclusion."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,floor
import hashlib,json,time,sys
from flint import arb,ctx
P=Path(__file__).resolve().parent
def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator
def bracket(z,scale=10**15):
    k=floor(float(z)*scale);lo,hi=F(k-1,scale),F(k+2,scale)
    assert ar(lo)<z<ar(hi)
    return {'ball':str(z),'lower':str(lo),'upper':str(hi)}
def main():
    began=time.monotonic();ctx.prec=384
    raw=json.loads((P/'corner_screen.json').read_text());assert raw['completed_W_power']==9 and raw['beta']=='1/2'
    mom={k:list(map(F,v)) for k,v in raw['moments'].items()}
    old=json.loads((P.parent/'astra_upper015'/'corner_screen.json').read_text())['moments']
    for name,array in old.items():assert mom[name][:len(array)]==list(map(F,array))
    target=json.loads((P.parent/'astra_upper014'/'budget.json').read_text())['cases'][0]['g_required_for_cost_32_using_exact_expression']
    tlo,thi=F(target['rational_lower']),F(target['rational_upper'])
    def weighted(name,n):return sum((F(comb(n,k))*2**(n-k)*mom[name][k] for k in range(n+1)),F(0))/12**n
    c=1/(32*arb.pi()**3);rows=[]
    for m in (2,3):
        E=weighted('one',2*m);A=weighted('d',2*m);V=weighted('b',3*m);assert min(E,A,V)>0
        cap=c*ar(E*A)**(arb(3)/2)/ar(V)**2
        excluded=cap<ar(tlo);unexcluded=cap>ar(thi);assert excluded or unexcluded
        rows.append({'m':m,'E_u_2m':str(E),'E_d_u_2m':str(A),'E_B_u_3m':str(V),'U_squared_upper':str(E*A),
                     'all_a_corner_ceiling':bracket(cap),'excludes_014_required_g':bool(excluded),'passes_only_necessary_screen':bool(unexcluded)})
    candidates=[r for r in rows if r['passes_only_necessary_screen']]
    selected=max(candidates,key=lambda r:(F(r['all_a_corner_ceiling']['lower']),-r['m']))['m'] if candidates else None
    out={'status':'two_registered_ridge_power_necessary_corner_screens','new_global_g':False,'new_C5_bound':False,'rows':rows,'selected_m_for_only_next_stage':selected,
         'all_015_low_degree_moments_exactly_match':True,'max_graph_degree':13,
         'source_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [P/'registration.md',P/'corner_screen.cpp',P/'corner_screen.json',Path(__file__),P.parent/'astra_upper015'/'graph_moments_beta_half.hpp']},
         'runtime':{'python':sys.executable,'arb_bits':ctx.prec,'seconds':time.monotonic()-began}}
    (P/'screen_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='source_sha256'},indent=2))
if __name__=='__main__':main()
