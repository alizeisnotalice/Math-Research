"""One rationalized Gauss-node Hermite minorant for only selected m=3."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,floor
import hashlib,json,sys,time
import sympy as sp
from flint import arb,ctx
P=Path(__file__).resolve().parent
def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator
def main():
    began=time.monotonic();ctx.prec=384
    screen=json.loads((P/'screen_verification.json').read_text());m=screen['selected_m_for_only_next_stage'];assert m==3
    raw=json.loads((P/'hermite_moments.json').read_text());assert raw['completed_W_power']==6
    assert raw['backend_source_sha256']=='bfdf502cac067d2a1d9465c96ac15076c2e897cee89f8121bfed40977fc3fbf0'
    screening_moments=json.loads((P/'corner_screen.json').read_text())['moments']
    assert list(map(F,raw['moments']['1']))==list(map(F,screening_moments['d']))
    mu=[]
    for k in range(1,5):
        vals=list(map(F,raw['moments'][str(k)]))
        mu.append(sum((F(comb(6,n))*2**(6-n)*vals[n] for n in range(7)),F(0))/12**6)
    row=next(r for r in screen['rows'] if r['m']==m);assert mu[0]==F(row['E_d_u_2m'])
    V=F(row['E_B_u_3m']);variance=mu[0]*mu[2]-mu[1]**2;assert variance>0
    A=(mu[0]*mu[3]-mu[1]*mu[2])/variance;B=(mu[1]*mu[3]-mu[2]**2)/variance
    delta=A*A-4*B;assert delta>0
    exact_nodes=[(ar(A)-ar(delta).sqrt())/2,(ar(A)+ar(delta).sqrt())/2]
    z1,z2=exact_nodes
    weights=[(ar(mu[0])*z2-ar(mu[1]))/(z2-z1),(ar(mu[1])-ar(mu[0])*z1)/(z2-z1)]
    assert 0<z1<z2<1 and all(w>0 for w in weights)
    for j in range(4):
        quadrature=sum((w*z**j for w,z in zip(weights,exact_nodes)),arb(0))
        assert quadrature.contains(ar(mu[j]))
    ideal_U=sum((w/z.sqrt() for w,z in zip(weights,exact_nodes)),arb(0))
    ideal_cap=ideal_U**3/ar(V)**2/(32*arb.pi()**3)
    assert ideal_cap<ar(F(337,10**7))
    nodes=[]
    for z in exact_nodes:
        k=floor(float(z)*10000+.5);r=F(k,10000)
        assert ar(r-F(1,20000))<z<ar(r+F(1,20000));nodes.append(r)
    r1,r2=nodes;assert 0<r1<r2<1
    matrix=sp.Matrix([[1,r1,r1**2,r1**3],[0,1,2*r1,3*r1**2],[1,r2,r2**2,r2**3],[0,1,2*r2,3*r2**2]])
    inv=matrix.inv();combinations=[];h_brackets=[]
    for r in nodes:
        h=ar(1/r).sqrt();n=floor(float(h)*10**14);lo,hi=F(n-1,10**14),F(n+2,10**14)
        assert 0<lo and lo*lo<1/r<hi*hi;h_brackets.append((lo,hi))
    coeff=[]
    for j in range(4):
        q1=F(str(inv[j,0]))-F(str(inv[j,1]))/(2*r1)
        q2=F(str(inv[j,2]))-F(str(inv[j,3]))/(2*r2)
        combinations.append((q1,q2))
        coeff.append(sum((q*(ends[0] if q>=0 else ends[1]) for q,ends in zip((q1,q2),h_brackets)),F(0)))
    U=sum((a*b for a,b in zip(coeff,mu)),F(0));assert U>0
    ratio=U/V;scaled=ratio*10**5;a=F(scaled.numerator//scaled.denominator,10**5)
    assert a<=ratio<a+F(1,10**5)
    data014=json.loads((P.parent/'astra_upper014'/'budget.json').read_text())['cases'][0]
    cL=F(data014['c']['rational_lower']);target_upper=F(data014['g_required_for_cost_32_using_exact_expression']['rational_upper'])
    value=3*a*a*U-2*a**3*V;g=cL*value
    out={'status':'one_selected_m3_actual_Hermite_corner_candidate','new_global_g':False,'new_C5_bound':False,
         'm':m,'weighted_d_powers_1_to_4':[str(q) for q in mu],
         'all_seven_dW_moments_match_stage1':True,
         'gauss_monic_quadratic':['1',str(-A),str(B)],'gauss_nodes_Arb':[str(z) for z in exact_nodes],
         'gauss_weights_Arb':[str(z) for z in weights],
         'all_d_times_cubic_minorants_U_upper':str(ideal_U),
         'all_d_times_cubic_minorants_corner_ceiling':str(ideal_cap),
         'strict_simple_ceiling_for_all_d_times_cubic_minorants':'337/10000000',
         'rational_nodes':[str(r) for r in nodes],'inverse_sqrt_rational_brackets':[[str(q) for q in z] for z in h_brackets],
         'Hermite_coefficients_as_combinations':[[str(q) for q in z] for z in combinations],
         'lower_coefficients_for_d_powers_1_to_4':[str(q) for q in coeff],
         'U_lower_exact':str(U),'V_exact':str(V),'a':str(a),'c_lower':str(cL),
         'corner_polynomial_value':str(value),'corner_g_lower_exact':str(g),'corner_g_lower_Arb':str(ar(g)),
         'passes_014_required_g_necessary_only':g>target_upper,
         'source_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [P/'registration.md',P/'stage2_selection.md',P/'hermite_moments.cpp',P/'hermite_moments.json',P/'screen_verification.json',Path(__file__)]},
         'runtime':{'python':sys.executable,'arb_bits':ctx.prec,'seconds':time.monotonic()-began}}
    (P/'hermite_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['m','gauss_nodes_Arb','rational_nodes','U_lower_exact','V_exact','a','corner_g_lower_Arb','all_d_times_cubic_minorants_corner_ceiling','passes_014_required_g_necessary_only','runtime']},indent=2))
if __name__=='__main__':main()
