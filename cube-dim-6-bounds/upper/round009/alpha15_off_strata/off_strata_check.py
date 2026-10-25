"""Stable rank2 deletion and four-distance-layer scalar budgets.

Reads frozen peak intervals. No old peak recomputation, MC, or spatial grid.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,itertools,platform,time
import sympy as sp
from flint import arb,ctx
ctx.prec=256
HERE=Path(__file__).resolve().parent
START=time.monotonic()
REG=json.loads((HERE/'registration.json').read_text())
SRC=HERE.parents[2]/'audit/round008/alpha15_product/check_results.json'
OLD=HERE.parents[1]/'round008'/'alpha15_coercivity'/'coercivity_results.json'
# HERE=upper/round009/alpha15_off_strata, parents[1]=upper.
S=json.loads(SRC.read_text());O=json.loads(OLD.read_text())


def aq(x):
    x=F(x);return arb(x.numerator)/x.denominator


def less(x,y):
    if x is y:return False
    if x<y:return True
    if y<x:return False
    raise ArithmeticError('unresolved interval comparison')


def mn(x,y):return x if less(x,y) else y
def mx(x,y):return y if less(x,y) else x


def power_integral(a,e,l,u):
    if l is u:return arb(0)
    assert less(l,u)
    exponent=e+1
    if exponent==0:return a*(u/l).log()
    if l.is_zero():
        assert exponent>0
        return a*u**aq(exponent)/aq(exponent)
    return a*(u**aq(exponent)-l**aq(exponent))/aq(exponent)


def density_pieces(peak,one,two):
    z=arb(0)
    if two is None:return [(z,(peak/one)**(arb(3)/5),one,F(-7,3)),((peak/one)**(arb(3)/5),None,peak,F(-4))]
    t21=(one/two)**(arb(3)/5);t10=(peak/one)**(arb(3)/5)
    if less(t21,t10):return [(z,t21,two,F(-2,3)),(t21,t10,one,F(-7,3)),(t10,None,peak,F(-4))]
    t20=(peak/two)**(arb(3)/10)
    return [(z,t20,two,F(-2,3)),(t20,None,peak,F(-4))]


def layer_budget(peak,one,two,volume,eps,A):
    total=arb(0)
    for lo,hi,a,e in density_pieces(peak,one,two):
        if not less(lo,eps):break
        end=eps if hi is None else mn(hi,eps)
        # Lq <= min(V^(3/4) P, (At)^(3/4) P^(1/4)).
        tc=(a*volume/A)**(1/aq(1-e))
        split=mn(end,mx(lo,tc))
        if less(lo,split):total+=power_integral(A**(arb(3)/4)*a**(arb(1)/4),e/4+F(3,4),lo,split)
        if less(split,end):total+=power_integral(volume**(arb(3)/4)*a,e,split,end)
    return total


def global_mass_budget(peak,ones,twos,eps,A):
    # All breakpoints of the eight distinct density monomials are known
    # analytically. Between them the min/max branch is fixed.
    coeff=[peak]+ones+[v for v in twos if v is not None]
    expo=[F(-4)]+[F(-7,3)]*4+[F(-2,3)]*3
    choices=[[0,1],[0,2,5],[0,3,6],[0,4,7]]
    roots=[]
    for i,j in itertools.combinations(range(8),2):
        if expo[i]==expo[j]:continue
        t=(coeff[j]/coeff[i])**(1/aq(expo[i]-expo[j]))
        if less(t,eps):roots.append(t)
    roots.sort(key=lambda z:float(z.mid()))
    for l,u in zip(roots,roots[1:]):assert less(l,u)
    endpoints=[arb(0)]+roots+[eps]
    total=arb(0)
    for lo,hi in zip(endpoints,endpoints[1:]):
        mid=hi/2 if lo.is_zero() else (lo*hi).sqrt()
        vals=[a*mid**aq(e) for a,e in zip(coeff,expo)]
        active=[]
        for ids in choices:
            chosen=ids[0]
            for i in ids[1:]:
                if less(vals[i],vals[chosen]):chosen=i
            active.append(chosen)
        chosen=active[0]
        for i in active[1:]:
            if i!=chosen and less(vals[chosen],vals[i]):chosen=i
        total+=power_integral(A**(arb(3)/4)*coeff[chosen]**(arb(1)/4),expo[chosen]/4+F(3,4),lo,hi)
    return total,len(roots)


def main():
    inputs=REG['inputs'];assert REG['input_sha256']==hashlib.sha256(json.dumps(inputs,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    fc=S['final_constants'];A=arb(fc['A']);M=arb(fc['M']);peak=arb(fc['C6rad']);Q1=arb(fc['Cquotient']);B0=arb(fc['B0'])
    P4=arb(O['eight_peak_cases']['d4_dimension6']['peak_upper']);J2=16*arb(O['pair_coarea_sum']);K2=J2*P4
    b32=A*32**(arb(3)/2);b31=A*32**(arb(5)/4)*31**(arb(1)/4)
    B2=K2*b32*b31
    assert B0.overlaps(32*Q1*b32)
    # Four exact geometry checks; cone fraction is 1/(5+1).
    D=sp.eye(5)-sp.ones(5,5)/6
    assert D.det()==sp.Rational(1,6) and D.inv()==sp.eye(5)+sp.ones(5,5)
    u=sp.symbols('u',real=True);assert sp.integrate(u**5,(u,0,1))==sp.Rational(1,6)
    assert all(D.inv()[i,i]==2 for i in range(5))
    b=sp.symbols('b',real=True);assert sp.simplify(1+5*b*b-(1+5*b)**2/6-sp.Rational(5,6)*(1-b)**2)==0
    assert max(sp.Rational(k*(6-k),6) for k in range(1,7))==sp.Rational(3,2)
    K=4*arb(6).sqrt()*arb.pi()**2/45
    distances=[arb(0),arb(1)/32,arb(1)/8,arb(1)/2,(arb(3)/2).sqrt()]
    volumes=[K*distances[i+1]**5-K*distances[i]**5 for i in range(3)]
    volumes.append(1-K*distances[3]**5)
    assert sum(volumes).overlaps(arb(1)) and all(v>0 for v in volumes)
    rmins=[]
    for i in range(4):
        if i==3:rmins.append(arb(1));continue
        lower_mean=1-(arb(6)/5).sqrt()*distances[i+1]
        assert lower_mean>0
        rmins.append((1+5*lower_mean**2).sqrt())
    ones=[B0*r**(-arb(5)/2) for r in rmins]
    twos=[None]+[B2*rmins[i]**(-arb(5)/2)*distances[i]**(-arb(5)/2) for i in range(1,4)]
    cases=[];aggregates=[]
    for epss in inputs['epsilon']:
        eps=aq(epss);layer=[]
        for i in range(4):
            ub=layer_budget(peak,ones[i],twos[i],volumes[i],eps,A)
            layer.append(ub)
            cases.append({'layer':i+1,'epsilon':epss,'T_epsilon_Lq_upper':ub.str(65)})
        disjoint=sum(v**(arb(4)/3) for v in layer)**(arb(3)/4)
        global_mass,nroots=global_mass_budget(peak,ones,twos,eps,A)
        baseline=next(arb(v['T_epsilon_norm_upper'])for v in S['finite_epsilon_budgets'] if v['epsilon']==epss)
        selected=mn(mn(disjoint,global_mass),baseline)
        aggregates.append({'epsilon':epss,'disjoint_layer_bound':disjoint.str(65),
                           'global_mass_with_layer_maximum_bound':global_mass.str(65),
                           'frozen_wholeface_baseline':baseline.str(65),'selected_upper':selected.str(65),
                           'global_analytic_breakpoint_count':nroots,
                           'passes_1e_minus5':bool(selected<arb(1)/10**5),
                           'budget_improvement':(baseline/selected).str(65)})
    productD=arb(15)/4*B0**(arb(1)/5)*B2**(arb(4)/5)
    delta=arb(1)/1000
    cap_rmin=(1+5*(1-(arb(6)/5).sqrt()*delta)**2).sqrt()
    capI=arb(15)/7*K*productD**(arb(4)/3)*cap_rmin**(-arb(10)/3)*delta**(arb(7)/3)
    capnorm=capI**(arb(3)/4)
    assert capnorm<arb(4525)/10**9 and capI<arb(749)/10**10
    short001=arb(aggregates[1]['selected_upper'])
    assert short001<arb(18624)/10**9 and short001>arb(1)/10**5
    h001=next(arb(v['H_epsilon_norm_upper'])for v in S['finite_epsilon_budgets'] if v['epsilon']=='1/1000')
    derivedG=h001+short001
    assert derivedG<arb(1000514)/10**9
    result={'id':REG['id'],'status':'STABLE_RANK2_DISTANCE_INTERFACE_AND_EXACT_S1_GEOMETRY; NO_C6_BOUND',
      'input_sha256':REG['input_sha256'],'source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'registration_sha256':hashlib.sha256((HERE/'registration.json').read_bytes()).hexdigest(),
      'input_file_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [SRC,OLD]},
      'constants_read_from_frozen_sources':True,
      'rank2_constants':{k:v.str(65)for k,v in dict(b32=b32,b31=b31,J2=J2,K2=K2,B2=B2,B0=B0,D_product=productD).items()},
      'exact_geometry':{'rho_matrix_determinant':'1/6','inverse':'I+11T','positive_orthant_cone_fraction':'1/6',
                         'rho_ball_volume_coefficient':K.str(65),'volume_formula_valid_for':'0<=delta<=1/sqrt2',
                         'rho_max_squared':'3/2','four_volumes':[v.str(65)for v in volumes],
                         'four_radius_lower_bounds':[v.str(65)for v in rmins]},
      'twelve_layer_epsilon_cases':cases,'three_derived_wholeface_bounds':aggregates,
      'fixed_true_S1_cap':{'rho_radius':'1/1000','G_I_upper':capI.str(65),'G_Lq_upper':capnorm.str(65),
                          'passes_G_cap_norm_1e_minus5':bool(capnorm<arb(1)/10**5),
                          'strict_rational_norm_upper':'4525/1000000000','strict_rational_I_upper':'749/10000000000'},
      'derived_whole_G_using_frozen_H001':{'norm_upper':derivedG.str(65),'strict_rational_upper':'1000514/1000000000',
                                         'source':'F008 frozen H_(1/1000) norm plus new T_(1/1000) bound; no new whole-time quadrature'},
      'runtime':{'wall_seconds':time.monotonic()-START,'python':platform.python_version(),'arb_bits':ctx.prec,
                 'new_peak_evaluations':0,'MC_samples':0,'spatial_grid_nodes':0,'matlab_calls':0,'geometric_checks':4}}
    assert len(cases)==12
    (HERE/'off_strata_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'constants':result['rank2_constants'],'aggregates':aggregates,'cap':result['fixed_true_S1_cap'],'runtime':result['runtime']},indent=2))


if __name__=='__main__':main()
