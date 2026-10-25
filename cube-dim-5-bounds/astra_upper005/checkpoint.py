"""A continuous linear-family ceiling and three necessary quadratic corner checks."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, floor
import ast
import hashlib
import importlib.util
import json
import sys
import time
from flint import arb, ctx

ROOT=Path(__file__).resolve().parent


def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pmul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def ppow(a,n):
    out=[F(1)]
    for _ in range(n):out=pmul(out,a)
    return out


def bernstein(power):
    n=len(power)-1
    return [sum((power[j]*F(comb(k,j),comb(n,j)) for j in range(k+1)),F(0)) for k in range(n+1)]


def certify_positive(power):
    initial=bernstein(power)
    stack=[(F(0),F(1),0,initial)]
    leaves=[]
    while stack:
        l,r,depth,coeff=stack.pop()
        if min(coeff)>0:
            leaves.append({'left':str(l),'right':str(r),'depth':depth,'minimum_coefficient':str(min(coeff))})
        elif depth==10:
            return {'certified':False,'failed_interval':[str(l),str(r)]}
        else:
            left,right=[coeff[0]],[coeff[-1]]
            row=coeff
            for _ in range(len(coeff)-1):
                row=[(a+b)/2 for a,b in zip(row,row[1:])]
                left.append(row[0]);right.append(row[-1])
            mid=(l+r)/2
            stack.extend([(mid,r,depth+1,list(reversed(right))),(l,mid,depth+1,left)])
    leaves.sort(key=lambda q:F(q['left']))
    assert leaves[0]['left']=='0' and leaves[-1]['right']=='1'
    assert all(a['right']==b['left'] for a,b in zip(leaves,leaves[1:]))
    return {'certified':True,'degree':len(power)-1,'partition':leaves,'power_coefficients':[str(x) for x in power]}


def main():
    began=time.monotonic();ctx.prec=256
    p003=ROOT.parent/'astra_upper003'/'verification.json'
    p004=ROOT.parent/'astra_upper004'/'corner_checkpoint.json'
    data003=json.loads(p003.read_text());data004=json.loads(p004.read_text())
    constants={k:F(data003['constants'][k]) for k in ['A_lower','B_upper','c_lower','A_upper','B_lower','c_upper']}
    d1=list(map(F,data004['moments']['E_d_W_power']))
    d2=list(map(F,data004['moments']['E_d_squared_W_power']))
    b=list(map(F,data004['moments']['E_b_corner_W_power']))
    U=[F(comb(7,k),2**k)*(constants['A_upper']*d1[k]-constants['B_lower']*d2[k]) for k in range(8)]
    V=[F(comb(11,k),2**k)*b[k] for k in range(12)]
    vp,up=ppow(V,7),ppow(U,11)
    residual=[F(0)]*max(len(vp),len(up))
    for i,x in enumerate(vp):residual[i]+=F(1,20000)**4*x
    for i,x in enumerate(up):residual[i]-=constants['c_upper']**4*x
    exclusion=certify_positive(residual)
    exclusion.update({'candidate_family':'Z=a(1+eta W), all a>=0, 0<=eta<=1/2, fixed quadratic ell',
                      'target_ceiling':'1/20000','variable':'eta=u/2'})
    source=ROOT.parent/'astra_upper004'/'corner_checkpoint.py'
    module=load_module('upper004_corner_locked',source)
    module.ORDER=22
    one=[F(0)]*23;two=[F(0)]*23;bm=[F(0)]*23
    for key,weight in data004['geometry']['basis_group_coefficients'].items():
        profile=ast.literal_eval(key)
        one=[a+weight*z for a,z in zip(one,module.mixed(profile,(0,0,0)))]
    for key,weight in data004['geometry']['pair_pattern_weights'].items():
        ones,twos=ast.literal_eval(key)
        two=[a+weight*z for a,z in zip(two,module.mixed(ones,twos))]
    for p in range(5):
        weight=5*comb(4,p)*4**p*12**(4-p)
        bm=[a+weight*z for a,z in zip(bm,module.mixed((0,p,4-p),(0,0,0)))]
    assert one[:12]==d1 and two[:12]==d2 and bm[:12]==b
    ell=[constants['A_lower']*x-constants['B_upper']*y for x,y in zip(one,two)]
    c=3/(2**(arb(15)/2)*arb.pi()**(arb(5)/2))
    candidates=[]
    for kappa in [F(1,512),F(1,256),F(1,128)]:
        f=[F(1),F(1,8),-kappa]
        endpoint_values=[sum(f[i]*F(w)**i for i in range(3)) for w in [-2,10]]
        assert min(endpoint_values)>0
        f7,f11=ppow(f,7),ppow(f,11)
        cu=sum((x*ell[i] for i,x in enumerate(f7)),F(0))
        cv=sum((x*bm[i] for i,x in enumerate(f11)),F(0))
        assert cu>0 and cv>0
        a_star=(module.ar(cu)/module.ar(cv))**(arb(1)/4)
        a=F(floor(float(a_star)*1000),1000)
        assert module.ar(a)<=a_star<module.ar(a+F(1,1000))
        val=F(11,4)*cu*a**7-F(7,4)*cv*a**11
        candidates.append({'kappa':str(kappa),'f_coefficients':[str(x) for x in f],
                           'f_endpoint_values':[str(x) for x in endpoint_values],
                           'a_selected':str(a),'a_star_Arb':str(a_star),'exact_U':str(cu),'exact_V':str(cv),
                           'corner_value':str(val),'corner_g_Arb':str(c*module.ar(val)),
                           'all_a_corner_ceiling_Arb':str(c*module.ar(cu)**(arb(11)/4)/module.ar(cv)**(arb(7)/4)),
                           'passes_necessary_screen':bool(c*module.ar(val)>arb(7)/156250),
                           'global_certificate_obtained':False})
    out={'status':'continuous_linear_family_exclusion_and_quadratic_necessary_corner_checks',
         'new_global_g':False,'new_C5_upper_bound':False,'linear_family_exclusion':exclusion,
         'moments_to_order_22':{'E_d_W_power':[str(x) for x in one],'E_d_squared_W_power':[str(x) for x in two],'E_b_corner_W_power':[str(x) for x in bm]},
         'constants':{k:str(v) for k,v in constants.items()},'quadratic_candidates':candidates,
         'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [p003,p004,source,ROOT/'registration.md',Path(__file__)]},
         'runtime':{'python':sys.executable,'arb_bits':ctx.prec,'seconds':time.monotonic()-began}}
    (ROOT/'checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'exclusion':{k:v for k,v in exclusion.items() if k not in ['power_coefficients','partition']},
                     'Bernstein_leaves':len(exclusion.get('partition',[])),'candidates':candidates,'runtime':out['runtime']},indent=2))


if __name__=='__main__':main()
