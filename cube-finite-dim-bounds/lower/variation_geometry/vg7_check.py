"""Twelve fixed exact identities; local toy, not a new maximal lower bound."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
import sympy as sy
from vg7_joint import parse_key,descriptor,compile_rational
from normal_intervals import signed_length

P=Path(__file__).resolve().parent


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x


def key_from_atoms(points):
    return tuple(z for point in points for z in point)*2


def main():
    passed=[]
    def check(name,truth):
        assert truth,name
        passed.append(name)
    points=[(0,0,0,0),(1,-1,0,0),(1,0,0,0)]
    key=key_from_atoms(points)
    phi=[(0,0,0,0),(0,F(7,8),0,0),(0,F(1,8),0,0)]
    check('24-key periodic counts',parse_key(key)[2]==(1,1,1))
    wider=list(key);wider[16:20]=[2,0,0,0]
    desc=descriptor(wider,[F(1,2),F(1,2)])
    check('captured mass differs from full cell mass',desc['counts']==(1,4,1)
          and desc['captured_mass_slope']==F(5,2) and desc['complete_cell_mass_slope']==1)
    n,s,t=sy.symbols('n s t',positive=True)
    radius=(n+s*t)**sy.Rational(1,4)/2
    check('symbolic common-radius identity',
          sy.simplify(sy.diff(radius,t).subs(t,0)-s/(8*n**sy.Rational(3,4)))==0
          and sy.simplify(sy.diff(radius,t).subs(t,0)-radius.subs(t,0)*s/(4*n))==0)
    eta=F(1,256);d=[F(3,16),0,0,0]
    joint=compile_rational(key,[F(1,2),F(1,2)],phi,[1]*4,d,eta)
    key1=list(key);key1[8:12]=[0]*4;key1[20:24]=[-1]*4
    key2=list(key);key2[4:8]=[0]*4;key2[16:20]=[-1]*4
    single1=compile_rational(key1,[F(1,2),0],phi,[1]*4,d,eta)
    single2=compile_rational(key2,[0,F(1,2)],phi,[1]*4,d,eta)
    check('joint shared radius and normal velocities',joint['r1']==F(1,8)
          and joint['lower'][0].velocity==F(-1,32) and joint['upper'][0].velocity==F(1,32))
    check('both single directions have empty born interval',
          single1['lower'][0].velocity==single2['lower'][0].velocity==F(1,32)
          and single1['upper'][0].velocity==single2['upper'][0].velocity==F(-1,32))
    L=F(-3,32);R=F(3,32)
    D0=signed_length([L],[R]);D1=signed_length([L],[R],[(F(1,32),F(-1,32))])
    D2=D1;D12=signed_length([L],[R],[(F(-1,32),F(1,32))])
    check('contact mixed increment is nonzero',D0==D1==D2==F(-3,16)
          and D12==F(-1,8) and D12-D1-D2+D0==F(1,16))
    area=F(1,8);V0=4*area
    G0=area*D0+V0*d[0]
    G1=area*D1+V0*(d[0]-F(1,2))
    G2=G1;G12=area*D12+V0*(d[0]-1)
    check('complete mass and Jacobian preserve mixed increment',G12-G1-G2+G0==F(1,128))
    tau=[(0,0,0,0),(F(1,8),0,0,0),(F(-1,8),0,0,0)]
    switching=compile_rational(key,[F(1,2)]*2,phi,[1]*4,d,eta,tau)
    check('tied max winner moves center and kills slab',
          switching['lower'][0].velocity==F(3,32)
          and switching['upper'][0].velocity==F(1,32))
    empty=(0,)*12+(-1,)*12
    check('canonical empty key',parse_key(empty)[2]==(0,0,0))
    pure=list(key);pure[0:4]=[0]*4;pure[12:16]=[-1]*4
    pdesc=descriptor(pure,[F(1,2)]*2)
    check('pure common cluster has different scaling',pdesc['regime']=='pure'
          and pdesc['captured_mass_slope']==1 and 'radius_velocity_factor' not in pdesc)
    root=F(2049,2048);sigma=root**4-1;h=1+F(3,16)*sigma
    old_front=F(1,2)/h-F(1,2);born=root/(2*h)-F(1,2)
    length0=4+2*old_front;length12=length0+2*born
    Vbase=area*length0;Vjoint=area*length12
    check('one exact finite perturbation validates joint-only gain',0<sigma<eta
          and h**4>1+sigma/2 and 0<born<-old_front and Vjoint-Vbase==2*area*born)
    wrapped=list(key);wrapped[5]+=1;wrapped[17]+=1
    rephi=[phi[0],(0,F(-1,8),0,0),phi[2]]
    same=compile_rational(wrapped,[F(1,2)]*2,rephi,[1]*4,d,eta)
    check('integer phase-wrap compensation preserves receiver',
          same['lower']==joint['lower'] and same['upper']==joint['upper'])
    assert len(passed)==12
    out=dict(status='12 exact identities PASS; local toy only, not a new bound',
             identities=passed,key=key,weights=[F(1,2)]*2,period_direction=d,eta=eta,
             joint_radius_velocity=joint['r1'],joint_normal_interval=[F(-1,32),F(1,32)],
             single_normal_interval=[F(1,32),F(-1,32)],
             normal_derivatives=dict(base=D0,single1=D1,single2=D2,joint=D12),
             mixed_normal_increment=F(1,16),tangential_patch_area=area,
             quotient_derivatives=dict(base=G0,single1=G1,single2=G2,joint=G12),
             mixed_full_ratio_increment=F(1,128),
             switching_normal_interval=[F(3,32),F(1,32)],
             finite=dict(sigma=sigma,common_root=root,period_factor=h,
                         normal_old_front=old_front,normal_born=born,
                         normalized_volume_no_satellite_gain=Vbase,normalized_joint_volume=Vjoint,
                         exact_joint_extra_volume=Vjoint-Vbase,
                         complete_joint_ratio=h*Vjoint/(1+sigma)),
             source_sha256={name:hashlib.sha256((P/name).read_bytes()).hexdigest()
                            for name in ['vg7_joint.py','vg7_check.py','vg6_bounds.py','normal_intervals.py']})
    (P/'vg7_results.json').write_text(json.dumps(encode(out),indent=2)+'\n')
    print(json.dumps(encode(out),indent=2))


if __name__=='__main__':main()
