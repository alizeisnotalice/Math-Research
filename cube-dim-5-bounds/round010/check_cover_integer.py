"""Pure integer corner/log replay of the first048 partial cover, no kernel calls."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import comb
import argparse,hashlib,json,time,os
P=Path(__file__).resolve().parents[1];SCALE=1<<384;N=96
sha=lambda b:hashlib.sha256(b).hexdigest()
def ceildiv(a,b):return -(-a//b)
def log_atanh(z):
    assert 0<=z<=F(1,3)
    zl=z.numerator*SCALE//z.denominator;zu=ceildiv(z.numerator*SCALE,z.denominator)
    z2l=zl*zl//SCALE;z2u=ceildiv(zu*zu,SCALE);pl,pu=zl,zu;lo=hi=0
    for j in range(N):
        lo+=2*pl//(2*j+1);hi+=ceildiv(2*pu,2*j+1)
        pl=pl*z2l//SCALE;pu=ceildiv(pu*z2u,SCALE)
    hi+=ceildiv(2*pu*SCALE,(2*N+1)*(SCALE-z2u))
    return lo,hi
LOG2=log_atanh(F(1,3))
def log_fixed(x):
    assert x>0
    e=x.numerator.bit_length()-x.denominator.bit_length()
    y=x/(F(2)**e)
    if y<1:e-=1;y*=2
    if y>=2:e+=1;y/=2
    assert 1<=y<2
    l,u=log_atanh((y-1)/(y+1))
    return (l+e*LOG2[0],u+e*LOG2[1]) if e>=0 else (l+e*LOG2[1],u+e*LOG2[0])
def dyadic(x):
    x=F(x);d=x.denominator;assert d&(d-1)==0
    return x.numerator,d.bit_length()-1
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--state',required=True,type=Path)
    ap.add_argument('--expected-sha');ap.add_argument('--parent-state',required=True,type=Path);args=ap.parse_args();began=time.monotonic()
    raw=args.state.read_bytes();h=sha(raw)
    if args.expected_sha:assert h==args.expected_sha
    a=json.loads(raw);oldraw=args.parent_state.read_bytes();old=json.loads(oldraw)
    assert sha(oldraw)=='b23f4ee65bc7ce1d8bdee7f0050b4e492b453241a6594cce8e79f96e7bcfbb2f'
    roots=a['roots'];assert roots==old['roots'] and len(roots)==comb(35,4)==52360
    assert all(31>=r[0]>=r[1]>=r[2]>=r[3]>=0 for r in roots) and len({tuple(r) for r in roots})==52360
    assert a['epsilon']=='1/64' and a['target']=='1/10000' and not a['splits'] and not a['inflight'] and not a['resume_jobs']
    passed=a['passed'];queue=a['queue'];assert len(passed)==4224 and len(queue)==48136
    assert [x['node'] for x in passed]==[[i,''] for i in range(4224)]
    assert queue==[[i,''] for i in range(4224,52360)]
    for key,v in old['cache'].items():assert a['cache'][key]==v
    worst=None;corners=0;logcount=0
    for row in passed:
        i,path=row['node'];assert path=='';c=row['certificate'];point=list(map(F,c['anchor'].split(',')))
        assert len(point)==5 and point[0]==1
        lo=[F(t,32) for t in roots[i]];hi=[x+F(1,32) for x in lo]
        assert all(l<x<r for l,x,r in zip(lo,point[1:],hi))
        record=a['cache'][c['anchor']];vals=[list(map(F,z)) for z in record['values']];HL,HU=vals[0]
        assert 0<HL<=HU and HL>F(1,10000)
        delta=F(c['Delta_lower']);ll,lu=log_fixed(HL/F(1,10000));assert delta.numerator*SCALE<=ll*delta.denominator;logcount+=1
        gd=[list(map(F,z)) for z in c['loggradient_directed']]
        for (gl,gu),(dl,du) in zip(gd,vals[1:]):
            ratios=[d/v for d in (dl,du) for v in (HL,HU)]
            assert gl<=min(ratios)<=max(ratios)<=gu
        dn,dp=dyadic(delta);margins=[]
        for bits in product((0,1),repeat=4):
            v=[(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)]
            p=max(dyadic(x)[1] for x in v);vi=[int(x*(1<<p)) for x in v]
            slopes=[dyadic(gd[j][0 if vi[j]>=0 else 1]) for j in range(4)]
            K=max(dp,2*p,*(sp+p for sn,sp in slopes))
            num=9*dn*(1<<(K-dp))
            num+=sum(9*sn*z*(1<<(K-sp-p)) for (sn,sp),z in zip(slopes,vi))
            num-=1024*(9*sum(z*z for z in vi)-sum(vi)**2)*(1<<(K-2*p))
            assert num>0
            margins.append((F(num,9*(1<<K)),list(bits)));corners+=1
        value,bits=min(margins);assert value==F(c['minimum_log_margin']) and bits==c['worst_corner_bits']
        if worst is None or value<worst[0]:worst=(value,i,bits)
        assert time.monotonic()-began<30 and corners<=100000
    assert corners==67584 and logcount==4224
    out={'status':'PASS_primary_integer_partial_cover_certificate','state_sha256':h,
         'old_parent_state_sha256':sha(oldraw),'root_count':52360,'passed_leaves':4224,'remaining_leaves':48136,
         'all_corner_margins_strictly_positive':True,'integer_corners':corners,'independent_fixed_point_logs':logcount,
         'fixed_point_bits':384,'log_series_terms':N,'old128_cache_unchanged':True,
         'minimum_log_margin':str(worst[0]),'worst_root':worst[1],'worst_corner':worst[2],
         'independence':'No Arb, author certificate function, H, gradient or kernel replay; underlying H balls still depend on frozen046 API.',
         'global_g_certified':False,'new_C5_bound':False,
         'source_sha256':{str(p.relative_to(P)):sha(p.read_bytes()) for p in [Path(__file__),P/'round010/cover_integer_registration.md',P/'upper_cover048/report.md',P/'audit_harnack047/report.md']},
         'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began}}
    (P/'round010/cover_integer.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','passed_leaves','remaining_leaves','integer_corners','runtime']},indent=2))
if __name__=='__main__':main()
