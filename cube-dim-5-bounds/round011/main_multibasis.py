"""Integer multibasis replay; no producer certificate or kernel execution."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import importlib.util,json,hashlib,time,os,argparse
P=Path(__file__).resolve().parents[1]
EXPECTED='b79c134fdcec283c120d464cf78f131c669c905ac5cdacda309a34f07a9fb29b'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('fixed010',P/'round010/check_cover_integer.py')
fixed=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixed)
def margin(v,gd,delta,multi):
    dn,dp=fixed.dyadic(delta);p=max(fixed.dyadic(x)[1] for x in v)
    vi=[int(x*(1<<p)) for x in v]
    slopes=[fixed.dyadic(gd[j][0 if vi[j]>=0 else 1]) for j in range(4)]
    K=max(dp,2*p,*(sp+p for sn,sp in slopes))
    num=9*dn*(1<<(K-dp))
    num+=sum(9*sn*z*(1<<(K-sp-p)) for (sn,sp),z in zip(slopes,vi))
    summed=sum(abs(z) for z in vi) if multi else sum(vi)
    num-=1024*(9*sum(z*z for z in vi)-summed*summed)*(1<<(K-2*p))
    return F(num,9*(1<<K))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--state',type=Path,required=True)
    ap.add_argument('--producer-results',type=Path,required=True);a=ap.parse_args();began=time.monotonic()
    assert sha(a.state)==EXPECTED;d=json.loads(a.state.read_text());producer=json.loads(a.producer_results.read_text())
    assert d['epsilon']=='1/64' and d['target']=='1/10000' and len(d['splits'])==162
    pc={tuple(c['node']):c for c in producer['certificates']};output=[]
    for row in d['splits']:
        i,path=row['node'];assert path=='';c=row['failed_certificate']
        x=list(map(F,c['anchor'].split(',')));assert len(x)==5 and x[0]==1
        lo=[F(z,32)-v for z,v in zip(d['roots'][i],x[1:])];hi=[v+F(1,32) for v in lo]
        assert all(l<0<u for l,u in zip(lo,hi))
        vals=[list(map(F,z)) for z in d['cache'][c['anchor']]['values']];HL,HU=vals[0]
        delta=F(c['Delta_lower']);ll,lu=fixed.log_fixed(HL*10000)
        assert delta.numerator*fixed.SCALE<=ll*delta.denominator
        gd=[list(map(F,z)) for z in c['loggradient_directed']]
        for (gl,gu),(dl,du) in zip(gd,vals[1:]):
            rr=[n/z for n in (dl,du) for z in (HL,HU)]
            assert gl<=min(rr)<=max(rr)<=gu
        old=min(margin(list(v),gd,delta,False) for v in product(*zip(lo,hi)))
        assert old==F(c['minimum_log_margin']) and old<=0
        new=min(margin(list(v),gd,delta,True) for v in product(*[(l,F(0),u) for l,u in zip(lo,hi)]))
        assert new>0 and new>old
        q=pc[(i,path)];assert new==F(q['new_minimum'])
        output.append({'node':[i,path],'new_minimum':str(new),'old_minimum':str(old)})
        assert time.monotonic()-began<30
    assert len(output)==162 and len(pc)==162
    out={'status':'PASS_primary_integer_81_point_qualification','state_sha256':EXPECTED,
         'producer_results_sha256':sha(a.producer_results),'boxes':162,'new_integer_points':13122,
         'old_integer_points':2592,'independent_fixed_logs':162,'new_kernel_calls':0,
         'new_gradient_calls':0,'Arb_calls':0,'all_strict':True,'certificates':output,
         'global_g':False,'new_C5':False,'independence':'Uses main010 fixed log arithmetic; underlying saved H and gradients shared with046. No056 functions executed.',
         'source_hashes':{str(q.relative_to(P)):sha(q) for q in [Path(__file__),P/'round011/main_multibasis_registration.md',P/'round010/check_cover_integer.py']},
         'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began}}
    (P/'round011/main_multibasis.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='certificates'},indent=2))
if __name__=='__main__':main()
