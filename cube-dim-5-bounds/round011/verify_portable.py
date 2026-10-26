"""Portable exact local replay of061, no private path or kernel dependency."""
import argparse,gzip,json,time,hashlib,os
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import main_multibasis as core
P=Path(__file__).resolve().parents[1]
LOCK='e25c6cbcad0fabcdc8df5ade205de098c8c952e18f73a9d57a1d888cc40291f4'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,default=P/'upper_portable065/replay_001/records.json.gz')
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert not args.output.exists()
    began=time.monotonic();raw=args.bundle.read_bytes();assert hashlib.sha256(raw).hexdigest()==LOCK
    data=json.loads(gzip.decompress(raw));rows=data['records'];assert len(rows)==2048
    fixed=core.fixed;fixed.SCALE=1<<896;fixed.N=224;fixed.LOG2=fixed.log_atanh(F(1,3))
    unique=set();worst=None
    for row in rows:
        root=row['root_integers'];i,path=row['node'];assert (i,path) not in unique;unique.add((i,path))
        assert len(root)==4 and all(type(x) is int and 0<=x<=31 for x in root)
        assert root==sorted(root,reverse=True) and all(x in '01' for x in path)
        lo=[];hi=[]
        for j,r in enumerate(root):
            s=path[j::4];m=1<<len(s);b=int(s,2) if s else 0
            lo.append(F(r*m+b,32*m));hi.append(F(r*m+b+1,32*m))
        assert lo==list(map(F,row['box_lower'])) and hi==list(map(F,row['box_upper']))
        rec=row['cache_record'];c=row['certificate'];point=list(map(F,rec['point']))
        assert point[0]==1 and ','.join(map(str,point))==c['anchor']
        assert all(l<x<u for l,x,u in zip(lo,point[1:],hi))
        vals=[list(map(F,z)) for z in rec['values']];assert len(vals)==5 and all(l<=u for l,u in vals)
        HL,HU=vals[0];assert HL>F(1,10000);delta=F(c['Delta_lower'])
        l,u=fixed.log_fixed(HL*10000);assert delta.numerator*fixed.SCALE<=l*delta.denominator
        gd=[list(map(F,z)) for z in c['loggradient_directed']];assert len(gd)==4
        for (gl,gu),(dl,du) in zip(gd,vals[1:]):
            q=[n/d for n in (dl,du) for d in (HL,HU)];assert gl<=min(q)<=max(q)<=gu
        a=[l-x for l,x in zip(lo,point[1:])];b=[u-x for u,x in zip(hi,point[1:])]
        old=min(core.margin(list(v),gd,delta,False) for v in product(*zip(a,b)))
        assert old==F(row['old_minimum_log_margin']) and old<=0
        new=min(core.margin(list(v),gd,delta,True) for v in product(*[(l,F(0),u) for l,u in zip(a,b)]))
        assert new==F(c['minimum_log_margin']) and new>0 and c['checked_points']==81
        if worst is None or new<worst:worst=new
        assert time.monotonic()-began<120
    out={'status':'PASS_portable_primary_integer_2048_regions','bundle_sha256':LOCK,'regions':2048,
         'new_points':165888,'independent_Delta_logs':2048,'log_fixed_bits':896,'atanh_terms':224,
         'minimum_log_margin':str(worst),'private_paths_required':False,'new_kernel_calls':0,
         'global_g':False,'new_C5':False,'shared_dependency':'Saved046 H/gradient intervals; main010 fixed log and011 integer formula reused. Does not repeat kernel values or globaltree audit.',
         'source_hashes':{str(q.relative_to(P)):hashlib.sha256(q.read_bytes()).hexdigest() for q in [Path(__file__),P/'round011/main_multibasis.py',P/'round010/check_cover_integer.py',P/'round011/portable_verification_registration.md']},
         'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began}}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
