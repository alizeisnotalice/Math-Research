"""One registered new-frame point; never rerun after execution.json exists."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,subprocess,itertools,time,datetime,sys
sys.dont_write_bytecode=True
H=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T')
A=ROOT/'Math-Research-cube6-theory/cube-dim-6-bounds/upper/round020/pair_orbit_frame'
C=ROOT/'Math-Research-cube6-audit/cube-dim-6-bounds/upper/round014/generic_short_engine'
F=ROOT/'Math-Research-cube6-literature/cube-dim-6-bounds/audit/round020/pair_frame_review'
BIN=C/'scratch/generic_short';FAN=A/'pair_orbit_beta1_integer_fan.txt'
NS=[100000000000000,99999998999999,99999997999988,90268521812170,53714088833376,157107163435]
DEN=10**14;T=Q(32,875)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def interval(d):
    assert set(d)=={'lower_mantissa','upper_mantissa','common_binary_exponent'}
    e=int(d['common_binary_exponent']);assert abs(e)<10000
    q=Q(2)**e;l=int(d['lower_mantissa'])*q;u=int(d['upper_mantissa'])*q
    assert l<=u;return l,u
def directions():
    v=[(1,*s)for s in itertools.product((-1,1),repeat=5)]
    for i in range(6):
        for j in range(i+1,6):
            for s in [1,-1]:
                a=[0]*6;a[i]=1;a[j]=s;v.append(tuple(a))
    assert len(v)==62;return v
def check(raw):
    assert raw['status']=='generic_nonzero_supported' and raw['precision_bits']==512
    assert raw['full_group_terms']==967680 and raw['computed_sign_pairs']==483840
    assert raw['zero_node_terms']==0 and 0<=raw['confluent_terms']<=967680
    k,p=interval(raw['K']),interval(raw['p_T']);g=[interval(x)for x in raw['gradient']]
    assert len(g)==6 and k[0]>0 and p[0]>0
    assert interval(raw['Euler_x_grad_plus5K_plusTp'])[0]<=0<=interval(raw['Euler_x_grad_plus5K_plusTp'])[1]
    for z in [k,p,*g]:assert z[1]-z[0]<max(Q(1),abs(sum(z)/2))*Q(2)**-200
    x=[Q(n,DEN)for n in NS];el=5*k[0]+T*p[0];eu=5*k[1]+T*p[1]
    for q,z in zip(x,g):el+=q*z[0];eu+=q*z[1]
    assert el<=0<=eu
    v=directions();n2=sum(q*q for q in x);ds=[]
    for j,b in enumerate(v):
        bb=sum(t*t for t in b);xb=sum(t*q for t,q in zip(b,x))
        for i,a in enumerate(v[:j]):
            aa=sum(t*t for t in a);ab=sum(t*s for t,s in zip(a,b));xa=sum(t*q for t,q in zip(a,x))
            det=aa*bb-ab*ab;assert det>0
            ds.append((n2-(bb*xa*xa-2*ab*xa*xb+aa*xb*xb)/det,i,j))
    assert len(ds)==1891;d2,i,j=min(ds);assert d2>0
    return {'status':'POINT_VALIDATED_ONLY','K_exact_interval':list(map(str,k)),'T':str(T),
            'threshold64':str(T/4096),'threshold625':str(T/4000),
            'strictly_excludes64':k[1]<=T/4096,'strictly_excludes625':k[1]<T/4000,
            'necessary64_passes':k[0]>T/4096,
            'cost_interval':[str(T/(64*k[1])),str(T/(64*k[0]))],
            'q_derivative_interval':[str((T*p[0]-k[1])/T**2),str((T*p[1]-k[0])/T**2)],
            'independent_Euler_interval':[str(el),str(eu)],
            'true_new_frame_pair_count':1891,'S2_distance_squared':str(d2),'closest_pair':[i,j],
            'new_C6_bound':False,'global_m_certified':False}
def main():
    assert not (H/'execution.json').exists(),'already executed: no automatic rerun'
    expected={C/'generic_short.cpp':'247a1f6c6967ad52f34320b999818214488700651b46cb8e5feb3eac7064953b',BIN:'67119aa799041ed75b66f043444640a2a62556f58517cbdf3bd6f7cccd564f04',FAN:'4202b6c685567e51fbd53275add62b6ed978176c834b9fb3f3aaac1f1e5f8a0e',A/'report.md':'e4553f5cde36f38a7eb5c6c781a07b66005084ba2edf9294355a87a8b538c6d4',F/'report.md':'83a3b14e3ca67e6d8852f99487721f7ba22391eadf769347973f6cd37555a835'}
    for p,h in expected.items():assert p.is_file()and sha(p)==h
    libs=json.loads((ROOT/'Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round015/spatial_minimum_search/preparation.json').read_text())['inputs']
    for row in libs:
        if Path(row['path']).name in ['libflint.24.0.dylib','libgmp.10.dylib','libmpfr.6.dylib']:
            assert sha(row['path'])==row['sha256'];expected[Path(row['path'])]=row['sha256']
    # Exact grammar and all 126 rate updates are checked before the trusted binary.
    tokens=FAN.read_text().split();assert all(str(int(t))==t for t in tokens);data=list(map(int,tokens));assert data[0]==21
    cur=1;rays=[];vs=directions()
    for cone in range(21):
        det=data[cur];cur+=1;assert det>0
        for m in range(6):
            q=data[cur];r=tuple(data[cur+1:cur+7]);cur+=7
            assert r[0]>=r[1]>=r[2]>=r[3]>=r[4]>=r[5]>=0
            assert q==sum(abs(sum(a*b for a,b in zip(v,r)))for v in vs)
            assert DEN*q<2**63 and sum(abs(a)*DEN for a in r)<2**63
            rays.append(r)
    assert cur==len(data)
    reg={'id':'primary021-pair-orbit-unique-probe','time':datetime.datetime.now().astimezone().isoformat(),
         'frame':'body32 plus30 pair directions, beta1, independent standard Cauchy at one common t',
         'directions':[list(v)for v in vs],'T':str(T),'numerators':NS,'denominator':DEN,
         'native_limit':1,'process_wall_limit_seconds':120,'total_limit_seconds':1200,
         'stop':'after sole probe regardless of outcome; no new point/T/frame/optimizer/G/H/precision retry',
         'source_sha256':{str(p):h for p,h in expected.items()},'driver_sha256':sha(Path(__file__)),
         'point_certificate_inherited':False,'interface':'F020 + C014/F014 generic fan formula',
         'new_C6_bound_authorized':False}
    reg['frame_sha256']=hashlib.sha256(json.dumps({'directions':reg['directions'],'beta':'1','Fourier':'(2*pi)^(-6)','fan':sha(FAN)},sort_keys=True,separators=(',',':')).encode()).hexdigest()
    dump(H/'registration.json',reg)
    cmd=[str(BIN),'full',str(FAN),'32','875',str(DEN),*map(str,NS)]
    ex={'status':'REGISTERED_BEFORE_SINGLE_NATIVE','command':cmd,'registration_sha256':sha(H/'registration.json'),'start_epoch':time.time(),'parent_pid':__import__('os').getpid()};dump(H/'execution.json',ex)
    start=time.monotonic()
    with (H/'native.stdout.json').open('w')as out,(H/'native.stderr.txt').open('w')as err:
        child=subprocess.Popen(cmd,stdout=out,stderr=err,cwd=H);ex['native_pid']=child.pid;dump(H/'execution.json',ex)
        try:code=child.wait(timeout=120)
        except BaseException:
            if child.poll()is None:child.kill()
            code=child.wait();ex.update(status='FAILED_OR_INTERRUPTED_CHILD_WAITED',returncode=code);dump(H/'execution.json',ex);raise
    ex.update(status='SINGLE_NATIVE_EXITED_STOPPED',returncode=code,elapsed_seconds=time.monotonic()-start,stdout_sha256=sha(H/'native.stdout.json'),stderr_sha256=sha(H/'native.stderr.txt'),all_children_waited=True);dump(H/'execution.json',ex)
    assert code==0,(H/'native.stderr.txt').read_text()
    result=check(json.loads((H/'native.stdout.json').read_text()));dump(H/'checked_results.json',result)
    print(json.dumps({k:result[k]for k in ['status','T','strictly_excludes64','necessary64_passes','S2_distance_squared','new_C6_bound']}));print('K display',float(sum(map(Q,result['K_exact_interval']))/2),'pointcost display',float(Q(result['cost_interval'][0])),'wall',ex['elapsed_seconds'])
if __name__=='__main__':main()
