#!/usr/bin/env python3
"""Twelve fixed structural configurations; no kernels or producer execution."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,combinations
from collections import defaultdict,Counter
from functools import reduce
import hashlib,json,math,os,platform,time
import sympy as sp

HERE=Path(__file__).resolve().parent
TOP=Path('/Users/zhengzhihao/Desktop/T')
SRC=TOP/'Math-Research-cube6-theory/cube-dim-6-bounds/upper/round020/pair_orbit_frame'
MAIN=TOP/'Math-Research-cube-dim-6/cube-dim-6-bounds'
INPUTS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    INPUTS[str(p)]=sha(p);return json.loads(p.read_text())
def pin(p,h):
    INPUTS[str(p)]=sha(p);assert INPUTS[str(p)]==h,str(p)
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def wedge0(x,y):return tuple(x[i]*y[j]-x[j]*y[i] for i,j in combinations(range(6),2))
def wedge(x,y):
    raw=wedge0(x,y);g=abs(reduce(math.gcd,raw,0));assert g
    raw=tuple(v//g for v in raw)
    return tuple(-v for v in raw) if next(v for v in raw if v)<0 else raw
def parse(path):
    it=iter(map(int,path.read_text().split()));rows=[]
    for _ in range(next(it)):
        d=next(it);rays=[]
        for _ in range(6):rays.append((next(it),tuple(next(it) for _ in range(6))))
        rows.append((d,rays))
    assert next(it,None) is None;return rows

def main():
    started=time.perf_counter();meta=read(SRC/'input_manifest.json');saved=read(SRC/'results.json')
    for path,h in meta['input_sha256'].items():pin(Path(path),h)
    pin(SRC/'build_frame.py',meta['source_sha256'])
    assert saved['source_sha256']==meta['source_sha256']
    for name,row in saved['exported_artifacts'].items():pin(SRC/name,row['sha256'])
    # Author's final report and manifest pins, received after the writer stopped.
    pin(SRC/'report.md','e4553f5cde36f38a7eb5c6c781a07b66005084ba2edf9294355a87a8b538c6d4')
    pin(SRC/'input_manifest.json','ec5825006729b8bbcb524778bba2da48910740a8da98bdf7f6b95a5d7cc3736b')
    pin(SRC/'results.json','128e8093912ec5f1d6a45f74e120061cf0e4cfdea2e60d6d4e430369043bbaaa')
    body=[(1,)+v for v in product((-1,1),repeat=5)]
    axes=[tuple(int(i==j) for j in range(6)) for i in range(6)]
    pair=[tuple(axes[i][k]+s*axes[j][k] for k in range(6))
          for i,j in combinations(range(6),2) for s in (1,-1)]
    frame=body+pair;checks=[]
    def record(name,details):checks.append({'configuration':len(checks)+1,'name':name,'PASS':True,'details':details})
    # 1: one beta=1 physical frame, Gram/rank and marginal time units.
    Gb=sp.Matrix(6,6,lambda i,j:sum(v[i]*v[j] for v in body))
    Gp=sp.Matrix(6,6,lambda i,j:sum(v[i]*v[j] for v in pair))
    assert Gb==32*sp.eye(6) and Gp==10*sp.eye(6) and (Gb+Gp).rank()==6
    assert len(frame)==62 and [sum(abs(v[i]) for v in frame) for i in range(6)]==[42]*6
    assert 42*Q(32,875)==32*Q(6,125)==Q(192,125)
    record('beta1_frame_Gram_rank_marginal_units',{'columns':62,'rank':6,'body_Gram':'32I','pair_Gram':'10I',
           'Gram':'42I','coordinate_rate':42,'normalized_T':'32/875','T_evaluated':False})
    # 2: symbolic linear identity on the full closed sorted chamber.
    xi=sp.symbols('x1:7');ell=(10,8,6,4,2,0)
    expr=sum((xi[i]+xi[j])+(xi[i]-xi[j]) for i,j in combinations(range(6),2))
    assert sp.expand(expr-sum(a*x for a,x in zip(ell,xi)))==0
    record('whole_chamber_pair_symbol',{'ell':ell,'symbolic_residual':'0','ties_and_zero_allowed':True})
    # 3: existing fan artifact transform; no old coverage verifier rerun.
    oldpath=MAIN/'upper/round003/direct_kernel/fan_integer_input.txt'
    old=parse(oldpath);new=parse(SRC/'pair_orbit_beta1_integer_fan.txt');fan=read(SRC/'pair_orbit_beta1_fan.json')
    coverage=read(MAIN/'upper/round002/fan/verify_exact_fan_results.json')
    bridge=read(MAIN/'upper/round004/primary/integer_fan_bridge.json')
    assert coverage['status']=='PASS' and coverage['full_cells']==21
    assert bridge['status']=='PASS_exact_ray_rate_determinant_scaling_bridge'
    assert bridge['integer_input_sha256']==INPUTS[str(oldpath)]
    assert bridge['source_sha256']==coverage['input_sha256']
    assert len(old)==len(new)==len(fan['cones'])==21
    rates=[]
    for (d0,r0s),(d1,r1s),jc in zip(old,new,fan['cones']):
        assert d0==d1==jc['determinant_integer'] and d1>0
        qs=[]
        for (q0,r0),(q1,r1),jr in zip(r0s,r1s,jc['rows']):
            assert r0==r1==tuple(jr['ray_integer'])
            assert all(r0[i]>=r0[i+1] for i in range(5)) and r0[-1]>=0 and r0[0]>0
            assert q0==sum(abs(dot(v,r0)) for v in body)==jr['old_rate_integer']
            inc=sum(abs(dot(v,r0)) for v in pair)
            assert inc==dot(ell,r0)==jr['pair_increment_integer']
            assert q1==q0+inc==jr['new_rate_integer']>q0
            assert all(r0[i]+r0[j]>=0 and r0[i]-r0[j]>=0 for i,j in combinations(range(6),2))
            qs.append(q1);rates.append(q1)
        assert Q(jc['new_alpha_exact'])==Q(d1,math.prod(qs))
    assert len(rates)==126
    record('frozen_fan_rate_update',{'cones':21,'ray_occurrences':126,'new_rate_range':[min(rates),max(rates)],
           'rays_and_determinants_unchanged':True,'old_baseline_rerun':False})
    # 4--8: five canonical plane configurations, counted with exact wedge keys.
    def members(a,b):
        key=wedge(a,b)
        return [i for i,w in enumerate(frame) if not any(wedge0(a,w)) or wedge(a,w)==key]
    v=(1,)*6;p=tuple(a+b for a,b in zip(axes[0],axes[1]))
    minus=tuple(a-b for a,b in zip(axes[0],axes[1]));w=tuple(a-2*b for a,b in zip(v,p))
    assert wedge(v,p)==wedge(v,w) and len(members(v,p))==3
    record('matching_body_pair',{'frame_lines':3,'old_body_alias':True})
    assert dot(v,minus)==0 and len(members(v,minus))==2
    record('orthogonal_body_pair',{'frame_lines':2,'body_dot_pair':0})
    assert wedge(p,minus)==wedge(axes[0],axes[1]) and len(members(p,minus))==2
    record('same_support_pair_pair',{'frame_lines':2,'coordinate_plane':True})
    q=tuple(a+b for a,b in zip(axes[0],axes[2]));third=tuple(a-b for a,b in zip(axes[1],axes[2]))
    assert tuple(a-b-c for a,b,c in zip(p,q,third))==(0,)*6 and len(members(p,q))==3
    record('one_overlap_pair_triangle',{'frame_lines':3,'dependent_triple':True})
    q=tuple(a-b for a,b in zip(axes[2],axes[3]))
    assert dot(p,q)==0 and len(members(p,q))==2
    record('disjoint_support_pair_pair',{'frame_lines':2})
    # 9: complete pair classification for this same physical frame.
    planes=defaultdict(list);pc=Counter()
    for i,j in combinations(range(62),2):
        a,b=frame[i],frame[j];key=wedge(a,b)
        if j<32:kind='body_body'
        elif i<32:kind='body_pair_matching' if dot(a,b) else 'body_pair_orthogonal'
        else:kind={2:'pair_same',3:'pair_triangle',4:'pair_disjoint'}[sum(bool(x or y) for x,y in zip(a,b))]
        pc[kind]+=1;planes[key].append((i,j,kind))
    assert dict(pc)=={'body_body':496,'body_pair_matching':480,'body_pair_orthogonal':480,'pair_same':15,'pair_triangle':240,'pair_disjoint':180}
    classes=Counter();multiplicity=Counter();deps=0
    for records in planes.values():
        kinds={r[2] for r in records};lines={i for r in records for i in r[:2]}
        assert len(records)==math.comb(len(lines),2)
        if 'body_body' in kinds:
            assert kinds<={'body_body','body_pair_matching'};kind='old_body'
        else:assert len(kinds)==1;kind=next(iter(kinds))
        classes[kind]+=1;multiplicity[len(lines)]+=1;deps+=math.comb(len(lines),3)
    assert dict(classes)=={'old_body':496,'body_pair_orthogonal':480,'pair_same':15,'pair_triangle':80,'pair_disjoint':180}
    assert len(planes)==1251 and sum(pc.values())==1891 and deps==320
    assert dict(multiplicity)=={2:931,3:320}
    record('complete_1891_pair_classification',{'direction_pairs':1891,'planes':1251,'new_planes':755,
           'pair_classes':dict(pc),'plane_classes':dict(classes),'plane_line_counts':dict(multiplicity),
           'dependent_triples':320,'triple_enumeration':False})
    # 10--12: the three already registered symbol vectors, with no additions.
    for x in [(1,1,1,1,1,1),(3,2,1,0,0,0),(-5,1,-3,2,0,-1)]:
        pv=sum(abs(dot(v,x)) for v in pair);bv=sum(abs(dot(v,x)) for v in body)
        assert pv==dot(ell,sorted(map(abs,x),reverse=True))
        record('registered_symbol_vector',{'x':x,'body':bv,'pair':pv,'beta1_total':bv+pv})
    assert len(checks)==12
    out={'status':'PASS_PAIR_FRAME_THEORY_FAN_AND_S2_STRUCTURE_ONLY','pid':os.getpid(),
         'python':platform.python_version(),'sympy':sp.__version__,'configurations':checks,'configuration_count':12,
         'new_kernel_calls':0,'MC_calls':0,'MATLAB_calls':0,'old_baseline_replayed':False,
         'new_global_m_certified':False,'new_C6_bound':False,'input_sha256':INPUTS,
         'script_sha256':sha(Path(__file__)),'elapsed_seconds':time.perf_counter()-started}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','pid','elapsed_seconds','configuration_count','new_kernel_calls']}))
if __name__=='__main__':main()
