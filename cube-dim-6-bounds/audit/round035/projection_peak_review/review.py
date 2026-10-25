#!/usr/bin/env python3
"""One independent full trace coverage audit; no producer code execution."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product,permutations
from math import factorial,prod
import gzip,hashlib,json,os,signal,sys,time
from flint import fmpz_mat,fmpq_mat,arb,ctx
sys.set_int_max_str_digits(0);ctx.prec=512
O=Path(__file__).resolve().parent
source=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round033/exact_projection_peaks')
U=source.parents[1]
FAN=U/'round003/direct_kernel/fan_integer_input.txt'
GEOM=U/'round002/fan/exact_cells_results.json'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
start=time.time();registration=json.loads((O/'input_manifest.json').read_text())
assert not (O/'execution_receipt.json').exists(),'one full enumeration only'
assert start<registration['started_unix']+1200
inputs=dict(registration['input_sha256'])
for row in json.loads((source/'input_manifest.json').read_text())['inputs']:inputs[row['path']]=row['sha256']
for path,h in inputs.items():assert sha(path)==h,path
assert sha(source/'selected_terms.json.gz')=='a2bcd0b5704570dc02643728f83cf8542c1efd5f5dc0da0689533bd7f50a1cc3'
inputs[str(Path(__file__))]=sha(Path(__file__))
receipt={'status':'RUNNING_SINGLE_INDEPENDENT_FULL_ENUMERATION','pid':os.getpid(),'start_unix':start,'full_enumeration_attempts':1,'native_G_H_calls':0,'completed_pattern_items':0}
def write(name,data):(O/name).write_text(json.dumps(data,indent=2)+'\n')
write('execution_receipt.json',receipt)
signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(TimeoutError('600 second full audit cap')))
signal.alarm(min(600,int(registration['started_unix']+1200-start)))

def fdet(mat):return int(fmpz_mat(mat).det())
def solve_and_det(M,b):
    # Rational Gaussian elimination, independent of producer Cramer formulas.
    n=len(M);a=[[F(x) for x in row]+[F(b[i])] for i,row in enumerate(M)];det=F(1)
    for col in range(n):
        pivot=next(i for i in range(col,n) if a[i][col])
        if pivot!=col:a[pivot],a[col]=a[col],a[pivot];det=-det
        v=a[col][col];det*=v;a[col]=[x/v for x in a[col]]
        for i in range(n):
            if i==col:continue
            v=a[i][col]
            if v:a[i]=[x-v*y for x,y in zip(a[i],a[col])]
    assert det.denominator==1
    return [a[i][-1] for i in range(n)],det.numerator

def signed_multiset_sequences(cols):
    types=sorted(Counter(cols));counts=[cols.count(c) for c in types]
    def rec(prefix):
        if len(prefix)==6:
            yield tuple(prefix);return
        for j,c in enumerate(types):
            if counts[j]:
                counts[j]-=1
                for s in (-1,1):yield from rec(prefix+[tuple(s*x for x in c)])
                counts[j]+=1
    return rec([])
try:
    res=json.loads((source/'results.json').read_text());expected={r['id']:r for r in res['cases']}
    saved_raw=gzip.decompress((source/'selected_terms.json.gz').read_bytes())
    assert hashlib.sha256(saved_raw).hexdigest()==res['uncompressed_terms_sha256']
    saved={x['case']:x['terms'] for x in json.loads(saved_raw)['cases']}
    words=iter(map(int,FAN.read_text().split()));assert next(words)==21
    geom=json.loads(GEOM.read_text());cones=[];vol=F(0);signwords=set()
    body=[(1,)+v for v in product((-1,1),repeat=5)]
    for ci in range(21):
        de=next(words);qs=[];rays=[]
        for j in range(6):
            qs.append(next(words));rays.append([next(words) for _ in range(6)])
        R=fmpz_mat([list(t) for t in zip(*rays)])
        assert abs(int(R.det()))==de>0
        normalized=[]
        for v,q in zip(rays,qs):
            assert v[0]>0 and v[-1]>=0 and all(v[i]>=v[i+1] for i in range(5))
            assert sum(abs(t) for t in v)<=12
            assert q==sum(abs(sum(a*b for a,b in zip(v,s))) for s in body)>0
            normalized.append([F(t,v[0]) for t in v])
        old=geom['cells'][ci]
        assert normalized==[list(map(F,v)) for v in old['rays_a0_equals_one']]
        assert [F(q,v[0]) for v,q in zip(rays,qs)]==list(map(F,old['rates']))
        # Section simplex is contained in the ordered chamber. Constant signs
        # distinguish interiors, while exact total volume gives completeness.
        center=[sum(v[j] for v in normalized)/6 for j in range(6)]
        signs=[]
        for b in body:
            vals=[sum(a*t for a,t in zip(b,v)) for v in rays]
            assert min(vals)>=0 or max(vals)<=0
            z=sum(F(a)*t for a,t in zip(b,center));assert z
            signs.append(1 if z>0 else -1)
        assert tuple(signs) not in signwords;signwords.add(tuple(signs))
        vol+=F(de,prod(v[0] for v in rays)*factorial(5))
        cones.append((de,qs,rays,R))
    assert next(words,None) is None and vol==F(1,120)
    assert sha(FAN)=='7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812'
    # Independent combinatorial orbit counts, not replay of 4960 triples.
    classes=[(1,1,1),(1,1,-1),(1,-1,1),(1,-1,-1)]
    actions=set()
    for p in permutations(range(3)):
        for s in product((-1,1),repeat=3):
            ims=[]
            for v in classes:
                w=[s[i]*v[p[i]] for i in range(3)];ims.append(classes.index(tuple(x*w[0] for x in w)))
            actions.add(tuple(ims))
    assert len(actions)==24
    parts=[p for p in product(range(7),repeat=4) if sum(p)==6 and p[0]>=p[1]>=p[2]>=p[3] and p[2]>0]
    assert len(parts)==5
    orbit={}
    for p in parts:
        labeled=0
        for n in set(permutations(p)):
            if n[0]:labeled+=factorial(5)//(factorial(n[0]-1)*prod(factorial(x) for x in n[1:]))
        assert (32*labeled)%6==0
        orbit[p]=32*labeled//6
    assert sum(orbit.values())==4960
    pair_counts={1:96,2:240,3:160};assert sum(pair_counts.values())==496
    configs=[('r1',[(1,)]*6,None)]+[(f'r2_k{k}',[(1,1)]*(6-k)+[(1,-1)]*k,None) for k in (1,2,3)]
    configs += [('r3_'+''.join(map(str,p)),[c for c,n in zip(classes,p) for _ in range(n)],orbit[p]) for p in sorted(parts,reverse=True)]
    assert set(x[0] for x in configs)==set(expected)==set(saved)
    total_unique=0;total_saved=0;answer=[]
    for name,cols,orbit_count in configs:
        old=expected[name];r=len(cols[0]);d=6-r;v=[1009**(r-i-1) for i in range(r)]
        L=[list(x) for x in zip(*cols)]
        gram=fdet([[sum(a*b for a,b in zip(x,y)) for y in L] for x in L])
        assert gram>0 and L==old['normal_rows'] and gram==old['Gram_determinant']
        assert old['trace_direction']==v and old['dimension']==d and old['r']==r
        mult=prod(factorial(n) for n in Counter(cols).values())
        sequences=list(signed_multiset_sequences(cols));assert len(sequences)==64*720//mult
        assert len(set(sequences))==len(sequences)
        patterns=[(tuple(zip(*pc)),fmpz_mat([list(t) for t in zip(*pc)])) for pc in sequences]
        unsigned=sorted({tuple(tuple(t*c[0] for t in c) for c in pc) for pc in sequences})
        assert len(unsigned)==720//mult==old['unique_coordinate_permutations']
        assert mult==old['permutation_multiplicity'] and orbit_count==old['triple_orbit_count']
        if orbit_count is not None:assert orbit_count=={64:480,96:1920,128:480,160:640,192:1440}[gram]
        hist=Counter();positives={};trace=F(0);zero_solutions=0
        # Cone-major traversal; every projected matrix is formed independently
        # by exact FLINT multiplication rather than the producer's row table.
        for ci,(de,qs,rays,R) in enumerate(cones):
            for normal_rows,mat in patterns:
                Z=mat*R
                rr=[[int(Z[i,j]) for j in range(6)] for i in range(r)]
                assert all(min(row)>=0 or max(row)<=0 for row in rr)
                assert all(abs(x)<=12 for row in rr for x in row)
                J=[j for j in range(6) if any(rr[i][j] for i in range(r))]
                m=len(J);assert m>=r;hist[m]+=mult;total_unique+=1
                if m!=r:continue
                M=[[rr[i][j] for j in J] for i in range(r)]
                sol,dm=solve_and_det(M,v)
                assert all(sol),'generic ray must have no solved zero'
                if not all(x>0 for x in sol):continue
                w=F(de*mult,abs(dm)*prod(qs[j] for j in range(6) if j not in J))
                key=(normal_rows,ci);assert key not in positives
                positives[key]=(J,dm,sol,w);trace+=w
        assert sum(hist.values())==967680
        assert dict(hist)=={int(k):v for k,v in old['m_histogram_logical'].items()}
        assert trace==F(old['R_exact']) and len(positives)==old['selected_terms']
        # All 260 saved entries: reconstruct identifiers, compare determinant,
        # Cramer metadata via solved rational coordinates, and exact weight.
        seen=set()
        for term in saved[name]:
            pc=unsigned[term['permutation_index']];mask=term['sign_mask'];ci=term['cone_index']
            rows=tuple(tuple(pc[k][i]*(1 if mask&(1<<k) else -1) for k in range(6)) for i in range(r))
            key=(rows,ci);assert key not in seen;seen.add(key)
            J,dm,sol,w=positives[key]
            assert J==term['nonzero_normal_columns'] and dm==term['detM']
            assert [x*dm for x in sol]==list(map(F,term['cramer_numerators']))
            assert w==F(term['weight']);total_saved+=1
        assert seen==set(positives)
        # Generic-ray proof: all nonzero rank2 normals are bounded by12;
        # all rank3 cross minors by288. Base1009 cannot cancel.
        for n in old['generic_facets']:
            assert len(n)==r and any(n) and max(map(abs,n))<=288
            assert sum(a*b for a,b in zip(n,v))!=0
        Rq=arb(str(trace.numerator))/trace.denominator
        kap=Rq/(2*arb.pi())**d;peak=arb(gram).sqrt()*kap
        assert kap.overlaps(arb(old['kappa_arb512'])) and peak.overlaps(arb(old['peak_arb512']))
        answer.append({'id':name,'r':r,'d':d,'Gram':gram,'R_exact':str(trace),
                       'kappa_exact':f'({trace})/(2*pi)^{d}','peak_exact':f'sqrt({gram})*({trace})/(2*pi)^{d}',
                       'unique_patterns':len(patterns)*21,'logical_items':sum(hist.values()),
                       'm_histogram_logical':dict(hist),'saved_selected_terms':len(seen),
                       'multiplicity':mult,'triple_orbit_count':orbit_count,
                       'peak_arb512':str(peak),'peak_lower_exact':str(peak.lower().fmpq()),
                       'peak_upper_exact':str(peak.upper().fmpq())})
        receipt['completed_pattern_items']=total_unique;receipt['completed_cases']=[a['id'] for a in answer];write('execution_receipt.json',receipt)
        print(name,'PASS',len(positives),'selected',str(trace),flush=True)
    assert total_unique==701568 and total_saved==260 and sum(a['logical_items'] for a in answer)==8709120
    result={'status':'PASS_INDEPENDENT_FULL_COVERAGE_AND_ALL_260_SAVED_TERMS',
            'cases':answer,'unique_pattern_cone_items':total_unique,'logical_items':8709120,
            'selected_items_checked':total_saved,'single_full_enumeration':True,
            'fan_cover':{'section_volume':'1/120','contained_in_chamber':True,'pairwise_distinct_interior_signwords':21,'all_body_forms_constant_weak_sign':True},
            'projective_S4_actions':24,'triple_orbits':{str(p):n for p,n in orbit.items()},
            'generic_base':1009,'normal_entry_bound':12,'rank3_minor_bound':288,
            'input_sha256':inputs,'source_sha256':sha(Path(__file__)),
            'native_G_H_calls':0,'spatial_grid':0,'optimization':0,'matlab':0,
            'pid':os.getpid(),'elapsed_seconds':time.time()-start,
            'elapsed_since_registration':time.time()-registration['started_unix']}
    write('results.json',result);receipt['status']='PASS_STOPPED'
except BaseException as e:
    receipt['status']='FAILED_STOPPED_NO_RETRY';receipt['error']=repr(e);raise
finally:
    receipt['elapsed_seconds']=time.time()-start;write('execution_receipt.json',receipt);signal.alarm(0)
print(json.dumps({k:result[k] for k in ('status','unique_pattern_cone_items','logical_items','selected_items_checked','elapsed_seconds')}))
