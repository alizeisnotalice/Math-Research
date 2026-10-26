"""Frozen056 only: integer multibasis81 requalification, no author functions."""
import argparse,hashlib,itertools,json,math,os,sys,time
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent;RUN='audit_multibasis060_integer_001'
MANIFEST='0f3261b09dabc9557903547c0670125392e2bf106a8f4e7c3d952005b3978ba4'
STATE='b79c134fdcec283c120d464cf78f131c669c905ac5cdacda309a34f07a9fb29b'
STATE_REG='585a5de78b8f159ca03f138013894041dba350d144c3217b7bbda09fd5c7c9e9'
SOURCE046='df3fb4165aea29e4aa3d1fb7ecacdb9cc5ad4a531bc03e3cc43a0820671fe028'
SOURCE043='a298ca56ba980fce4fe2542929e03befb97b528afd3d390ed9d15651da68322a'
SOURCE047='0b680c1bd7e4aea30c027fe32cc3e51436f7dda6dd8b74ddf521c2d74bb8ca16'
CHECKPOINT=Path('/tmp/audit_multibasis060_checkpoint.json')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(o):return json.dumps(o,sort_keys=True,separators=(',',':')).encode()
def digest(o):return hashlib.sha256(canon(o)).hexdigest()
def project(p):
    p=p.resolve();return p/'cube-dim-5-bounds' if (p/'cube-dim-5-bounds').is_dir() else p
def save(p,o):p.write_text(json.dumps(o,indent=2)+'\n')
def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def det(a):
    return sum((-1)**sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))*
               math.prod(a[i][p[i]] for i in range(5)) for p in itertools.permutations(range(5)))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--author-root',type=Path,required=True);ap.add_argument('--upper-root',type=Path,required=True)
    ap.add_argument('--state-root',type=Path,required=True)
    ap.add_argument('--run',action='store_true');ap.add_argument('--arm',action='store_true');args=ap.parse_args()
    author=project(args.author_root)/'upper_multibasis056';upper=project(args.upper_root)
    manifest=author/'manifest.json';assert sha(manifest)==MANIFEST
    locked={'upper_multibasis056/manifest.json':MANIFEST}
    m=json.loads(manifest.read_text())
    for row in m['white_list']:
        p=author/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes']
        locked['upper_multibasis056/'+row['path']]=row['sha256']
    reg=json.loads((author/'registration.json').read_text())
    sp=args.state_root.resolve()/'state.json';sr=args.state_root.resolve()/'registration.json'
    sources=[(sp,STATE),(sr,STATE_REG),(upper/'upper_cover046/arb_cover.py',SOURCE046),
             (upper/'upper_harnack043/report.md',SOURCE043),
             (project(args.author_root)/'audit_harnack047/report.md',SOURCE047)]
    for p,h in sources:assert sha(p)==h;locked['/'.join(p.parts[-2:])]=h
    assert reg['input_state_sha256']==STATE and reg['input_registration_sha256']==STATE_REG
    assert reg['source046_sha256']==SOURCE046
    if not args.run:
        print(json.dumps({'run_handle':RUN,'status':'locked_inputs_only','calculation_executed':False}));return
    pid=os.getpid();print(json.dumps({'run_handle':RUN,'pid':pid,'status':'armed' if args.arm else 'started'}),flush=True)
    if args.arm:input()
    began=time.monotonic();rows=[]
    def checkpoint(stage):
        elapsed=time.monotonic()-began
        if len(rows)%16==0 or elapsed>120:
            save(CHECKPOINT,{'run_handle':RUN,'pid':pid,'stage':stage,'completed_boxes':len(rows),
                 'rows':rows,'locked_inputs_sha256':locked,'code_sha256':sha(Path(__file__)),'seconds':elapsed})
        if elapsed>120:raise TimeoutError('finite162 checkpoint; no automatic restart')
    # Independent integer inverse/Gram path; no symbolic matrix library.
    columns=[(1,1,1,1,1),(1,-1,-1,-1,1),(1,-1,-1,1,-1),
             (1,-1,1,-1,-1),(1,1,-1,-1,-1)]
    N=[list(row) for row in zip(*columns)]
    M=[[2,1,1,1,1],[1,-1,-1,-1,2],[1,-1,-1,2,-1],
       [1,-1,2,-1,-1],[1,2,-1,-1,-1]]
    sixI=[[6*int(i==j) for j in range(5)] for i in range(5)]
    assert det(N)==-48 and mm(N,M)==sixI and mm(M,N)==sixI
    signs=list(itertools.product((-1,1),repeat=4));inverses=[]
    for s in signs:
        d=(1,*s);Ns=[[d[i]*N[i][j] for j in range(5)] for i in range(5)]
        Ms=[[M[i][j]*d[j] for j in range(5)] for i in range(5)]
        assert mm(Ns,Ms)==sixI and mm(Ms,Ns)==sixI and abs(det(Ns))==48
        assert all(col[0]==1 and set(col)<=set((-1,1)) for col in zip(*Ns))
        T=[row[1:] for row in Ms]
        gram=mm([list(row) for row in zip(*T)],T)
        assert gram==[[9*int(i==j)-s[i]*s[j] for j in range(4)] for i in range(4)]
        inverses.append(T)
    # Nonconcavity is a concrete algebraic warning, not a hypothesis used below.
    def psi(v):return sum(x*x for x in v)/F(4)-sum(map(abs,v))**2/F(36)
    gap=-psi((F(0),F(5),F(0),F(0)))+psi((F(1),F(5),F(0),F(0)))
    assert gap==F(-1,18)
    state=json.loads(sp.read_text());oldreg=json.loads(sr.read_text())
    assert state['fingerprint']==reg['input_fingerprint']
    assert state['source052_hashes']==oldreg['sources']==reg['source052_hashes']
    assert state['source_hashes']['arb_cover.py']==SOURCE046
    assert state['epsilon']==oldreg['epsilon']=='1/64' and state['target']==oldreg['target']=='1/10000'
    assert state['batch_number']==2 and len(state['splits'])==162 and not state['inflight'] and oldreg['parent_splits']==0
    assert digest(state['roots'])==state['initial_root_sha256']
    data=json.loads((author/'results.json').read_text());reference={tuple(r['node']):r for r in data['certificates']}
    assert len(reference)==162
    seen=set();allpoints=0;allbasis=0;allold=0;allstrict=0
    for split in state['splits']:
        node=tuple(split['node']);assert node not in seen and node[1]=='';seen.add(node)
        cert=split['failed_certificate'];assert cert['passed'] is False
        root=state['roots'][node[0]];assert len(root)==4 and all(isinstance(i,int) and 0<=i<32 for i in root)
        cached=state['cache'][cert['anchor']];point=tuple(map(F,cached['point']))
        assert point[0]==1 and ','.join(map(str,point))==cert['anchor']
        lo=[F(i,32) for i in root];hi=[v+F(1,32) for v in lo]
        assert all(l<=a<=h for l,a,h in zip(lo,point[1:],hi))
        allstrict+=all(l<a<h for l,a,h in zip(lo,point[1:],hi))
        values=[tuple(map(F,r)) for r in cached['values']]
        assert len(values)==5 and all(l<=h for l,h in values) and F(1,10000)<values[0][0]
        gradient=[tuple(map(F,r)) for r in cert['loggradient_directed']]
        assert len(gradient)==4 and all(l<=h for l,h in gradient)
        # Positive denominator cross-products enclose every endpoint ratio.
        for (derlo,derhi),(gl,gh) in zip(values[1:],gradient):
            for h in values[0]:
                for deriv in (derlo,derhi):assert gl*h<=deriv<=gh*h
        Delta=F(cert['Delta_lower']);assert Delta>0
        lows=[l-a for l,a in zip(lo,point[1:])];highs=[h-a for h,a in zip(hi,point[1:])]
        vd=math.lcm(*(v.denominator for v in (*lows,*highs)))
        lints=[int(v*vd) for v in lows];hints=[int(v*vd) for v in highs]
        assert all(l<=0<=h for l,h in zip(lints,hints))
        # Clear all real arithmetic once; the 81*16 inner loop is integers only.
        cd=math.lcm(Delta.denominator,36*vd*vd,*(g.denominator*vd for pair in gradient for g in pair))
        deltaint=int(Delta*cd);gints=[]
        assert F(4096*cd,36*vd*vd).denominator==1
        qfactor=(4096*cd)//(36*vd*vd)
        for gl,gh in gradient:
            pair=gl*cd/vd,gh*cd/vd;assert all(g.denominator==1 for g in pair);gints.append(tuple(map(int,pair)))
        def score(v,T=None):
            lin=sum(min(gl*x,gh*x) for x,(gl,gh) in zip(v,gints))
            if T is not None:
                curvature=sum(sum(a*x for a,x in zip(row,v))**2 for row in T)
            else:
                curvature=min(sum(sum(a*x for a,x in zip(row,v))**2 for row in T) for T in inverses)
                assert curvature==9*sum(x*x for x in v)-sum(map(abs,v))**2
            return deltaint+lin-qfactor*curvature
        baseT=[row[1:] for row in M];old=[]
        for bits in itertools.product((0,1),repeat=4):
            v=tuple((hints[j] if b else lints[j]) for j,b in enumerate(bits));old.append((score(v,baseT),bits))
        oldmin,oldbits=min(old);assert F(oldmin,cd)==F(cert['minimum_log_margin'])<=0
        assert list(oldbits)==cert['worst_corner_bits'];allold+=16
        # Generate vertices via closed orthant-cell intersections, not a direct ternary product.
        candidates=set()
        for s in signs:
            sides=[(0,h) if sig==1 else (l,0) for sig,l,h in zip(s,lints,hints)]
            for v in itertools.product(*sides):candidates.add(v)
        assert candidates==set(itertools.product(*(set((l,0,h)) for l,h in zip(lints,hints))))
        scores=[(score(v),v) for v in candidates];newmin,worst=min(scores)
        assert newmin>0 and newmin>=oldmin
        ref=reference[node]
        assert ref['root']==root and ref['anchor']==cached['point'] and ref['Delta_lower']==cert['Delta_lower']
        assert ref['loggradient_directed']==cert['loggradient_directed']
        assert F(ref['new_minimum'])==F(newmin,cd) and F(ref['old_16_minimum'])==F(oldmin,cd)
        assert ref['new_point_count']==len(candidates) and ref['new_certificate_passed'] is True
        assert tuple(map(F,ref['new_worst_offset']))==tuple(F(v,vd) for v in worst)
        assert F(ref['exact_improvement'])==F(newmin-oldmin,cd)
        allpoints+=len(candidates);allbasis+=len(candidates)*16
        rows.append({'node':list(node),'root':root,'candidate_count':len(candidates),
             'old_minimum':str(F(oldmin,cd)),'new_minimum':str(F(newmin,cd)),
             'worst_offset':[str(F(v,vd)) for v in worst],
             'improvement':str(F(newmin-oldmin,cd)),'new_certificate_passed':True,
             'offset_common_denominator_bits':vd.bit_length(),'margin_common_denominator_bits':cd.bit_length()})
        checkpoint('integer_162_prefix')
    assert seen==set(reference) and data['summary']['newly_certified_parents']==len(rows)==162
    assert sha(sp)==STATE and sha(sr)==STATE_REG
    minimum=min(F(r['new_minimum']) for r in rows);improvement=min(F(r['improvement']) for r in rows)
    assert minimum==F(data['summary']['minimum_new_log_margin'])
    assert improvement==F(data['summary']['minimum_exact_improvement'])
    out={'run_handle':RUN,'pid':pid,'status':'PASS_056_selected162_closed_parent_boxes_integer_multibasis81',
         'algebra':{'signed_legal_bases':16,'inverse_denominator':6,'tail_Gram_denominator':36,
                    'tail_Gram_integer':'9I-ssT','positive_eigenvalues':['1/4','5/36'],
                    'global_negative_curvature_midpoint_concavity_gap':str(gap)},
         'summary':{'closed_parent_boxes':len(rows),'passed':162,'old_corners_checked':allold,
                    'new_vertices_checked':allpoints,'integer_signed_basis_curvatures':allbasis,
                    'strict_inside_anchors':allstrict,'all_gradient_quotient_endpoint_crossproducts_valid':True,
                    'all81_minima_and_worst_offsets_equal_author':True,'minimum_new_margin':str(minimum),
                    'minimum_improvement':str(improvement),'new_H_or_gradient_calls':0,
                    'new_log_exp_Gamma_Arb_calls':0,'source_state_unchanged':True},
         'certificates':rows,'locked_inputs_sha256':locked,'source_sha256':sha(Path(__file__)),
         'registration_sha256':sha(HERE/'registration.md'),
         'dependencies':'Saved H/derivative balls and Delta lower bounds are accepted locked numerical inputs after source read-audit, not independently recomputed. Positive-mixture theory is proved in report.',
         'scope':'Exactly162 saved split parents only. No052 controller/cache/queue mutation, no global cover/g/C5 conclusion.',
         'runtime':{'python':sys.executable,'seconds':time.monotonic()-began,'author_functions_called':0,
                    'matrix_or_interval_libraries_imported':0}}
    save(HERE/'result.json',out)
    save(CHECKPOINT,{'run_handle':RUN,'pid':pid,'stage':'complete','completed_boxes':len(rows),'seconds':time.monotonic()-began})
    print(json.dumps({k:out[k] for k in ('run_handle','pid','status','algebra','summary','runtime')},indent=2))

if __name__=='__main__':main()
