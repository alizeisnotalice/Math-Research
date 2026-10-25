"""Finite pair-orbit and all-weight budget audit; default hashes only."""
import argparse,hashlib,json,math,sys,time
from collections import Counter
from fractions import Fraction as F
from itertools import permutations,product
from pathlib import Path
from flint import arb,ctx

HERE=Path(__file__).resolve().parent;RUN='audit_I026_fixed_001'
LOCKS={'upper_I022':'668f4b36fa225c865da727a1bf92f211eccfb06a35aad0ebc44cfa7071531602',
       'upper_I024':'96f65eaa792afc290d3b83aa2284f64170386de045e54f04070b96335f93b764'}
TRANSFER='3d2acf640e1be0ee11d10fd6e29719ddc4291117975407e0175c8142573d5465'
V=[(1,)+v for v in product((-1,1),repeat=4)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def project(p):
    p=p.resolve();return p/'cube-dim-5-bounds' if (p/'cube-dim-5-bounds').is_dir() else p
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def det(a):
    n=len(a)
    return sum((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))*math.prod(a[i][p[i]] for i in range(n)) for p in permutations(range(n)))
def null(rows):
    a=[list(map(F,r)) for r in rows];m=len(a);pivots=[];r=0
    for c in range(5):
        pivot=next((i for i in range(r,m) if a[i][c]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
        for i in range(m):
            if i!=r:
                v=a[i][c];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        pivots.append(c);r+=1
        if r==m:break
    assert r==m;out=[]
    for c in sorted(set(range(5))-set(pivots)):
        z=[F(0)]*5;z[c]=1
        for i,p in enumerate(pivots):z[p]=-a[i][c]
        out.append(z)
    return out
def coordinates(rows):
    N=null(rows);H=[[dot(a,b) for b in N] for a in N]
    return [tuple(dot(n,v) for n in N) for v in V],det(H)
def classes(z):
    out={}
    for i,v in enumerate(z):
        if not any(v):continue
        axis=next(j for j,x in enumerate(v) if x);key=tuple(x/v[axis] for x in v)
        out.setdefault(key,[]).append(i)
    return [tuple(row) for row in out.values()]
def hist(values):
    total=sum(values);assert total>0
    return tuple(sorted(Counter(x/total for x in values if x).items()))
def generators():
    index={v:i for i,v in enumerate(V)};out=[]
    for kind in range(9):
        mapping=[]
        for v in V:
            w=list(v)
            if kind<5:w[kind]*=-1
            else:
                c=kind-5;w[c],w[c+1]=w[c+1],w[c]
            mapping.append(index[tuple(x*w[0] for x in w)])
        out.append(mapping)
    return out
GEN=generators()
def orbit(state):
    seen={state};todo=[state]
    while todo:
        x=todo.pop()
        for g in GEN:
            image=tuple(tuple(sorted(g[i] for i in part)) for part in x)
            if image not in seen:seen.add(image);todo.append(image)
    return seen
def pair_geometry(J):
    rows=[V[i] for i in J];D=det([[dot(a,b) for b in rows] for a in rows]);assert D>0
    z,H3=coordinates(rows);assert H3>0 and {i for i,v in enumerate(z) if not any(v)}==set(J)
    Xclasses=classes(z);groups=[];spot=False
    for X in Xclasses:
        r=X[0];axis=next(i for i,x in enumerate(z[r]) if x)
        ax={i:(z[i][axis]/z[r][axis])**2 for i in X}
        rows3=rows+[V[r]];D3=det([[dot(a,b) for b in rows3] for a in rows3]);assert D3>0
        residual,H2=coordinates(rows3);assert H2>0
        assert {i for i,v in enumerate(residual) if not any(v)}==set(J)|set(X)
        for Y in classes(residual):
            s=Y[0];axis=next(i for i,x in enumerate(residual[s]) if x)
            ay={j:(residual[j][axis]/residual[s][axis])**2 for j in Y}
            outside=tuple(i for i in range(16) if i not in set(J)|set(X)|set(Y))
            assert len(J)+len(X)+len(Y)+len(outside)==16
            weights={}
            for k in outside:
                area=residual[s][0]*residual[k][1]-residual[s][1]*residual[k][0]
                weights[k]=F(D3,H2)*area**2
                assert weights[k]>0 and weights[k].denominator==1
            for i in X:
                for j in Y:
                    for k in outside:
                        full=F(D,H3)*det([z[i],z[j],z[k]])**2
                        assert full==ax[i]*ay[j]*weights[k]
            if not spot:
                i,j,k=X[0],Y[0],outside[0]
                assert det(rows+[V[i],V[j],V[k]])**2==ax[i]*ay[j]*weights[k];spot=True
            C=sum(ax.values())*sum(ay.values())*sum(weights.values());assert C>0 and C.denominator==1
            key=tuple(sorted((hist(list(ax.values())),hist(list(ay.values())),hist(list(weights.values())))))
            groups.append((X,Y,int(C),key))
    assert sum(g[2] for g in groups)==6*D*16**3
    return D,[len(X) for X in Xclasses],groups
def ar(x):
    x=F(x);return arb(x.numerator)/x.denominator
def parse(s):
    a,b=map(F,s[1:-1].split(' +/- '));return a-b,a+b
def bounds(x):return F(str(x.lower().fmpq())),F(str(x.upper().fmpq()))
def iroot(n,k):
    x=1<<((n.bit_length()+k-1)//k)
    while True:
        y=((k-1)*x+n//x**(k-1))//k
        if y>=x:break
        x=y
    while x**k>n:x-=1
    while (x+1)**k<=n:x+=1
    return x
def rational_root_ball(ball,k):
    lo,hi=bounds(ball);scale=10**110
    a=F(iroot(lo.numerator*scale**k//lo.denominator,k),scale)
    b=F(iroot(hi.numerator*scale**k//hi.denominator,k)+1,scale)
    assert a**k<=lo<hi<=b**k;return a,b
def encoded_bounds(row):return F(row['lower_fraction']),F(row['upper_fraction'])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',required=True,type=Path);parser.add_argument('--model-root',required=True,type=Path)
    parser.add_argument('--run',action='store_true');args=parser.parse_args()
    source,model=project(args.source_root),project(args.model_root);locked={};manifests={}
    for directory,expected in LOCKS.items():
        p=source/directory/'manifest.json';assert sha(p)==expected;locked[directory+'/manifest.json']=expected
        manifests[directory]=json.loads(p.read_text())
        for row in manifests[directory]['white_list']:
            p=source/directory/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes']
            locked[directory+'/'+row['path']]=sha(p)
        for row in manifests[directory]['source_inputs']:
            suffix=Path(row['source']).parts[-2:]
            options=[root.joinpath(*suffix) for root in (source,model)]
            p=next((p for p in options if p.exists() and sha(p)==row['sha256']),None);assert p is not None
            locked['/'.join(suffix)]=sha(p)
    assert sha(model/'alternative_alpha/transfer.md')==TRANSFER;locked['alternative_alpha/transfer.md']=TRANSFER
    old=json.loads((source/'upper_I022/results.json').read_text());opt=json.loads((source/'upper_I024/results.json').read_text())
    if not args.run:
        print(json.dumps({'run_handle':RUN,'status':'locked_inputs_only','calculation_executed':False}));return
    began=time.monotonic();pairs=set();covered={};types={};orbits=[];pair_receipts=[]
    for J in ((0,1),(0,3)):
        D,sizes,groups=pair_geometry(J);images=orbit((J,));assert not(pairs&images);pairs|=images
        pair_receipts.append({'J':J,'D':D,'orbit_size':len(images),'first_line_sizes':sorted(sizes)})
        for X,Y,C,key in groups:
            state=(J,X,Y)
            if state in covered:assert covered[state]==(C,key);continue
            flags=orbit(state);assert not(set(flags)&set(covered))
            for state2 in flags:covered[state2]=(C,key)
            row=types.setdefault(key,{'count':0,'sum_C':0});row['count']+=len(flags);row['sum_C']+=len(flags)*C
            orbits.append({'J':J,'X':X,'Y':Y,'orbit_size':len(flags),'C':C})
    mass=60*16**5;assert len(pairs)==120 and len(covered)==4560 and sum(r['sum_C'] for r in types.values())==mass
    expected={}
    for row in old['geometry']['probability_types']:
        key=tuple(sorted(tuple(sorted((F(p),n) for p,n in h.items())) for h in row['probability_histograms']))
        expected[key]={'count':row['count'],'sum_C':row['sum_C']};assert F(row['theta'])==F(row['sum_C'],mass)
    assert expected==types and len(types)==5
    ctx.prec=512;beta,b=F(3,4),F(3,20);sqrt2,sqrt3=arb(2).sqrt(),arb(3).sqrt()
    R=[(7+8*sqrt2+3*sqrt3)/arb(2)**ar(F(15,4)),arb(3)**ar(F(1,4))*(1+1/sqrt2),1+sqrt2,1+sqrt2,2*sqrt2]
    theta=[F(1,8),F(3,16),F(3,16),F(1,8),F(3,8)]
    for row,r,t in zip(old['geometry']['probability_types'],R,theta):
        assert F(row['theta'])==t;raw=arb(1)
        for h in row['probability_histograms']:raw*=sum((n*ar(p)**ar(beta) for p,n in h.items()),arb(0))
        assert (raw-r).contains(0)
    Z=sum((ar(t)*r**ar(1/beta) for r,t in zip(R,theta)),arb(0))
    short=ar(F(3,2))+(5+3*ar(F(3,4))**ar(F(1,3)))/16*(1+sqrt2)**ar(F(4,3))+(7+8*sqrt2+3*sqrt3)**ar(F(4,3))/256
    assert (Z-short).contains(0)
    Zlo,Zhi=bounds(Z);savedZlo,savedZhi=encoded_bounds(opt['Z_exact_interval']);assert savedZlo<=Zlo<Zhi<=savedZhi
    m=ar(b/beta).gamma()/(ar(beta)*ar(b).gamma());gamma7=ar(beta)*ar(beta).gamma()
    N=gamma7**10*m**50;C7=ar(F(3**10,2**65*60**7))/arb.pi()**11
    gain10=arb(1)
    for r,t in zip(R,theta):gain10*=r**ar(-2*t)
    I22=rational_root_ball(C7*N*gain10,7);J22=rational_root_ball(2**10*N*gain10,10)
    Iopt=rational_root_ball(C7*N/Z**ar(F(3,2)),7);Jopt=rational_root_ball(2**10*N/Z**ar(F(3,2)),10)
    for actual,saved in ((I22,parse(old['arb_384bit']['complete_I_upper'])),(J22,parse(old['arb_384bit']['J_upper'])),
                         (Iopt,encoded_bounds(opt['I_optimal_budget_interval'])),(Jopt,encoded_bounds(opt['J_optimal_budget_interval']))):
        assert saved[0]<=actual[0]<actual[1]<=saved[1]
    assert I22[1]<F(old['strict_rational_I_upper']) and Iopt[1]<F(opt['strict_rational_I_upper'])<I22[0]
    assert Iopt[0]>F(118,10**7)>F(11,10**6)
    qstar=[ar(t)*r**ar(F(4,3))/Z for t,r in zip((F(1,8),F(3,16),F(5,16),F(3,8)),(R[0],R[1],R[2],R[4]))]
    assert (sum(qstar,arb(0))-1).contains(0)
    for q,row in zip(qstar,opt['merged_type_qopt_intervals']):
        lo,hi=bounds(q);a,c=encoded_bounds(row['qopt_mass']);assert a<=lo<=hi<=c
    out={'run_handle':RUN,'status':'pass_I022_complete_bound_and_I024_all_flag_budget_optimum',
         'geometry':{'pair_orbits':pair_receipts,'flag_orbits':orbits,'pair_count':len(pairs),'flag_count':len(covered),
                      'probability_types':len(types),'mass_60W':mass,'scope_full_16_for_every_flag':True},
         'I022_formula_interval':[str(x) for x in I22],'J022_formula_interval':[str(x) for x in J22],
         'I022_strict_upper':old['strict_rational_I_upper'],
         'I024_Z_ball_512bit':str(Z),'I024_optimal_budget_I_interval':[str(x) for x in Iopt],
         'I024_optimal_budget_J_interval':[str(x) for x in Jopt],'I024_strict_upper':opt['strict_rational_I_upper'],
         'family_optimum_not_actual_I_lower_bound':True,'all_4560_weights_including_boundary_covered_by_KL':True,
         'k0_k1_empty_layers_duplicate_unordered_pair_flags':True,
         'shared_library_limit':'Arb/Gamma shared; geometry is independent RREF/nullspace/orbits, Mellin arguments unshifted, I^7/J^10 and Fraction roots differ from production. No cross-library Gamma certificate claimed.',
         'locked_inputs_sha256':locked,'source_sha256':sha(Path(__file__)),
         'registration_sha256':sha(HERE/'registration.md'),
         'runtime':{'python':sys.executable,'Arb_bits':ctx.prec,'seconds':time.monotonic()-began,'author_120_pair_loop_replayed':False}}
    (HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('run_handle','status','I022_strict_upper','I024_strict_upper','runtime')},indent=2))

if __name__=='__main__':main()
