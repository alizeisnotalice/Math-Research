"""One orbit/nullspace and seventh-power audit; default hashes only."""
import argparse,hashlib,json,math,sys,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,permutations,product
from pathlib import Path
from flint import arb,ctx

HERE=Path(__file__).resolve().parent;RUN='audit_I023_fixed_001'
MANIFEST='16eeeb3e07fc15813c7ece29d415ff79c219fa116b6a5add2f4926860c1ad42b'
TRANSFER='3d2acf640e1be0ee11d10fd6e29719ddc4291117975407e0175c8142573d5465'
V=[(1,)+v for v in product((-1,1),repeat=4)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def proj(p):
    p=p.resolve();return p/'cube-dim-5-bounds' if (p/'cube-dim-5-bounds').is_dir() else p
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def determinant(a):
    n=len(a)
    return sum((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))*math.prod(a[i][p[i]] for i in range(n)) for p in permutations(range(n)))
def nullspace(rows):
    a=[list(map(F,r)) for r in rows];pivots=[];r=0
    for col in range(5):
        pivot=next((i for i in range(r,3) if a[i][col]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];c=a[r][col];a[r]=[x/c for x in a[r]]
        for i in range(3):
            if i!=r:
                c=a[i][col];a[i]=[x-c*y for x,y in zip(a[i],a[r])]
        pivots.append(col);r+=1
        if r==3:break
    assert r==3
    basis=[]
    for col in set(range(5))-set(pivots):
        z=[F(0)]*5;z[col]=1
        for i,p in enumerate(pivots):z[p]=-a[i][col]
        basis.append(z)
    assert len(basis)==2
    return basis
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
                i=kind-5;w[i],w[i+1]=w[i+1],w[i]
            mapping.append(index[tuple(x*w[0] for x in w)])
        out.append(mapping)
    return out
GEN=generators()
def orbit(state):
    seen={state};todo=[state]
    while todo:
        current=todo.pop()
        for g in GEN:
            image=tuple(tuple(sorted(g[i] for i in part)) for part in current)
            if image not in seen:seen.add(image);todo.append(image)
    return seen
def representative(J):
    rows=[V[i] for i in J];D=determinant([[dot(a,b) for b in rows] for a in rows]);assert D>0
    N=nullspace(rows);H=[[dot(a,b) for b in N] for a in N];detH=determinant(H);assert detH>0
    z=[tuple(dot(n,v) for n in N) for v in V]
    nonzero=tuple(i for i,p in enumerate(z) if any(p));assert not(set(J)&set(nonzero))
    classes={}
    for i in nonzero:
        x,y=z[i];key=(F(1),y/x) if x else (F(0),F(1))
        classes.setdefault(key,[]).append(i)
    w=[[F(0)]*16 for _ in range(16)]
    for i,k in combinations(range(16),2):
        area=z[i][0]*z[k][1]-z[i][1]*z[k][0]
        value=F(D,detH)*area**2;assert value>=0 and value.denominator==1
        w[i][k]=w[k][i]=value
    assert sum(w[i][k] for i,k in combinations(range(16),2))==D*16**2
    i,k=next((i,k) for i,k in combinations(range(16),2) if w[i][k])
    assert determinant(rows+[V[i],V[k]])**2==w[i][k]
    groups=[]
    for line in classes.values():
        r=line[0];axis=next(j for j in range(2) if z[r][j])
        ratios=[(z[u][axis]/z[r][axis])**2 for u in line]
        outside=[u for u in nonzero if u not in line]
        weights=[w[r][u] for u in outside]
        assert all(w[u][k]==t*w[r][k] for u,t in zip(line,ratios) for k in outside)
        C=sum(ratios)*sum(weights);assert C>0 and C.denominator==1
        key=(3+len(nonzero),tuple(sorted((hist(ratios),hist(weights)))))
        groups.append((tuple(line),int(C),key))
    assert sum(C for _,C,_ in groups)==2*D*16**2
    return D,nonzero,groups
def ar(x):
    x=F(x);return arb(x.numerator)/x.denominator
def parse(s):
    a,b=map(F,s[1:-1].split(' +/- '));return a-b,a+b
def iroot(n,k):
    x=1<<((n.bit_length()+k-1)//k)
    while True:
        y=((k-1)*x+n//x**(k-1))//k
        if y>=x:break
        x=y
    while x**k>n:x-=1
    while (x+1)**k<=n:x+=1
    return x

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,required=True);parser.add_argument('--model-root',type=Path,required=True)
    parser.add_argument('--run',action='store_true');args=parser.parse_args()
    source,model=proj(args.source_root),proj(args.model_root);author=source/'upper_I021';manifest=author/'manifest.json'
    assert sha(manifest)==MANIFEST;locked={'upper_I021/manifest.json':MANIFEST}
    info=json.loads(manifest.read_text())
    for row in info['white_list']:
        p=author/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes'];locked['upper_I021/'+row['path']]=sha(p)
    resolved={}
    for row in info['source_inputs']:
        suffix=Path(row['source']).parts[-2:]
        options=[root.joinpath(*suffix) for root in (source,model)]
        p=next((p for p in options if p.exists() and sha(p)==row['sha256']),None);assert p is not None
        resolved[suffix]=p;locked['/'.join(suffix)]=sha(p)
    transfer=model/'alternative_alpha/transfer.md';assert sha(transfer)==TRANSFER;locked['alternative_alpha/transfer.md']=TRANSFER
    result=json.loads((author/'results.json').read_text())
    if not args.run:
        print(json.dumps({'run_handle':RUN,'status':'locked_inputs_only','calculation_executed':False}));return
    began=time.monotonic();triple_info={};representatives=[]
    for J in ((0,1,2),(0,1,6),(0,3,5)):
        D,nz,groups=representative(J);scope=tuple(sorted(set(J)|set(nz)))
        images=orbit((J,scope));assert len({j for j,_ in images})==len(images)
        for j,s in images:
            assert j not in triple_info;triple_info[j]=(s,D,len(nz),len(groups))
        representatives.append((J,D,nz,groups,len(images)))
    assert len(triple_info)==math.comb(16,3)==560
    types={};covered={};occupancy=[0]*16;flag_orbits=[]
    for J,D,nz,groups,count in representatives:
        for line,C,key in groups:
            state=(J,line)
            if state in covered:
                assert covered[state]==(C,key);continue
            images=orbit(state);assert not(set(images)&set(covered))
            for j,line_image in images:
                assert j in triple_info
                covered[j,line_image]=(C,key)
                for i in triple_info[j][0]:occupancy[i]+=C
            row=types.setdefault(key,{'count':0,'sum_C':0});row['count']+=len(images);row['sum_C']+=len(images)*C
            flag_orbits.append({'J':J,'line':line,'orbit_size':len(images),'C':C})
    W=16**5;mass=20*W;assert len(covered)==2240 and sum(r['sum_C'] for r in types.values())==mass
    assert len(set(occupancy))==1;rho=F(occupancy[0],mass);assert rho==F(247,256)
    expected={}
    for row in result['geometry']['probability_types']:
        pair=tuple(sorted(tuple(sorted((F(p),n) for p,n in h.items())) for h in row['probability_histograms']))
        key=(row['scope_size'],pair);expected[key]={'count':row['group_count'],'sum_C':row['sum_C']}
        assert F(row['theta'])==F(row['sum_C'],mass)
    assert types==expected and len(types)==5
    triple_types=Counter((D,nz,groups) for _,D,nz,groups in triple_info.values())
    assert sorted((k,v) for k,v in triple_types.items())==[((48,12,3),160),((64,12,4),240),((112,13,5),160)]
    assert sum(D for _,D,_,_ in triple_info.values())==10*16**3
    ctx.prec=512;beta=F(3,4);b=F(3,20);s=b*rho;q=F(10,7)
    sqrt2,sqrt3=arb(2).sqrt(),arb(3).sqrt()
    R=[arb(2)**ar(F(1,4))+arb(2)**ar(F(-1,4)),arb(2)**ar(F(5,4)),
       (7+8*sqrt2+3*sqrt3)/arb(2)**ar(F(15,4)),(1+1/sqrt2)*arb(3)**ar(F(1,4)),1+sqrt2]
    theta=[F(3,16),F(3,8),F(1,8),F(3,16),F(1,8)]
    for row,r in zip(result['geometry']['probability_types'],R):
        raw=arb(1)
        for h in row['probability_histograms']:raw*=sum((n*ar(p)**ar(beta) for p,n in h.items()),arb(0))
        assert (raw-r).contains(0)
    m=ar(s/beta).gamma()/(ar(beta)*ar(s).gamma()) # unshifted Mellin arguments
    gamma7=ar(beta)*ar(beta).gamma()
    I7=ar(F(3**10,2**65*60**7))*gamma7**10*m**ar(50/rho)/arb.pi()**11
    for r,t in zip(R,theta):I7*=r**ar(-2*t)
    lower7,upper7=F(str(I7.lower().fmpq())),F(str(I7.upper().fmpq()))
    scale=10**95;scale7=scale**7
    ilo=F(iroot(lower7.numerator*scale7//lower7.denominator,7),scale)
    ihi=F(iroot(upper7.numerator*scale7//upper7.denominator,7)+1,scale)
    assert ilo**7<=lower7<upper7<=ihi**7
    oldlo,oldhi=parse(result['arb_384bit']['complete_I_upper']);assert oldlo<=ilo<ihi<=oldhi
    cap=F(result['strict_rational_I_upper']);assert ihi<cap
    oldI=json.loads(resolved[('upper_I020','results.json')].read_text())
    old020low=parse(oldI['arb_384bit']['fractional_Holder_complete_I_upper'])[0]
    assert ihi<old020low and ilo>F(11,10**6)
    g=F(2821,50000000);transfer_cap7=F(10,3)**10*cap**7/g**10
    assert transfer_cap7>32**7
    out={'run_handle':RUN,'status':'pass_full_I021_bound_with_independent_orbit_geometry',
         'geometry':{'representative_triples':3,'triple_orbit_sizes':[r[-1] for r in representatives],
                     'triple_count':len(triple_info),'flag_count':len(covered),'flag_orbits':flag_orbits,
                     'probability_types':len(types),'mass_20W':mass,'coordinate_occupancy':occupancy,'rho':str(rho),
                     'five_types_match_locked_input':True,'three_Leibniz_completions_pass':True},
         'Gamma_Mellin':{'beta':str(beta),'b':str(b),'s_b_rho':str(s),'s_over_beta':str(s/beta),
                        'I_seventh_power_ball_512bit':str(I7),'shifted_Gamma_not_used':True},
         'formula_I_enclosure':[str(ilo),str(ihi)],'strict_rational_I_upper':str(cap),
         'formula_improves_locked020':True,'not_a_lower_bound_for_actual_I':True,
         'transfer':{'same_G_general_L1_all_scales':True,'current_g':'2821/50000000','displayed_certificate_improves_32':False},
         'shared_library_limit':'Arb/flint is shared for Gamma and algebraic interval functions. Geometry uses independent nullspace/orbit code; seventh-power formula and Fraction root enclosure differ from production. No cross-library Gamma certification claimed.',
         'locked_inputs_sha256':locked,'source_sha256':sha(Path(__file__)),
         'registration_sha256':sha(HERE/'registration.md'),
         'runtime':{'python':sys.executable,'arb_bits':ctx.prec,'seconds':time.monotonic()-began,'production_560_Gram_loop_replayed':False}}
    (HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({key:out[key] for key in ('run_handle','status','formula_I_enclosure','strict_rational_I_upper','runtime')},indent=2))

if __name__=='__main__':main()
