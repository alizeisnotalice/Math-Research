"""Second n3 review: own rays, recursive Hermite, exact algebra and endpoints.
Only data/Arb/SymPy are shared; no prior checker/evaluator is imported.
Default is finite: all g patches, eight diagnostic I patches, all saved chords.
Use --scalar-only to omit fresh nodes; --all-I replays all 8,256 I patches.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product, permutations, combinations
from math import comb, prod
import argparse, hashlib, json, time, sys
import sympy as sy
from flint import arb, arb_mat, ctx

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
ctx.prec=512
def aq(q):
    q=Q(q); return arb(q.numerator)/q.denominator
def fq(a): return Q(str(a.fmpq()))
def endpoints(a): return [str(fq(a.lower())),str(fq(a.upper()))]
def unpack(x):
    assert len(x)==4 and int(x[2])==0
    return Q(int(x[0]))*Q(2)**int(x[1])
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def det(v): return int(sy.Matrix.hstack(*map(sy.Matrix,v)).det())
def key(v): return tuple(sorted(v))
def dot(x,y): return sum((a*b for a,b in zip(x,y)),Q(0))

def geometry():
    body=[(1,a,b) for a,b in product((-1,1),repeat=2)]
    face=[]
    for i,j in combinations(range(3),2):
        for s in (-1,1):
            v=[0]*3;v[i]=1;v[j]=s;face.append(tuple(v))
    directions=body+face;weights=[Q(1)]*4+[Q(1,3)]*6
    psi=lambda v:sum((w*abs(dot(v,d)) for d,w in zip(directions,weights)),Q(0))
    bases=[((1,0,0),(1,1,0),(2,1,1)),((1,1,0),(2,1,1),(1,1,1))]
    cones={}
    for p in permutations(range(3)):
        for signs in product((-1,1),repeat=3):
            for base in bases:
                rr=tuple(tuple(signs[i]*r[p[i]] for i in range(3)) for r in base)
                assert abs(det(rr))==1
                for d in directions:
                    v=[dot(r,d) for r in rr]; assert min(v)>=0 or max(v)<=0
                kk=tuple(psi(r) for r in rr)
                assert all(k>=3*(abs(r[1])+abs(r[2])) for r,k in zip(rr,kk))
                assert key(rr) not in cones
                cones[key(rr)]=(rr,kk,Q(1)/prod(kk))
    assert len(cones)==96
    pub=read(ROOT/'upper/continuous_bridge/cb11_cones.json')['cones']
    assert set(cones)=={key(tuple(map(tuple,c['rays']))) for c in pub}
    for c in pub:
        rr=tuple(map(tuple,c['rays']));assert det(rr)==c['det']
        for r,(a,b) in zip(rr,c['rates_body_face_coefficients']):
            assert a==sum(abs(dot(r,v)) for v in body)
            assert b==sum(abs(dot(r,v)) for v in face)
    assert sum((w for rr,kk,w in cones.values()),Q(0))/8==Q(15,544)
    active=[(v,w) for v,w in zip(directions,weights) if v[0]]
    ds=[max(abs(det((v,*bb))) for bb in combinations(body,2)) for v,w in active]
    assert ds==[4]*8
    assert len(active)**2*sum((w/d for (v,w),d in zip(active,ds)),Q(0))==Q(256,3)
    assert sum((w for v,w in active),Q(0))==Q(16,3)
    # Exact proof identities for cones, integrated DD, repeated and zero jets.
    a,b,c,e=sy.symbols('a b c e',real=True)
    assert sy.Matrix.hstack(*map(sy.Matrix,bases[0])).inv()*sy.Matrix([a,b,c])==sy.Matrix([a-b-c,b-c,c])
    assert sy.Matrix.hstack(*map(sy.Matrix,bases[1])).inv()*sy.Matrix([a,b,c])==sy.Matrix([b-c,a-b,b+c-a])
    bs=(a,b,c)
    derivative=sum(e/(e*e+x*x)/prod(x-y for y in bs if y!=x) for x in bs)
    heat=e*(e*e-a*b-a*c-b*c)/prod(e*e+x*x for x in bs)
    assert sy.factor(derivative+heat)==0
    return list(cones.values()),{'cones':96,'rays':len({r for rr,_,_ in cones.values() for r in rr}),'peak_pi3': '15/544','active_basis_dets':ds,'Cf_pi3':'256/3','A':'16/3','DD_derivative_identity':True}

def evaluator(cones,eps):
    """Full96 terms; sorted Newton recursion, confluent jets at exact equality."""
    eps2=eps*eps
    rayset={(r,k) for rr,kk,w in cones for r,k in zip(rr,kk)}
    def point(y,z):
        proj={(r,k):(r[0]+r[1]*y+r[2]*z)/k for r,k in rayset}
        values={};dds={}
        def jet(b,n):
            if n==0:
                if abs(b) not in values: values[abs(b)]=aq(eps2+b*b).log()/2
                return values[abs(b)]
            if n==1:return aq(b/(eps2+b*b))
            return aq((eps2-b*b)/(2*(eps2+b*b)**2))
        def dd(t):
            if t not in dds:
                if t[0]==t[-1]: v=jet(t[0],len(t)-1)
                else:v=(dd(t[1:])-dd(t[:-1]))*aq(1/(t[-1]-t[0]))
                dds[t]=v
            return dds[t]
        total=arb(0)
        for rr,kk,w in cones:
            t=tuple(sorted(proj[r,k] for r,k in zip(rr,kk)))
            total+=aq(w)*dd(t)
        return total/(8*arb.pi()**3)
    return point

def transforms(d,grid,eps,rho):
    N=d+1;h=Q(1,2*grid);delta=eps-h*(rho-1/rho)/6
    assert delta>0
    M=aq(Q(15,1088)/delta**2)/arb.pi()**3
    error=(2*d+2)*4*M*aq(rho)**(-d)/aq(rho-1)
    D=aq(Q(15,1088)/((eps-eps/3)**2*eps))/arb.pi()**3
    displacement=2*D*aq(h/Q(2**128))
    roots=[Q(0)]*N;eta=Q(0)
    for j in range(d//2):
        x=(arb.pi()*Q(2*j+1).numerator/(2*N)).cos()
        f=(x*2**128).floor().unique_fmpz();assert f is not None
        roots[j]=Q(int(f),2**128);roots[d-j]=-roots[j]
        assert abs(x-aq(roots[j]))<arb(2)**-128
        eta=max(eta,fq(abs(x-aq(roots[j])).upper()))
    # A sharper valid displacement prevents an unjustified expectation that
    # differently rounded interval-arithmetic paths must automatically nest.
    assert 0<eta<Q(1,2**128)
    displacement=2*D*aq(h*eta)
    u=sy.Symbol('u');P=[sy.Poly(sy.chebyshevt(k,2*u-1),u) for k in range(N)]
    B=[[sum((Q(P[k].nth(l))*Q(comb(j,l),comb(d,l)) for l in range(min(j,k)+1)),Q(0)) for k in range(N)] for j in range(N)]
    mu=[sum((Q(P[k].nth(l))/Q(l+1) for l in range(k+1)),Q(0)) for k in range(N)]
    for k in range(N):
        back=sy.Poly(sum(sy.Rational(B[j][k].numerator,B[j][k].denominator)*comb(d,j)*u**j*(1-u)**(d-j) for j in range(N)),u)
        assert back==P[k]
        assert mu[k]==(Q(1) if k==0 else Q(0) if k%2 else Q(-1,k*k-1))
        assert sum((B[j][k] for j in range(N)),Q(0))/N==mu[k]
    A=arb_mat([[(arb.pi()*k*(2*j+1)/(2*N)).cos()*aq(Q(1 if k==0 else 2,N)) for j in range(N)] for k in range(N)])
    T=arb_mat([[aq(x) for x in row] for row in B])*A
    Qm=arb_mat([[(arb.pi()*k*(2*j+1)/(2*N)).cos() for k in range(N)] for j in range(N)])
    identity=A*Qm
    assert all(identity[i,j].contains(int(i==j)) and identity[i,j].rad()<arb(2)**-450 for i in range(N) for j in range(N))
    MA=arb_mat([[sum((aq(mu[k])*A[k,j] for k in range(N)),arb(0)) for j in range(N)]])
    return roots,T,MA,error,displacement,{'degree':d,'delta':str(delta),'M':endpoints(M),'error':endpoints(error),'node_displacement':endpoints(displacement),'certified_max_root_displacement':str(eta),'polynomial_identities':N,'DCT_entries':N*N}

def patch(point,setup,i,j,grid):
    roots,T,MA,error,displacement,_=setup;N=len(roots)
    ys=[Q(2*i+1,2*grid)+r/(2*grid) for r in roots]
    zs=[Q(2*j+1,2*grid)+r/(2*grid) for r in roots]
    V=arb_mat([[point(y,z)+arb(0,displacement.upper()) for z in zs] for y in ys])
    assert all(V[a,b].is_finite() for a in range(N) for b in range(N))
    C=T*V*T.transpose();mean=(MA*V*MA.transpose())[0,0]
    L=max(Q(0),fq((min(C[a,b].lower() for a in range(N) for b in range(N))-error).lower()))
    U=fq((max(C[a,b].upper() for a in range(N) for b in range(N))+error).upper())
    m=min(U,fq((mean+error).upper()));assert L<=m<=U
    return L,U,m

def saved_scalars(raw):
    s=ROOT/'upper/continuous_bridge/cb11_integral_cheb_final_staging'
    prov=read(s.parent/'cb11_integral_cheb_final_provenance.json')
    records={};total=Q(0);weights=0;cache=[];fresh=0
    cfg={'schema':1,'cone_sha256':sha(s/'cb11_cones.json'),'epsilon':'9/5000','tau':'1/3','degree':32,'grid':128,'node_bits':128,'precision_bits':384,'rho':'5/2'}
    checker=sha(s/'cb11_integral_cheb_check.py')
    for part in range(12):
        name=f'cb11_integral_cheb_part_{part:02}.json';data=read(s/name)
        assert data['config']==cfg and data['checker_sha256']==checker
        assert data['part']==part and data['parts']==12 and len(data['patches'])==688
        if raw:
            f=raw/name;assert sha(f)==prov['raw_parts'][name]['sha256'];original=read(f)
            assert original['config']==cfg and original['checker_sha256']==checker
        for number,r in enumerate(data['patches']):
            i,j=r['patch'];assert 0<=j<=i<128 and (i*(i+1)//2+j)%12==part and (i,j) not in records
            L,U,m,c=map(unpack,[r['L'],r['U'],r['mean_H_upper'],r['IH_contribution_upper']])
            assert Q(r['IH_contribution_exact'])==c and 0<=L<=m<=U and U>0
            w=1 if i==j else 2;weights+=w
            a,b=aq(L).sqrt(),aq(U).sqrt()
            own=(a**3+(b*b+a*b+a*a)/(a+b)*aq(m-L))*w/128**2
            assert fq(own.upper())<=c
            total+=c;records[i,j]=(L,U,m,c)
            if raw:
                rr=original['patches'][number];assert all(rr[k]==v for k,v in r.items())
                assert rr['weight']==w and rr['nodes']==1089
                assert list(map(Q,rr['whole_square']))==[Q(i,128),Q(i+1,128),Q(j,128),Q(j+1,128)]
                assert unpack(rr['chord_mean_upper'])*w/128**2==c
                if rr['cache_hit']:cache.append([i,j])
                else:fresh+=1089
    assert len(records)==8256 and weights==128**2
    if raw:
        assert sorted(cache)==sorted([[0,0],[127,0],[127,57],[127,127]]) and fresh==8986428
    allr=read(s/'cb11_integral_cheb_all_result.json')
    assert total==Q(allr['IH_upper_exact'])==Q(prov['IH_upper_exact'])
    E=aq(Q(1,16)).root(3)*aq(Q(3,4))*(aq(Q(256,3))/arb.pi()**3).root(3)*(aq(Q(16,3))/arb.pi()).root(3)**2*aq(Q(9,5000)).root(3)**4
    v=aq(total).root(3)**2+E;I=v*v.sqrt()
    assert fq(I.upper())<=unpack(allr['I_upper'])<=Q(23129127,40000000000)
    g=Q(2621323989,500000000000);cap=Q(1979,250);I0=Q(23129127,40000000000)
    cert=read(ROOT/'audit/033_results.json');gp={tuple(r['patch']):Q(r['H_lower_exact']) for r in cert['patches']}
    assert set(gp)=={(i,j) for i in range(16) for j in range(i+1)}
    assert min(gp.values())==Q(cert['minimum_H_lower_exact'])>g
    cost2=27*I0**2/g**3;assert cap*cap>cost2
    result={'chords':len(records),'weighted_area':'1','B_exact':str(total),'short_time_E':endpoints(E),'I':endpoints(I),'g':str(g),'I_outward':str(I0),'cost_squared':str(cost2),'strict_square_gap':str(cap*cap-cost2),'cost':str(aq(cost2).sqrt()),'raw_records':8256 if raw else 0,'original_fresh_nodes':fresh if raw else None,'original_cached_patches':cache}
    return records,gp,result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--scalar-only',action='store_true');ap.add_argument('--all-I',action='store_true');ap.add_argument('--raw-archive',type=Path);args=ap.parse_args()
    start=time.monotonic();snap=read(OUT/'input_snapshot.json')['files']
    for p,h in snap.items():assert sha(ROOT/p)==h,p
    cones,geom=geometry();records,gp,scalars=saved_scalars(args.raw_archive)
    gsetup=transforms(20,16,Q(1,50),Q(3));isetup=transforms(32,128,Q(9,5000),Q(5,2))
    result={'geometry':geom,'scalars':scalars,'g_constants':gsetup[-1],'I_constants':isetup[-1],'fresh_g_nodes':0,'fresh_I_nodes':0,'fresh_patches':[],'shared_libraries':{'Arb':'python-flint512bit','SymPy':sy.__version__},'source_lock_files':len(snap)}
    if not args.scalar_only:
        pg=evaluator(cones,Q(1,50));pi=evaluator(cones,Q(9,5000))
        # Full minimum rerun has only59976 nodes; eliminates dependence on old g endpoints.
        low=None
        for count,(i,j) in enumerate(sorted(gp),1):
            L,U,m=patch(pg,gsetup,i,j,16);assert L>Q(2621323989,500000000000)
            low=L if low is None else min(low,L);result['fresh_g_nodes']+=441
            if count%16==0 or count==136:print(json.dumps({'g_patches':count,'minimum':str(aq(low)),'elapsed':time.monotonic()-start}),flush=True)
        result['fresh_g_minimum_exact']=str(low)
        # Include original four cached cells, fresh diagonal/interior, largest chord and range.
        wanted={(0,0),(127,0),(127,57),(127,127),(71,34),(64,64),max(records,key=lambda ij:records[ij][3]),max(records,key=lambda ij:records[ij][1]-records[ij][0])}
        if args.all_I:wanted=set(records)
        for i,j in sorted(wanted):
            L,U,m=patch(pi,isetup,i,j,128);oldL,oldU,oldm,c=records[i,j]
            # Direct enclosure dominance, not numerical proximity or a sample pass.
            assert L>=oldL and U<=oldU and m<=oldm,(i,j,str(aq(L-oldL)),str(aq(oldU-U)),str(aq(oldm-m)))
            result['fresh_patches'].append({'patch':[i,j],'L':str(L),'U':str(U),'mean_upper':str(m),'dominates_recorded_range_and_mean':True})
            result['fresh_I_nodes']+=1089
            print(json.dumps({'I_patch':[i,j],'dominance':True,'elapsed':time.monotonic()-start}),flush=True)
    for p,h in snap.items():assert sha(ROOT/p)==h,p
    result.update(status='PASS reviewed n3 numeric/geometric certificate; shared maximal theorem audited separately',checker_sha256=sha(Path(__file__)),seconds=time.monotonic()-start)
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','fresh_g_nodes','fresh_I_nodes','seconds']},indent=2))
if __name__=='__main__': main()
