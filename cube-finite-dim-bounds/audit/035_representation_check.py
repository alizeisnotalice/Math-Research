"""Independent4D cone/star/global/g reconstruction, no author producers."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,permutations,product
from collections import Counter
from math import prod,comb
import argparse,hashlib,json,time
import sympy as sy
from flint import fmpq_mpoly_ctx,fmpq,arb,ctx as actx
HERE=Path(__file__).resolve().parent

def run(root):
    start=time.monotonic();sources=json.loads((HERE/'035_source_snapshot.json').read_text())
    for name,h in sources.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source drift')
    V=[(1,)+v for v in product((-1,1),repeat=3)]
    rays0={'A':(1,0,0,0),'B':(1,1,0,0),'C':(1,1,1,0),'D':(1,1,1,1),'E':(2,1,1,0),'F':(3,1,1,1)}
    labels=[('A','B','E','F'),('E','B','C','D'),('E','B','D','F')]
    physical=[tuple(rays0[k] for k in row) for row in labels]
    a,b,c,d=sy.symbols('a b c d',real=True);inverse_expected=[(a-b-c-d,b-c,c-d,d),(a-b,b-c,b+c-a-d,d),(c-d,b-c,(-a+b+c+d)/2,(a-b-c+d)/2)]
    basedets=[]
    for rs,inv in zip(physical,inverse_expected):
        R=sy.Matrix.hstack(*(sy.Matrix(v) for v in rs));basedets.append(abs(int(R.det())))
        assert R.inv()*sy.Matrix([a,b,c,d])==sy.Matrix(inv)
    assert basedets==[1,1,2]
    cones={}
    for p in permutations(range(4)):
        for signs in product((-1,1),repeat=4):
            for rs,det in zip(physical,basedets):
                rr=tuple(tuple(signs[i]*v[p[i]] for i in range(4)) for v in rs)
                for v in V:
                    z=[sum(x*y for x,y in zip(v,r)) for r in rr]
                    assert min(z)>=0 or max(z)<=0
                key=tuple(sorted(rr));assert key not in cones;cones[key]=(rr,det)
    assert len(cones)==1152
    rays=sorted({r for rr,_ in cones.values() for r in rr});assert len(rays)==240
    rates={r:Q(sum(abs(sum(x*y for x,y in zip(v,r))) for v in V)) for r in rays}
    assert [rates[rays0[k]] for k in ('A','B','C','D','E','F')]==[8,8,12,12,16,24]
    for axis in range(4):
        if axis<3:trans=lambda r,j=axis:tuple(r[j+1] if k==j else r[j] if k==j+1 else r[k] for k in range(4))
        else:trans=lambda r:tuple(-v if k==0 else v for k,v in enumerate(r))
        image={tuple(sorted(trans(r) for r in rr)):det for rr,det in cones.values()}
        assert image=={key:det for key,(_,det) in cones.items()}
        assert all(rates[trans(r)]==rates[r] for r in rays)
    volume=384*sum((Q(det)/(24*prod(rates[r] for r in rr)) for rr,det in zip(physical,basedets)),Q(0))
    peak=24*volume/Q(16);assert volume==Q(11,4608) and peak==Q(11,3072)
    ranks={k:sy.Matrix([v for v in V if sum(x*y for x,y in zip(v,r))==0]).rank() if any(sum(x*y for x,y in zip(v,r))==0 for v in V) else 0 for k,r in rays0.items()}
    assert [ranks[k] for k in ('A','B','C','D','E','F')]==[0,3,0,3,2,1]
    assert sum(ranks[next(k for k,v in rays0.items() if sorted(map(abs,v))==sorted(map(abs,r)))]==3 for r in rays)==40
    t=sy.symbols('t',positive=True)
    assert sy.integrate(1/((t*t+1)*(t*t+4)),(t,0,sy.oo))==sy.pi/12
    assert sy.integrate(1/(t*t+1)**2,(t,0,sy.oo))==sy.pi/4
    for nodes,wanted in [((-1,-2,1,2),Q(-1,6)),((-1,-1,1,1),Q(-1,2))]:
        zz=sorted(map(Q,nodes));values=[Q(1 if z>0 else -1) for z in zz]
        for order in range(1,4):values=[Q(0) if zz[i+order]==zz[i] else (values[i+1]-values[i])/(zz[i+order]-zz[i]) for i in range(4-order)]
        assert values[0]==wanted
    CTX=fmpq_mpoly_ctx.get(('a','b','c','d'),ordering='lex');xs=CTX.gens();zero=CTX.constant(0);one=CTX.constant(1)
    def const(q):q=Q(q);return CTX.constant(fmpq(q.numerator,q.denominator))
    def lin(v):return sum((const(v[i])*xs[i] for i in range(4)),zero)
    def poly(dic):return CTX.from_dict({tuple(map(int,e.split(','))):fmpq(v) for e,v in dic.items()})
    def add(n,d,p,q):
        g=d.gcd(q);dl=d//g;ql=q//g;nn=n*ql+p*dl;dd=dl*q
        if not nn:return zero,one
        common=nn.gcd(dd);return nn//common,dd//common
    stored_stars=json.loads((root/'upper/u11_minimum_star_results.json').read_text())['stars']
    computed={};star_records=[]
    for name,r in rays0.items():
        nn,dd=zero,one;count=0
        for rr,det in cones.values():
            if r not in rr:continue
            denominator=one
            for v in rr:
                if v!=r:denominator*=lin(tuple(Q(r[k])/rates[r]-Q(v[k])/rates[v] for k in range(4)))
            nn,dd=add(nn,dd,const(Q(det)/prod(rates[v] for v in rr)),denominator);count+=1
        computed[name]=(nn,dd)
        if name in ('B','D'):
            row=stored_stars[0 if name=='B' else 1];assert nn*poly(row['denominator_dict'])==dd*poly(row['numerator_dict'])
        else:assert nn==0
        star_records.append({'representative':name,'terms':count,'zero':not bool(nn),'exact_identity':True})
        print('035 star',name,'PASS',flush=True)
    # Independently substitute the two representative stars into all20 lines.
    terms=[]
    for i,j in combinations(range(4),2):
        for sign in (-1,1):
            rest=[k for k in range(4) if k not in (i,j)];sub=[xs[i],sign*xs[j],xs[rest[0]],xs[rest[1]]]
            ray=[0]*4;ray[i]=1;ray[j]=sign;n,dn=computed['B']
            terms.append((tuple(ray),n.compose(*sub),dn.compose(*sub)))
    for signs in product((-1,1),repeat=3):
        ray=(1,)+signs;sub=[ray[k]*xs[k] for k in range(4)];n,dn=computed['D']
        terms.append((ray,n.compose(*sub),dn.compose(*sub)))
    raw=json.loads((root/'upper/u11_minimum_global_results.json').read_text())['physical_cases']
    cert=json.loads((root/'upper/u11_minimum_certificate_results.json').read_text())['physical_cases']
    gtarget=Q(21,25000);piup=Q(22,7);ktarget=gtarget*piup**3;assert ktarget==Q(3993,153125)
    actx.prec=256;assert arb.pi()<arb(22)/7
    case_records=[]
    for idx,(row,witness,rs) in enumerate(zip(raw,cert,physical)):
        n,dn=zero,one;signs=[]
        for ray,star_n,star_d in terms:
            co=[sum(ray[k]*v[k] for k in range(4)) for v in rs]
            assert any(co) and (min(co)>=0 or max(co)<=0)
            sign=1 if max(co)>0 else -1;signs.append(sign)
            n,dn=add(n,dn,const(Q(-sign,16))*star_n,star_d)
        N,D=poly(row['numerator_dict']),poly(row['denominator_dict']);assert n*D==dn*N
        sub=[sum((rs[j][i]*xs[j] for j in range(4)),zero) for i in range(4)]
        dc,df=D.factor();factor_count=0
        for f,m in df:
            assert f.total_degree()==1
            coefficients=list(f.compose(*sub).to_dict().values());assert coefficients
            sg=-1 if all(v<=0 for v in coefficients) else 1
            assert all(sg*v>=0 for v in coefficients) and any(sg*v>0 for v in coefficients)
            dc*=sg**m;factor_count+=int(m)
        assert dc>0 and factor_count==20
        P=(xs[0]**3*N-const(ktarget)*D).compose(*sub);coeff=P.to_dict()
        assert coeff and all(v>=0 for v in coeff.values())
        assert coeff=={tuple(e['powers']):fmpq(e['value']) for e in witness['P_coefficients']}
        assert len(coeff)==[1739,1702,1702][idx]
        case_records.append({'case':row['case'],'global_rational_identity':True,'all20_signs_on_entire_open_cone':signs,'D_positive_factor_multiplicity':factor_count,'P_nonnegative_terms':len(coeff)})
        print('035 global/g case',row['case'],'PASS',flush=True)
    # Independent full-cone confluent pole accumulation at prelocked degenerate points.
    degenerate=[]
    for point,expected in [((1,0,0,0),Q(3,32)),((1,1,0,0),Q(1,24)),((1,1,1,0),Q(17,640))]:
        residues={}
        for rr,det in cones.values():
            nodes=[sum(Q(a)*b for a,b in zip(point,r))/rates[r] for r in rr];groups=Counter(nodes);alpha=Q(det)/prod(rates[r] for r in rr)
            for b,m in groups.items():
                p=[Q(1)]+[Q(0)]*(m-1)
                for a,ma in groups.items():
                    if a==b:continue
                    series=[Q((-1)**l*comb(ma+l-1,l))/(b-a)**(ma+l) for l in range(m)]
                    p=[sum((p[h]*series[l-h] for h in range(l+1)),Q(0)) for l in range(m)]
                for k in range(1,m+1):residues[b,k]=residues.get((b,k),Q(0))+alpha*p[m-k]
        assert all(value==0 for (b,k),value in residues.items() if b==0)
        C=-sum((value*(1 if b>0 else -1) for (b,k),value in residues.items() if b and k==1),Q(0))/32
        assert C==expected
        degenerate.append({'x':point,'zero_poles_cancel_exactly':True,'G_times_pi3':str(C)})
    pairdets=[]
    for i,j in combinations(range(8),2):
        ds=[abs(int(sy.Matrix.hstack(*(sy.Matrix(V[k]) for k in (i,j,a,b))).det())) for a,b in combinations([k for k in range(8) if k not in (i,j)],2)]
        assert any(ds);pairdets.append(min(d for d in ds if d))
    assert min(pairdets)==8
    return {'status':'PASS independently reconstructed1152cone/wholeGreen/star/global/g, fixed g>=21/25000',
            'cone_count':len(cones),'ray_count':len(rays),'base_determinants':basedets,'spectrum_volume':str(volume),'peak_times_pi4_t4':str(peak),'true_vertex_count':40,'ray_zero_ranks':ranks,'star_checks':star_records,'case_checks':case_records,'degenerate_fullsum_checks':degenerate,'all28_pair_completions_min_det':min(pairdets),'fixed_target':'21/25000','pi_upper':'22/7','k':str(ktarget),'source_sha256':sources,'seconds':time.monotonic()-start,'shared_libraries':'SymPy/Fraction/flint fmpq_mpoly/Arb256bit;no author producer/checker/module executed;inputstar/global dictionaries only exact cross-multiplication targets'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');a=ap.parse_args();out=run(a.project);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if a.save:(HERE/'035_representation_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='source_sha256'},indent=2))
