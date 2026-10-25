"""Independent fixed Green quadrature with altered cuts/series; no author import."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,argparse,time
import flint
from flint import arb,acb,ctx
HERE=Path(__file__).resolve().parent
CASES=[(3,F(1,2)),(4,F(1,2)),(3,F(47716536955561023,10**17)),(4,F(23547658784814787,5*10**16))]
CUTS=[F(1,2000),F(1,100),F(1,8),F(1),F(4),F(36)]

def realq(x):return arb(x.numerator)/x.denominator
def complexq(x):return acb(x.numerator)/x.denominator

def tails(n):
    delta=F(1,2000);L=36;power=F(n)*F(47,100)**2/(4*delta);floor=power.numerator//power.denominator
    assert 45**2>2000
    lowg=F(45 if n==3 else 2000,2**floor)
    lowm=2*n*delta/F(2**500)
    highg=F(2**n,2**(81*n));highm=F(4**n,2**(81*n))
    # M9 missing min distance38; ten erfc terms missing odd21; split1/8.
    imageg=F(n*2**(n-1)*8,2**2888);imagem=F(n*2**(n-1)*4,2**882)
    # N44 missing odd89, pi>3 at t>=1/8 gives exponent>2227.
    assert F(9*89**2,32)>2227
    spectral=F(n*2**n,2**2227)
    totalg=lowg+highg+imageg+spectral;totalm=lowm+highm+imagem+spectral
    assert totalg<F(1,10**70) and totalm<F(1,10**70) and max(totalg,totalm)<F(1,10**40)
    return {'n':n,'Green_total_absolute_omission_upper':str(totalg),'mass_total_absolute_omission_upper':str(totalm),'both_less_than':'1e-70','added_radius':'1e-40'}
def interval(z):return [str(z.lower().fmpq()),str(z.upper().fmpq())]

def run(root):
    snapshot=json.loads((HERE/'023_source_snapshot.json').read_text())
    for name,h in snapshot.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source changed')
    old=json.loads((root/'results/round008_public_green.json').read_text())['rows']
    ctx.prec=128;pi=acb(arb.pi());one=acb(1);delta=realq(CUTS[0]);parts_records=[];cacheM={}
    def survival(t,analytic,kind):
        if kind=='images':
            z=t.sqrt(analytic=analytic);return one-2*sum((((-1)**j)*((2*j+1)/(2*z)).erfc() for j in range(10)),acb(0))
        return 4/pi*sum((acb((-1)**j)/(2*j+1)*(-pi*pi*(2*j+1)**2*t/4).exp() for j in range(44)),acb(0))
    def density(t,analytic,kind,space):
        if kind=='images':
            out=sum((((-(space-4*j)**2/(4*t)).exp()-(-(space-4*j-2)**2/(4*t)).exp()) for j in range(-9,10)),acb(0))
            return out/(4*pi*t).sqrt(analytic=analytic)
        return sum((((pi*(2*j+1)*space/2).cos()*(-pi*pi*(2*j+1)**2*t/4).exp()) for j in range(44)),acb(0))
    def integrate(n,space,mode):
        pieces=[]
        for left,right in zip(CUTS,CUTS[1:]):
            kind='images' if right<=F(1,8) else 'spectral'
            if mode=='mass':
                def fn(t,analytic):return survival(t,analytic,kind)**n
            else:
                def fn(t,analytic):return density(t,analytic,kind,space)**n
            z=acb.integral(fn,complexq(left),complexq(right),rel_tol=arb('1e-21'),abs_tol=arb('1e-23'),eval_limit=100000,depth_limit=40,deg_limit=96)
            assert z.is_finite(),(n,mode,left,right,z);pieces.append(z)
        parts_records.append({'n':n,'mode':mode,'a':str(space) if space is not None else None,'parts':[z.str(40) for z in pieces]})
        return sum(pieces,acb(0)).real
    rows=[]
    for i,(n,a) in enumerate(CASES):
        start=time.monotonic();assert F(47,100)<=a<=F(1,2)
        if n not in cacheM:cacheM[n]=integrate(n,None,'mass')+delta+arb(0,arb('1e-40'))
        M=cacheM[n];space=complexq(a);g=integrate(n,space,'green')+arb(0,arb('1e-40'));assert g>0
        cost=M/(2**n*realq(a)**n*g);lo,hi=map(F,interval(cost));oldlo,oldhi=map(F,old[i]['model_cost_interval']);assert oldlo<lo<=hi<oldhi
        assert str(a)==old[i]['a'] and hi<F(old[i]['exact_simple_model_upper'])
        radial=27*arb.pi()*arb(3).sqrt()/16 if n==3 else 2*arb.pi()**2
        assert cost<radial and lo>2**n
        rows.append({'n':n,'a_exact':str(a),'R_exact':str(1/a),'mass_interval':interval(M),'corner_green_interval':interval(g),'model_cost_interval':interval(cost),'model_cost_ball':cost.str(40),'inside_author_cost_interval':True,'author_simple_model_upper':old[i]['exact_simple_model_upper'],'strictly_below_radial_model':True,'strictly_above_current_global_target':True,'seconds':time.monotonic()-start})
    # Exact interval ordering proves model descent; no optimization here.
    assert F(rows[2]['model_cost_interval'][1])<F(rows[0]['model_cost_interval'][0])
    assert F(rows[3]['model_cost_interval'][1])<F(rows[1]['model_cost_interval'][0])
    # Acb endpoints are binary outward enclosures; fmpq conversion kept exact.
    return {'status':'PASS independent fixed-parameter rigorous Green model certification','source_sha256':snapshot,'parameters':{'precision_bits':128,'cuts':list(map(str,CUTS)),'image_M':9,'survival_terms':10,'spectral_terms':44,'components_integrated_separately':True},'independent_tail_bounds':[tails(3),tails(4)],'rows':rows,'certified_finite_integral_balls':parts_records,'libraries':{'python_flint':flint.__version__,'Fraction':'standard Python'},'shared':'same heat-kernel identities and Arb/Acb software kernel; no author code imports; altered partition/cutoffs/precision/component integration','scope':'four legal model upper bounds and strict model descent only; global C3<=8,C4<=16 unchanged'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');args=ap.parse_args();out=run(args.project);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('023_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='certified_finite_integral_balls'},indent=2))
