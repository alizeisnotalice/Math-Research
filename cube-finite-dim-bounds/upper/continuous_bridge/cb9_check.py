"""Preregistered 12 generic + 6 degenerate spatial points, no cost search."""
from pathlib import Path
from fractions import Fraction as F
import runpy,json,math,time,hashlib
import mpmath as mp
from scipy.integrate import quad
OUT=Path(__file__).resolve().parent
NS=runpy.run_path(str(OUT/'cb9_frame_green.py'));CONES=NS['CONES']
GENERIC=[(1,F(2,7),F(1,11)),(1,F(1,3),F(1,7)),(1,F(2,3),F(1,13)),(1,F(3,5),F(2,9)),(1,F(4,5),F(1,11)),(2,F(1,2),-F(1,7)),(1,-F(2,5),F(1,11)),(F(1,3),F(5,7),1),(F(3,4),-F(2,7),F(4,5)),(F(5,4),F(3,8),-F(2,11)),(F(1,5),-F(4,7),-F(3,2)),(F(7,6),F(5,11),F(2,13))]
DEGENERATE=[('origin',(0,0,0)),('body_diagonal',(1,1,1)),('axis',(1,0,0)),('coordinate_plane_edge',(1,1,0)),('zero_linear_form',(1,1,F(1,3))),('repeated_nonzero_forms',(1,F(2,3),F(1,3)))]

def m(v):return mp.mpf(v.numerator)/v.denominator if isinstance(v,F) else mp.mpf(v)
def green_dd(x):
    # Exact rational classification precedes arbitrary-precision arithmetic.
    total=mp.mpf(0)
    for cone in CONES:
        ell=[sum(F(a)*b for a,b in zip(r,x)) for r in cone]
        if any(a==0 for a in ell):raise ValueError('termwise integral invalid: zero form')
        if len(set(ell))==1:total-=1/(2*m(ell[0])**2)
        elif len(set(ell))==2:
            repeated=next(a for a in ell if ell.count(a)==2);other=next(a for a in ell if a!=repeated);a,b=m(repeated),m(other)
            total+=(mp.log(abs(b))-mp.log(abs(a))-(b-a)/a)/(b-a)**2
        else:
            total+=sum(mp.log(abs(m(ell[j])))/mp.fprod(m(ell[j]-ell[k]) for k in range(3) if k!=j) for j in range(3))
    return total/(16*mp.pi**3)

def positive_density(t,x):
    x=list(map(float,x));b=[0.,-(x[1]+x[2])/2,(x[0]-x[2])/2,(x[0]-x[1])/2]
    breaks=[-math.inf]+sorted(set(b))+[math.inf]
    def f(s):return 1/math.prod((s-v)**2+t*t for v in b)
    ans=sum(quad(f,a,b,epsabs=2e-11,epsrel=2e-11,limit=150)[0] for a,b in zip(breaks,breaks[1:]))
    return t**4*ans/(4*math.pi**4)

def run():
    start=time.monotonic();records=[];t=3/7
    with mp.workdps(60):
        for label,x in [('generic_'+str(i+1),x) for i,x in enumerate(GENERIC)]+DEGENERATE:
            exact=tuple(map(F,x));generic=label.startswith('generic')
            if generic:
                assert all(all(a!=0 for a in ll) and len(set(ll))==3 for c in CONES for ll in [[sum(F(a)*b for a,b in zip(r,exact)) for r in c]]),label
            pt=NS['heat_kernel'](t,x);pc=NS['heat_kernel'](t,x,method='cones');pp=positive_density(t,x)
            relp=max(abs(pt-pc),abs(pt-pp))/pt;assert relp<1e-10,(label,pt,pc,pp)
            value,error=NS['green'](x,return_error=True);high=NS['green'](x,method='mp',dps=60)
            row={'label':label,'x_exact':list(map(str,exact)),'t':'3/7','p_rational':pt,'p_all_complex_cones':pc,'p_positive_Cauchy_convolution':pp,'p_max_relative_difference':relp,'G_float':value if math.isfinite(value) else 'infinity','G_quad_error_estimate':error if math.isfinite(error) else 'infinity','G_mp60':str(high)}
            if math.isfinite(value):
                relg=float(abs(mp.mpf(value)-high)/high);assert relg<1e-9,(label,value,high);row['G_float_vs_mp_relative_difference']=relg
                try:
                    dd=green_dd(exact);reldd=abs(dd-high)/high;assert reldd<mp.mpf('1e-45'),(label,dd,high)
                    row.update(G_complete_cone_DD_mp60=str(dd),G_DD_relative_difference=str(reldd))
                except ValueError:row['G_cone_termwise_integral']='not valid because a zero form occurs; used total positive rational integral'
            else:assert mp.isinf(high) and label in ('origin','body_diagonal')
            if label=='axis':
                anchor=3/(8*mp.pi**3);assert abs(high-anchor)<mp.mpf('1e-55');row['exact_anchor']='3/(8*pi^3)'
            if label=='coordinate_plane_edge':
                anchor=(4*mp.log(2)-1)/(8*mp.pi**3);assert abs(high-anchor)<mp.mpf('1e-55');row['exact_anchor']='(4*log(2)-1)/(8*pi^3)'
            if label=='origin':
                assert abs(pt-5/(64*math.pi**3*t**3))<1e-15;row['exact_heat_anchor']='5/(64*pi^3*t^3)'
            records.append(row)
    result={'status':'PASS preregistered 12 generic + 6 degenerate points','spatial_points':len(records),'mp_dps':60,'seconds':time.monotonic()-start,'module_sha256':hashlib.sha256((OUT/'cb9_frame_green.py').read_bytes()).hexdigest(),'records':records,'scope':'numerical cross-check and exact anchor consistency; no certified g/J or weak-type improvement'}
    (OUT/'cb9_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2));print('max p rel',max(r['p_max_relative_difference'] for r in records));print('max G rel',max(r.get('G_float_vs_mp_relative_difference',0) for r in records))
if __name__=='__main__':run()
