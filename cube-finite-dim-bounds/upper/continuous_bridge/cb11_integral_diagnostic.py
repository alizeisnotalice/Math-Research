"""11H: finite rigorous interval upper diagnostics, fixed tau=1/3.

247 spatial boxes x (15 time boxes + exact infinite-tail item)=3952 terms.
Initial pointwise short-time budget; one registered global-Lp revision.
No point quadrature, minimum search or tau optimization.
"""
from pathlib import Path
from fractions import Fraction as F
import json,math,time,hashlib
from flint import arb,acb,ctx
OUT=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds')
TAU=F(1,3);TARGET=F(57,100000);RHO=F(1,10**8);DELTA=F(1,2**28)

def af(v):
    v=F(v);return arb(v.numerator)/v.denominator

def cover():
    corners=((F(1),F(0)),(F(0),F(1)),(F(1),F(1)))
    boxes=[];removed=[]
    for i in range(4):
      for j in range(4):
        box=(F(i,4),F(i+1,4),F(j,4),F(j+1,4))
        q=next((q for q in corners if q[0] in box[:2] and q[1] in box[2:] and (q[0]==1 or box[0]==0) and (q[1]==1 or box[2]==0)),None)
        # Only the three corner cells touching their own singular corner.
        if q is None:boxes.append(box);continue
        cur=box
        while cur[1]-cur[0]>DELTA:
            a,b,c,d=cur;m=(a+b)/2;n=(c+d)/2
            kids=[(a,m,c,n),(a,m,n,d),(m,b,c,n),(m,b,n,d)]
            cur=next(k for k in kids if q[0] in k[:2] and q[1] in k[2:])
            boxes.extend(k for k in kids if k!=cur)
        removed.append(cur)
    assert len(boxes)==247 and len(removed)==3
    assert sum((area(b) for b in boxes+removed),F(0))==1
    assert all((b[1]-b[0])==DELTA and (b[3]-b[2])==DELTA for b in removed)
    return boxes,removed,corners

def area(b):return (b[1]-b[0])*(b[3]-b[2])
def r2_lower(b,corners):
    def d(q,lo,hi):return lo-q if q<lo else q-hi if q>hi else F(0)
    return min(d(q[0],b[0],b[1])**2+d(q[1],b[2],b[3])**2 for q in corners)

def prepared():
    doc=json.loads((MAIN/'upper/continuous_bridge/cb11_cones.json').read_text());seen=set();out=[]
    for row in doc['cones']:
        rays=tuple(map(tuple,row['rays']));key=tuple(sorted(rays));opp=tuple(sorted(tuple(-v for v in r) for r in rays));pair=min(key,opp)
        if pair in seen:continue
        seen.add(pair);rates=[F(a)+TAU*b for a,b in row['rates_body_face_coefficients']];W=F(2)/math.prod(rates)
        out.append((rays,rates,af(W)))
    assert len(out)==48
    return out

def time_boxes(revision):
    if revision:pts=[F(2)**k for k in range(-26,5,2)]
    else:
        pts=[F(2)**k for k in range(-6,5)];out=[]
        for i,(a,b) in enumerate(zip(pts,pts[1:])):
            if i<5:out.extend([(a,(a+b)/2),((a+b)/2,b)])
            else:out.append((a,b))
        assert len(out)==15;return out
    out=list(zip(pts,pts[1:]));assert len(out)==15;return out

def H_interval(box,times,cones,pi,tail):
    y=arb(af((box[0]+box[1])/2),af((box[1]-box[0])/2));z=arb(af((box[2]+box[3])/2),af((box[3]-box[2])/2));nodes=[]
    for rays,rates,W in cones:
        b=[(r[0]+r[1]*y+r[2]*z)/af(k) for r,k in zip(rays,rates)];nodes.append((b,W))
    total=arb(0)
    for a,b in times:
        t=arb(af((a+b)/2),af((b-a)/2));p=acb(0)
        for knots,W in nodes:
            term=acb(W)
            for knot in knots:term/=acb(t,-knot)
            p+=term
        value=p.real/(8*pi**3);assert value.is_finite()
        total+=af(b-a)*value
    # The last registered time item encloses the entire infinite tail.
    total+=tail/2+arb(0,tail.abs_upper()/2)
    assert total.is_finite() and total.upper()>0
    return total

def power23(v):return (v.log()*af(F(2,3))).exp()
def power32(v):return v*v.sqrt()
def integer_upper(v):
    w=v.upper();f=float(w);u=F(math.ceil(f*(1+1e-6))+2) if f>1 else F(math.ceil(f*10**12)+2,10**12)
    assert af(u)>v;return u

def run_case(revision,boxes,corners,cap,cones,pi,K):
    start=time.monotonic();times=time_boxes(revision);epsilon=times[0][0];T=times[-1][1];peak=F(15,544);tail=af(peak)/(2*pi**3*af(T)**2)
    A=af(55200)/pi**3;upper=arb(0);h_upper=arb(0);rows=[]
    for idx,box in enumerate(boxes):
        r2=r2_lower(box,corners);assert r2>0
        H=H_interval(box,times,cones,pi,tail);hu=H.upper();ar=af(area(box));short=A*(1+af(400)*af(epsilon)**2/af(r2)).log()
        gu=(hu+short.upper()).upper();contrib=(ar*power32(gu)).upper();hc=(ar*power32(hu)).upper()
        upper+=contrib;h_upper+=hc
        rows.append({'box':list(map(str,box)),'area':str(area(box)),'r_Q_squared_lower':str(r2),'H_interval':str(H),'H_upper_display':float(hu),'short_time_pointwise_upper_ball':str(short),'initial_G_power_contribution_upper_display':float(contrib),'H_power_contribution_upper_display':float(hc),'requires_tighter_continuous_enclosure':True})
    initialI=upper+cap;BH=h_upper+cap
    short_power=K*af(epsilon)**2
    revisedI=power32((power23(BH)+power23(short_power)).upper())
    actual_upper=revisedI if revision else initialI
    assert actual_upper.is_finite()
    passed=actual_upper<af(TARGET);U=integer_upper(actual_upper)
    return {'representation':'global Lp short-time/Minkowski revision' if revision else 'pointwise short-time initial','epsilon_exact':str(epsilon),'time_T_exact':str(T),'finite_time_intervals':[[str(a),str(b)] for a,b in times],'space_boxes':len(boxes),'time_items_per_box':len(times)+1,'space_time_items':len(boxes)*(len(times)+1),'infinite_time_tail_G_upper_ball':str(tail),'cap_I_upper_ball':str(cap),'H_power_integral_upper_ball':str(BH),'uniform_pointwise_I_upper_ball':str(initialI),'global_short_time_Lp_power_upper_ball':str(short_power),'Minkowski_I_upper_ball':str(revisedI),'used_I_upper_expression_ball':str(actual_upper),'checked_rational_I_upper':str(U),'target_certified':bool(passed),'remaining_boxes':[] if passed else rows,'seconds':time.monotonic()-start,'scope':'upper-expression enclosure; its lower endpoint is not an actual integral lower bound'}

if __name__=='__main__':
    ctx.prec=128;pi=arb.pi();boxes,removed,corners=cover();assert af(2).sqrt()*af(DELTA)<af(RHO)
    data=json.loads((MAIN/'audit/030_results.json').read_text());case=next(c for c in data['cases'] if c['n']==3 and c['tau']=='1/3');bud=next(b for b in case['budgets'] if b['epsilon_exact']=='1/100000000')
    capq=F(bud['all_face_caps_total_interval'][1])/24;cap=af(capq)
    # q=3/2 radial log integral = Gamma(5/2) zeta(3/2)/2.
    A=af(55200)/pi**3;J=3*pi.sqrt()*af(F(3,2)).zeta()/8;K=1600*pi*power32(A)*J
    conditional=power32((power23(af(F(9,16000)))+power23(K*af(F(1,2**26))**2)).upper());conditional_ok=conditional<af(TARGET)
    cones=prepared();runs=[]
    row=run_case(False,boxes,corners,cap,cones,pi,K);runs.append(row);print(json.dumps({k:v for k,v in row.items() if k not in ('remaining_boxes','finite_time_intervals')}),flush=True)
    if not row['target_certified']:
        row=run_case(True,boxes,corners,cap,cones,pi,K);runs.append(row);print(json.dumps({k:v for k,v in row.items() if k not in ('remaining_boxes','finite_time_intervals')}),flush=True)
    result={'status':'target certified' if runs[-1]['target_certified'] else 'finite diagnostic incomplete; bounds not sharp enough','fixed_tau':'1/3','I_target_exact':str(TARGET),'space_domain':'single positive-face quadrant [0,1]^2, x1=1','discarded_corner_boxes':list(map(lambda b:list(map(str,b)),removed)),'cap_radius_exact':str(RHO),'cap_conversion':'030 all-face48 cap upper /24 = single quadrant upper','cap_I_upper_exact':str(capq),'short_time_A_ball':str(A),'short_time_radial_integral_J_ball':str(J),'short_time_global_K_ball':str(K),'conditional_H_target_exact':'9/16000','conditional_epsilon_exact':'1/67108864','conditional_I_upper_ball':str(conditional),'conditional_target_implication_certified':bool(conditional_ok),'runs':runs,'source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'shared_dependencies':'author cb11 cone representation; audited030 pair bound/cap expressions; Python Fraction, Arb/Acb128bit; no independent formula/interval audit yet','no_cost_or_minimum_search':True}
    (OUT/'cb11_integral_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('conditional',str(conditional),'conditional passes',conditional_ok,flush=True)
