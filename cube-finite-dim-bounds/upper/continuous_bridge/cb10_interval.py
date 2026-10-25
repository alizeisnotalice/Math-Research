"""Two fixed Arb/partition levels, true mu_c(e1), no nested quadrature."""
from pathlib import Path
from fractions import Fraction as F
from math import floor,ceil
import json,time,hashlib
from flint import arb,arb_series,ctx
OUT=Path(__file__).resolve().parent

def af(q):
    q=F(q);return arb(q.numerator)/q.denominator

def phi(r,sign,c,pi):
    Q=(-r*r*(3+sign*r)/(2+sign*r)**2*(r/2).log()+(sign*r-1)*(1+sign*r/2).log()+sign*r/(2+sign*r))/(8*pi**3)
    return c-Q

def blocks(delta,end):
    a=delta
    while a<end:
        b=min(2*a,end);yield a,b;a=b

def branch(sign,end,delta,N,c,pi):
    total=arb(0);positive_error=arb(0);switch_budget=arb(0);records=[];counts={'positive':0,'negative':0,'switch':0};point_radius=arb(0)
    for a,b in blocks(delta,end):
        bc={'a':str(a),'b':str(b),'panels':N,'positive':0,'negative':0,'switch':0};hs=(b-a)/N;h=af(hs)
        for j in range(N):
            lo=a+j*hs;hi=lo+hs;mid=af((lo+hi)/2);box=arb(mid,af(hs/2))
            jet=phi(arb_series([box,1,0],prec=3),sign,c,pi);ran=jet[0]
            assert ran.is_finite()
            if ran>0:
                pm=phi(mid,sign,c,pi);err=2*jet[2].abs_upper()*h**3/24
                panel=pm*h+arb(0,err);positive_error+=err;kind='positive';point_radius+=pm.rad()*h
            elif ran<0:panel=arb(0);kind='negative'
            else:
                upper=ran.upper()
                if upper<=0:panel=arb(0)
                else:
                    cap=h*upper;panel=cap/2+arb(0,cap/2);switch_budget+=cap
                kind='switch'
            total+=panel;counts[kind]+=1;bc[kind]+=1
        records.append(bc)
    # Both signs have the same conservative analytic integral-Q bound at zero.
    d=af(delta)
    tailB=((3+d)/(2-d)**2*d**3/3*((2/d).log()+af(F(1,3)))+(2+d)/(2-d)*d*d/2)/(8*pi**3)
    tail=c*d-tailB/2+arb(0,tailB.abs_upper()/2)
    assert c*d-tailB>0
    return total+tail,{'sign':sign,'endpoint':str(end),'counts':counts,'dyadic_blocks':records,'finite_part_ball':str(total),'tail_integral_ball':str(tail),'tail_Q_integral_upper_ball':str(tailB),'positive_midpoint_remainder_sum_ball':str(positive_error),'switch_range_integral_upper_budget_ball':str(switch_budget),'midpoint_roundoff_radius_sum_ball':str(point_radius)}

def runlevel(bits,N,delta):
    ctx.prec=bits;ctx.cap=3;start=time.monotonic();pi=arb.pi();T=4*arb(2).log()-1;c=T/(24*pi**3)
    # Reuse the separate Jensen lemma.  |e+s v|^2 <=11/3 on [-4/3,2/3].
    assert T<af(F(45,22)) and c>0 and af(F(3,8))>T/24
    plus,pr=branch(1,F(3,2),delta,N,c,pi);minus,mr=branch(-1,F(3,4),delta,N,c,pi)
    mu=4*(plus+minus)/pi;assert mu.is_finite() and mu>0
    den=10**12;lo=F(floor(float(mu.lower())*den)-2,den);hi=F(ceil(float(mu.upper())*den)+2,den)
    assert af(lo)<mu and mu<af(hi)
    return {'bits':bits,'panels_per_dyadic_block':N,'delta_exact':str(delta),'c_ball':str(c),'mu_ball':str(mu),'strict_rational_lower':str(lo),'strict_rational_upper':str(hi),'rational_interval_width':str(hi-lo),'display_lower':float(lo),'display_upper':float(hi),'branches':[pr,mr],'seconds':time.monotonic()-start,'quadrature_justification':'positive panels: midpoint with interval sup |phi_second| h^3/24; negative zero; switches paid full nonnegative range; inverse infinity-tail bound explicit; Arb contains true log/rational G'}
if __name__=='__main__':
    rows=[]
    for args in ((128,32,F(1,256)),(192,128,F(1,4096))):
        row=runlevel(*args);rows.append(row);print(json.dumps({k:v for k,v in row.items() if k!='branches'}),flush=True)
    result={'status':'PASS two fixed rigorous Arb levels','point':[1,0,0],'fixed_c':'(4log2-1)/(24pi^3)','mu_definition':'(4/pi) integral_R (c-G(e1+s(1,1,1)))_+/s^2 ds','levels':rows,'source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'only this density at fixed e1,c; no global source inequality or new C_n bound','dependencies':'special-curve residue/partial fraction identity; reused Jensen lemma for central zero region; Arb ball arithmetic and Taylor remainder theorem'}
    (OUT/'cb10_interval_results.json').write_text(json.dumps(result,indent=2)+'\n')
