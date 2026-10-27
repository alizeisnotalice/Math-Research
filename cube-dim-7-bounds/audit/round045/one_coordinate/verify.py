"""Independent local R045 identities and new Arb derivative-point checks.

This does not inspect the whole proposed Cn theorem or replay Terra's quadrature.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,sys
import sympy as sp
from flint import arb,fmpq,ctx

HERE=Path(__file__).resolve().parent
R,s,g,q,c,e,t,w=sp.symbols('R s g q c e t w',positive=True)
def FF(R,s):
    return -2*R*sp.log(1-s/(2*R))-s
F=FF(R,s)
assert sp.simplify(sp.diff(F,s)-s/(2*R-s))==0
assert sp.simplify(sp.diff(F,R)+2*sp.log(1-s/(2*R))+2*s/(2*R-s))==0
primitive=g*g/2-(2*R-g)*FF(R,g)
assert sp.simplify(sp.diff(primitive,g)-FF(R,g))==0
assert sp.simplify(primitive.subs(g,0))==0
Rq=q+c;gq=c+e*q;aq=(1-e)*q;rhoq=t*q
central=1+FF(Rq,gq)/rhoq
edge=1+(FF(Rq,gq)+FF(Rq,Rq-w)-(w-aq))/(2*rhoq)
outside=1+(FF(Rq,gq)-(w-aq))/(2*rhoq)
assert sp.simplify((edge-central).subs(w,aq))==0
jump=(1-e)/(2*rhoq*(1-gq/(2*Rq)))
assert sp.simplify((sp.diff(edge,q)-sp.diff(central,q)).subs(w,aq)-jump)==0
assert sp.simplify((edge-outside).subs(w,Rq))==0
assert sp.simplify((sp.diff(edge,q)-sp.diff(outside,q)).subs(w,Rq))==0

# Exact endpoint arithmetic for conservative estimates in the analytic lemma.
dd=Q(1,16)
assert 2/(1-dd)<4
assert 1+2*dd+10*dd*dd<2
assert 12*dd*dd<2
assert 11*dd/4<1
assert 1280*Q(1,2**20)==Q(5,4096)

ctx.prec=512
def ab(x):
    return arb(fmpq(x.numerator,x.denominator))
def fv(R,s):
    return -2*ab(R)*(1-ab(s/(2*R))).log()-ab(s)
def fs(R,s):
    return ab(s/(2*R-s))
def fr(R,s):
    return -2*(1-ab(s/(2*R))).log()-ab(2*s/(2*R-s))
d=Q(1,32);h=d/16;eps=h/4;tt=Q(128);r=1+d/2
q0=r-d+h;qq=q0-h;RR=qq+d-h;gg=d-h+eps*qq;aa=(1-eps)*qq;rho=tt*qq
assert qq>=Q(1,2) and RR>=Q(1,2) and rho>=4*RR and 0<gg<=d
assert h<=d/8 and eps<=d and 0<=r-RR<=2*h
F0=fv(RR,gg);F0p=fr(RR,gg)+ab(eps)*fs(RR,gg)
rows=[]
for ww in (Q(0),r-3*d,r-d,r):
    m=max(Q(0),ww-aa)
    ss=max(Q(0),gg-m)
    B=F0+fv(RR,ss)-ab(m)
    if ww<aa:
        Bp=2*F0p
    elif ww<RR:
        Bp=F0p+fr(RR,ss)+fs(RR,ss)+ab(1-eps)
    else:
        Bp=F0p+ab(1-eps)
    J=1+B/(2*ab(rho))
    Jp=(Bp-B/ab(qq))/(2*ab(rho))
    edge_indicator=int(ww>r-2*d)
    bound=Q(64,rho)*(d*d+edge_indicator)
    assert J>ab(Q(1,2)) and abs(Jp/J)<ab(bound)
    rows.append({'w':str(ww),'edge':edge_indicator,'J':str(J),'abs_log_derivative':str(abs(Jp/J)),'bound':str(bound)})
# Stable nonnegative series certificate for the same new local point.
u=gg/(2*RR)
prefix=2*RR*sum((u**k/k for k in range(2,9)),Q(0))
tail=2*RR*u**9/(9*(1-u))
assert F0>ab(prefix) and F0<ab(prefix+tail)

receipt={'status':'PASS independent local one-coordinate checks',
         'sympy':sp.__version__,'arb_precision_bits':ctx.prec,
         'symbolic_checks':['F partial derivatives','baseline exact mean primitive','capture-onset kink jump','C1 match at |w|=R'],
         'new_point':{'d':str(d),'h':str(h),'epsilon':str(eps),'r':str(r),'q':str(qq),'rho':str(rho),'checks':rows},
         'stable_positive_series_tail_verified':True,'uniform_log_product_shift':'<=5/4096 for h=2^-20/(Ld)',
         'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'scope':'Local closed form, a.e./AC q derivative, baseline moments and product bridge only. No Cn theorem acceptance, production/quadrature replay, n7, MATLAB, or Git mutation.'}
if len(sys.argv)>1:
    source=Path(sys.argv[1]);txt=source.read_text()
    section4=txt.split('## 4. Exact one-coordinate inverse row',1)[1].split('## 5. Typical targets',1)[0]
    bridge=txt.split('## 5. Typical targets',1)[1].split('The full inverse-weighted source shell coefficient',1)[0]
    receipt['source_identity']={'PROOF_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'section4_sha256':hashlib.sha256(section4.encode()).hexdigest(),
        'section5_through_product_bridge_sha256':hashlib.sha256(bridge.encode()).hexdigest()}
(HERE/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
