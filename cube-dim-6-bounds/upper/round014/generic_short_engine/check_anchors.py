from pathlib import Path
from fractions import Fraction as Q
from flint import arb,acb,ctx,fmpq
import json,hashlib
HERE=Path(__file__).resolve().parent
ctx.prec=640
reg=json.loads((HERE/'registration.json').read_text());native=json.loads((HERE/'anchors_native.json').read_text());T=arb(4)/125;t=fmpq(4,125);ii=acb(0,1)

def exact_bounds(row):
    e=int(row['common_binary_exponent']);return[Q(int(row[k]))*Q(2)**e for k in ['lower_mantissa','upper_mantissa']]
def dyadic_ref(v):
    lm,le=map(int,v.lower().man_exp());um,ue=map(int,v.upper().man_exp());e=min(le,ue)
    return {'lower_mantissa':str(lm<<(le-e)),'upper_mantissa':str(um<<(ue-e)),'common_binary_exponent':str(e)}
def enclose(v,row):
    record=dyadic_ref(v);lo,hi=exact_bounds(row);rl,ru=exact_bounds(record);assert lo<=rl<=ru<=hi
    return record
def jet(xq,k):
    x=arb(xq.numerator)/xq.denominator
    if k==0:return(arb.pi()/2 if xq>0 else-arb.pi()/2)-(x/T).atan()
    old=((1 if k&1 else-1)*(ii**k/acb(T,x)**k-(-ii)**k/acb(T,-x)**k)/(2*k*ii)).real
    return-old
def dd_phi(z):
    z=sorted(z);b=[jet(v,0)for v in z]
    for k in range(1,len(z)):
        b=[jet(z[j],k)if z[j]==z[j+k]else(b[j+1]-b[j])/(arb(z[j+k].numerator)/z[j+k].denominator-arb(z[j].numerator)/z[j].denominator)for j in range(len(z)-k)]
    return b[0]
def cmul(a,b):return(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def pole(x):
    D=t*t+x*x;return(t/D,x/D)
def rational_jet(x,k):
    p=(fmpq(1),fmpq(0));r=pole(x)
    for _ in range(k+1):p=cmul(p,r)
    return[p[1],p[0],-p[1],-p[0]][k%4]
def dd_rational(z):
    z=sorted(z);b=[rational_jet(v,0)for v in z]
    for k in range(1,len(z)):b=[rational_jet(z[j],k)if z[j]==z[j+k]else(b[j+1]-b[j])/(z[j+k]-z[j])for j in range(len(z)-k)]
    return b[0]

rows=[]
for i,(a,n)in enumerate(zip(reg['anchors'],native['anchors'])):
    if any(v[0]==0 for v in a):assert n=={'unsupported_zero_node':True};rows.append({'index':i,'zero_rejected':True});continue
    z=[Q(v[0],v[1])for v in a]
    row={'index':i,'phi_DD5_reference_dyadic':enclose(dd_phi(z),n['value']),'partial_DD6_reference_dyadic':[]}
    for j in range(6):row['partial_DD6_reference_dyadic'].append(enclose(dd_phi(z+[z[j]]),n['partials'][j]))
    zr=[fmpq(v.numerator,v.denominator)for v in z];p=dd_rational(zr);prod=(fmpq(1),fmpq(0))
    for x in zr:prod=cmul(prod,pole(x))
    assert p==prod[0]
    lo,hi=exact_bounds(n['p_cone']);q=Q(int(p.numerator),int(p.denominator));assert lo<=q<=hi
    row['p_DD5_equals_complexprod_exact']={'numerator':str(q.numerator),'denominator':str(q.denominator)};rows.append(row)
result={'status':'PASS_ALL6_PHI_DD_PARTIAL_AND_P_IDENTITY_CASES','T':'4/125','reference_precision_bits':640,'cases':rows,'native_source_sha256':hashlib.sha256((HERE/'generic_short.cpp').read_bytes()).hexdigest(),'reference_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'anchor_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'])
