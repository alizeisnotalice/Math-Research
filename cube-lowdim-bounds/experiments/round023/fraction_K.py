"""Independent exact rational B-splines and integer-root bounds; no Arb imports.
Budget fixed before run: seven existing K points times two certificates, 14
values, 60-second guarded run; no point search or precision escalation.
"""
from pathlib import Path
from fractions import Fraction as Q
from time import monotonic
import json,hashlib
ROOT=Path(__file__).resolve().parent
START=monotonic();D=10**70

def guard():
    if monotonic()-START>60:raise TimeoutError('60-second budget')
def root40(q):
    z=q*D**40;n=z.numerator//z.denominator
    lo,hi=0,1<<((n.bit_length()+39)//40)
    while lo+1<hi:
        guard();m=(lo+hi)//2
        if m**40<=n:lo=m
        else:hi=m
    assert Q(lo,D)**40<=q<Q(lo+1,D)**40
    return Q(lo,D),Q(lo+1,D)
def B(x):
    x=abs(x)
    if x>=2:return Q(0)
    if x<=1:return Q(2,3)-x*x+x**3/2
    return (2-x)**3/6
def orb(i,j):
    return {(a,b) for a in (i,-i) for b in (j,-j)}|{(b,a) for a in (i,-i) for b in (j,-j)}
def run():
    out=ROOT/'fraction_K_result.json'
    if out.exists():raise RuntimeError('Do not overwrite or repeat')
    points=[p for p in json.loads((ROOT/'points.json').read_text())['points'] if p['quantity']=='K'];assert len(points)==7
    values=[]
    for label,path in [('K1',ROOT.parent/'round01/scaled_certificate.json'),('K3',ROOT.parent/'round03/candidate.json')]:
        s=json.loads(path.read_text());a,R,h,c=(Q(s[k]) for k in ('alpha','R','h','c'));assert a==Q(47,40) and c>0
        rlo,rhi=root40(R**(-47))
        for p in points:
            guard();u,v=Q(p['u']),Q(p['v']);r=u+v;assert 0<r<=1<R
            spline=Q(0);terms=0
            for key,coef in s['a'].items():
                i,j=map(int,key.split(','))
                for x,y in orb(i,j):
                    bx=B(u/h-x)
                    if bx:
                        by=B(v/h-y)
                        if by:spline+=Q(coef)*bx*by;terms+=1
            lo,hi=root40(r**(-47));lower=c*(lo-rhi)+spline;upper=c*(hi-rlo)+spline
            state='violates_K_ge_1' if upper<1 else 'K_ge_1_proved_at_point' if lower>=1 else 'undetermined'
            values.append({'kernel':label,'point_id':p['id'],'u':p['u'],'v':p['v'],'spline_exact':str(spline),'active_terms':terms,'lower':str(lower),'upper':str(upper),'width':str(upper-lower),'status':state,'lower_display':float(lower),'upper_display':float(upper)})
    out.write_text(json.dumps({'status':'complete','values':values,'K_evaluations':len(values),'root_decimal_scale':str(D),'seconds':monotonic()-START,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Rigorous pointwise K bounds only; no AK or global proof.'},indent=2)+'\n')
    print([(x['kernel'],x['point_id'],x['status'],x['lower_display']) for x in values])
if __name__=='__main__':run()
