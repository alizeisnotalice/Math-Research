#!/usr/bin/env python3
"""Independent R016 exact ledger and ten 768-bit scalars. Input-only; JSON stdout.
No author module imported; no phase, kernel, spectral enumeration, sample or mesh.
Alternate path: positive roots, integer factorials and rational powers of pi.
"""
import argparse, hashlib, json, math
from pathlib import Path
from fractions import Fraction as F
from flint import arb, fmpq, ctx

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text())
def A(v):
    v=F(v); return arb(fmpq(v.numerator,v.denominator))
def dy(v): return F(int(v[0]))*F(2)**int(v[1])
def endpoint(x):
    m,e=x.man_exp(); return F(int(m))*F(2)**int(e)
def bounds(x): return endpoint(x.lower()),endpoint(x.upper())
def encode(x):
    l,u=bounds(x); return {'lower':str(l),'upper':str(u),'display':x.str(50)}
def ceilgrid(v):
    n=(v.numerator*10**18+v.denominator-1)//v.denominator
    return F(n,10**18)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for s in ('numeric','theory','primary','sharp','catalog','pilot'): ap.add_argument('--'+s,required=True,type=Path)
    args=ap.parse_args(); p=args.numeric; sources={}
    def lock(q): sources[str(q)]=sha(q); return sources[str(q)]
    lock(p/'manifest.json')
    for n,h in load(p/'manifest.json')['files'].items(): assert lock(p/n)==h,n
    lock(args.theory/'manifest.sha256')
    for line in (args.theory/'manifest.sha256').read_text().splitlines():
        h,n=line.split(maxsplit=1); assert lock(args.theory/n.lstrip('*'))==h,n
    for line in (p/'prerun-sha256.txt').read_text().splitlines():
        h,n=line.split('  ',1); q=Path(n)
        if not q.is_absolute(): target=p/q
        elif q.name=='catalog.json': target=args.catalog
        elif q.name==args.primary.name: target=args.primary
        elif q.name==args.sharp.name: target=args.sharp
        else: target=args.theory/q.name
        assert lock(target)==h,n
    for f,m,k in [('r010-cells.jsonl','r010-manifest.json','cells.jsonl'),('r006-certificate.json','r006-manifest.json','certificate.json')]:
        assert sha(p/f)==load(p/m)['files'][k]
    primary=load(args.primary); sharp=load(args.sharp)
    assert lock(args.pilot)==primary['source_sha256']
    assert lock(args.catalog)==sharp['source_sha256']
    assert sharp['max_squared_unitball_norm']=='1/2048' and sharp['psi_euclidean_lower_squared']==2048

    rows=[json.loads(s) for s in (p/'r010-cells.jsonl').read_text().splitlines()]
    assert len(rows)==135
    seen=set(); total=F()
    for row in rows:
        for item in row['simplices']:
            vertices=tuple(item['vertices']); assert len(vertices)==len(set(vertices))==7 and vertices not in seen
            seen.add(vertices); assert len(item['q'])==7 and all(v>0 for v in item['q'])
            weight=F(abs(item['detW']),math.prod(item['q'])); assert weight==F(item['alpha']) and weight>0
            total+=weight
    r=load(p/'result.json'); assert r['status']=='completed' and len(seen)==212
    assert str(total)==r['A_sum']==primary['A_sum']
    assert total==F(907157967345860903220972676779187,28127466839029475035590994122592223232000000000)
    assert 0<5040*total<=F(1,512)

    # Gamma(65/2)/Gamma(32)=C sqrt(pi). Half-integer formula is exact.
    C=F(math.factorial(64),4**32*math.factorial(32)*math.factorial(31))
    C7=F(math.factorial(70),4**35*math.factorial(35)*math.factorial(31))
    b=F(65*67*69,8); assert C7==C*b
    # Normalization exponents: raw I sixth power=(5040 A)/(840^6 pi^10).
    ctx.prec=768; pi=arb.pi(); pi3=pi**3
    p0=A(5040*total)/pi**7
    psharp=3*arb(2).sqrt()/(2**38*pi**4)
    mraw=(128*pi3*pi.sqrt()*p0).root(7)
    raw=(A(5040*total)/(840**6*pi**10)).root(6)
    sharpI=(psharp/(840**6*pi3)).root(6)
    eta=(A(C**6/b)*pi3).root(6)
    rb=(A(5040*total*C**6/b)/(840**6*pi**7)).root(6)
    rawup=ceilgrid(bounds(raw)[1]); rbup=ceilgrid(bounds(rb)[1]); g=F(1,320000)
    req=A(7**7*(rbup/128)**6).root(7)
    rawcost6=7**7*rawup**6/g**7; rbcost6=7**7*rbup**6/g**7
    rawcost=A(rawcost6).root(6); rbcost=A(rbcost6).root(6)
    values={'p0_fan':p0,'p0_sharp_upper':psharp,'M_raw_bound':mraw,'B_fan':raw,'B_sharp':sharpI,
            'eta':eta,'B_RB':rb,'g_required_B_RB':req,'conditional_cost_B_RB':rbcost,'conditional_cost_B_fan':rawcost}
    assert len(values)==10 and set(values)==set(r['enclosures'])
    containment={}
    for n,v in values.items():
        l,u=bounds(v); old=r['enclosures'][n]; ol,ou=dy(old['lower']),dy(old['upper'])
        assert 0<l<=u and ol<=l<=u<=ou,n
        containment[n]=True
    assert bounds(p0)[1]<bounds(psharp)[0]
    assert bounds(eta)[1]<1 and bounds(rb)[1]<bounds(raw)[0]
    assert (rb-eta*raw).contains(0)
    assert rawup==F(r['B_fan_rational_upper'])==F(4127002498830,10**18)
    assert rbup==F(r['B_RB_rational_upper'])==F(4018421141119,10**18)
    assert bounds(raw)[1]<rawup and bounds(rb)[1]<rbup
    old=load(p/'r006-certificate.json'); oldlower=dy(old['I_upper']['lower'])
    assert oldlower==dy(r['R006_deterministic_budget_expression_lower'])
    assert rbup<rawup<oldlower and bounds(sharpI)[1]<F(4940777144766,10**18)<oldlower
    assert bounds(req)[1]<g and rbcost6<rawcost6<128**6
    assert rawcost6==F(r['conditional_cost6_raw_exact']) and rbcost6==F(r['conditional_cost6_RB_exact'])
    assert not r['global_g_proved'] and not r['new_C7_bound']
    primary_raw=[F(v) for v in primary['I_bound']]
    assert primary_raw[0]<=bounds(raw)[0]<=bounds(raw)[1]<=primary_raw[1]
    print(json.dumps({'status':'PASS','precision_bits':768,'named_scalar_enclosures':10,
        'independent_path':'positive real root operators and exact half-integer factorial ratio; no log/exp or Gamma call',
        'A_sum':str(total),'cells':135,'simplices':212,'gamma_ratio_rational_coefficient_C':str(C),
        'gamma_seventh_recurrence_b':str(b),'Gamma71over32_equals_C_b_sqrtpi_exact':True,
        'enclosures':{n:encode(v) for n,v in values.items()},'contained_in_frozen_512bit_enclosures':containment,
        'B_fan_rational_upper':str(rawup),'B_RB_rational_upper':str(rbup),
        'raw_budget_gap_to_R006_expression_lower':str(oldlower-rawup),
        'R006_lower_endpoint_scope':'lower endpoint of the old deterministic UPPER-BUDGET expression; not actualI lower',
        'conditional_g':str(g),'conditional_cost6_raw_exact':str(rawcost6),'conditional_cost6_RB_exact':str(rbcost6),
        'conditional_raw_gap_to_128_power6':str(128**6-rawcost6),'conditional_RB_gap_to_128_power6':str(128**6-rbcost6),
        'new_phases':0,'DD_evaluations':0,'MC_samples':0,'MATLAB_calls':0,'source_sha256':sources,
        'scope':'New complete original-I budget only. Conditional g=1/320000 remains unproved. No new C7 bound or actualI/M lower.'},indent=2))

if __name__=='__main__': main()
