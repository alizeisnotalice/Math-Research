#!/usr/bin/env python3
"""R015 independent ledger/endpoint audit. Read-only inputs; JSON to stdout only.
No author module is imported. No phase, ray projection, DD or logarithm is evaluated.
"""
import argparse, ast, difflib, hashlib, json, math
from pathlib import Path
from fractions import Fraction as F
from flint import arb, fmpq, ctx

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text())
def dy(v): return F(int(v[0])) * F(2)**int(v[1])
def interval(v): return dy(v['lower']), dy(v['upper'])
def af(v): return arb(fmpq(v.numerator,v.denominator))
def exact_endpoint(v):
    m,e=v.man_exp(); return F(int(m))*F(2)**int(e)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for x in ('numeric','r012','point_theory','repository','extensions'): ap.add_argument('--'+x.replace('_','-'),required=True,type=Path)
    a=ap.parse_args(); p=a.numeric; oldp=a.r012; repo=a.repository
    sources={}
    def locked(path):
        sources[str(path)]=sha(path); return sources[str(path)]
    for base in (p,oldp):
        m=load(base/'manifest.json'); locked(base/'manifest.json')
        for name,h in m['files'].items(): assert locked(base/name)==h,(base,name)
    pm=a.point_theory/'manifest.sha256'; locked(pm)
    for line in pm.read_text().splitlines():
        h,n=line.split(maxsplit=1); assert locked(a.point_theory/n.lstrip('*'))==h
    r001=repo/'cube-dim-7-bounds/upper/round001/theory/proof.md'
    for line in (p/'prerun-sha256.txt').read_text().splitlines():
        h,n=line.split('  ',1); f=Path(n)
        if not f.is_absolute(): dest=p/f
        elif f.name=='extensions.md': dest=a.extensions
        elif f.name=='proof.md': dest=r001
        else: dest=oldp/f.name
        assert locked(dest)==h,n
    om=load(p/'r012-manifest.json'); assert sha(p/'r012-manifest.json')==sha(oldp/'manifest.json')
    for new,old in [('r012-result.json','result.json'),('cells.jsonl','cells.jsonl'),('r009-rays.json','r009-rays.json'),('r009-inputs.json','r009-inputs.json'),('r010-manifest.json','r010-manifest.json'),('testcases.json','testcases.json')]:
        assert sha(p/new)==om['files'][old]==sha(oldp/old)

    # An independent AST comparison isolates exactly the three changed function bodies.
    oldtext=(oldp/'kernel.py').read_text(); newtext=(p/'green_kernel.py').read_text()
    oldtree=ast.parse(oldtext); newtree=ast.parse(newtext)
    changed=[]
    assert len(oldtree.body)==len(newtree.body)
    for u,v in zip(oldtree.body,newtree.body):
        if ast.dump(u)!=ast.dump(v):
            assert isinstance(u,ast.FunctionDef) and isinstance(v,ast.FunctionDef) and u.name==v.name
            changed.append(u.name)
    assert changed==['dd_real','phase','generic_chunk'],changed
    hunks=list(difflib.unified_diff(oldtext.splitlines(),newtext.splitlines(),n=1))
    assert sum(s.startswith('@@') for s in hunks)==3
    assert "assert key[0]!=0;out=A(abs(key[0])).log()" in newtext
    assert "out=A(Fraction((-1)**(k-1),k*key[0]**k))" in newtext
    assert "assert eps==0;B=dx*L;beta=0" in newtext
    assert "assert cnt['zero_ray_knots']==0" in newtext
    assert 'phase(n,dx,Fraction(0),64*chunk+signid)' in newtext
    assert 'out=2*total/(2*arb.pi())**7' in newtext
    # Scope is actual run.py generic path. Retained old scalar/origin helpers are unused.
    run=(p/'run.py').read_text()
    assert 'kernel.generic_chunk(j,256)' in run and 'pending=iter(range(4,5040))' in run
    assert "assert not (D/'terminal.jsonl').exists()" in run
    assert 'for j in range(4)' in run and 'no_extra_normalization=True' in run

    r=load(p/'result.json'); old=load(oldp/'result.json'); assessment=load(p/'assessment.json'); pilot=load(p/'pilot.json')
    xhex=r['point_floathex']; assert xhex==old['point_floathex']==load(p/'testcases.json')['generic_case']['x_floathex']
    r003=repo/'cube-dim-7-bounds/upper/round003/truncated_numeric/results.json'; locked(r003)
    candidates=[v for v in load(r003)['minima'] if v['lowest_sobol_index']==2710]
    assert candidates and all(v['point_floathex']==xhex for v in candidates)
    x=[F.from_float(float.fromhex(h)) for h in xhex]
    assert x[0]==1 and all(0<v<1 for v in x[1:]) and len(set(map(abs,x)))==7
    dx=math.lcm(*(v.denominator for v in x)); nums=[int(v*dx) for v in x]
    rows=[json.loads(s) for s in (p/'terminal.jsonl').read_text().splitlines()]
    oldrows=[json.loads(s) for s in (oldp/'terminal.jsonl').read_text().splitlines()]
    assert len(oldrows)==5040 and {v['chunk'] for v in oldrows}==set(range(5040))
    assert all(v['status']=='completed' and v['counts']['zero_ray_knots']==0 for v in oldrows)
    assert old['halfphases']==322560 and old['counts']['zero_ray_knots']==0
    assert len(rows)==5040 and {v['chunk'] for v in rows}==set(range(5040))
    assert [v['chunk'] for v in rows[:4]]==list(range(4)) and pilot['status']=='passed' and pilot['chunks']==list(range(4))
    lo=F(); hi=F(); counts={k:0 for k in r['counts']}; negative=0; cross=0; pilotlo=F(); pilothi=F()
    for v in rows:
        j=v['chunk']; assert v['status']=='completed' and v['precision']==256
        assert v['start_phase']==64*j and v['end_phase']==64*(j+1)
        assert v['counts']['top_requests']==64*212 and v['counts']['zero_ray_knots']==0
        assert v['cpu_seconds']>0 and v['wall_seconds']>0 and 'reused' not in v
        l,u=interval(v['interval']); assert l<=u; lo+=l; hi+=u
        negative+=int(u<0); cross+=int(l<=0<=u)
        if j<4: pilotlo+=l; pilothi+=u
        for k in counts: counts[k]+=v['counts'][k]
    assert counts==r['counts'] and counts['top_requests']==68382720
    assert counts['unique_subdd']==659175346 and counts['log_evaluations']==14837760
    assert counts['hermite_evaluations']==counts['repeated_top_nodes']==0
    pl,pu=interval(pilot['interval']); assert pl<=pilotlo<=pilothi<=pu
    assert pu-pl<=F(1,10**12)*F(256,322560)
    assert r['status']=='completed' and r['epsilon']==0 and r['halfphases']==322560
    assert r['chunks']==5040 and r['pilot_chunks']==4 and r['precision']==256 and r['no_extra_normalization']
    gl,gu=interval(r['full_G']); assert 0<gl<=lo<=hi<=gu and gu-gl<=F(1,10**12)
    assert dy(r['width']['upper'])>=gu-gl
    # Only scalar outward conversions at768bits, starting from the independent exact ledger sums.
    ctx.prec=768
    il=exact_endpoint(af(lo).lower()); iu=exact_endpoint(af(hi).upper())
    assert gl<=il<=lo<=hi<=iu<=gu
    hu=dy(old['full_H']['upper']); assert hu<il
    thresholds={str(k):F(1,k)-iu for k in (219000,194000,179000)}
    assert all(v>0 for v in thresholds.values())
    assert F(assessment['G_upper_exact'])==gu
    for k in (219000,194000,179000):
        assert assessment['thresholds'][str(k)]['strict'] and gu<F(1,k)
    I0=F(7723627498344,10**18); cost6=7**7*I0**6/gu**7
    assert F(assessment['fixed_I0'])==I0 and F(assessment['fixed_cost_power6_lower_exact'])==cost6
    assert cost6>128**6 and assessment['fixed_cost_power6_strictly_gt_128_power6']
    out={'status':'PASS','kind':'independent exact signed-ledger and fixed-scalar audit; not a DD rerun',
         'precision_bits':768,'source_changed_functions':changed,'source_diff_hunks':3,
         'chunks':len(rows),'halfphases':322560,'new_G_chunks':5040,'reused_H_chunks':0,
         'negative_chunk_intervals':negative,'zero_crossing_chunk_intervals':cross,'counts':counts,
         'point_floathex':xhex,'point_integer_numerators':nums,'point_denominator':dx,
         'G_independent_lower':str(il),'G_independent_upper':str(iu),
         'G_exact_ledger_lower':str(lo),'G_exact_ledger_upper':str(hi),
         '768bit_conversion_exact':il==lo and iu==hi,'G_reported_lower':str(gl),'G_reported_upper':str(gu),
         'G_exact_ledger_width':str(hi-lo),'G_reported_width':str(gu-gl),'G_minus_H_gap_lower':str(il-hu),
         'threshold_positive_gaps_using_independent_upper':{k:str(v) for k,v in thresholds.items()},
         'fixed_I0':str(I0),'fixed_certificate_cost_power6_lower_using_reported_G_upper':str(cost6),
         'fixed_certificate_cost_power6_minus_128_power6':str(cost6-128**6),
         'new_phases':0,'DD_evaluations':0,'new_kernel_points':0,'MATLAB_calls':0,'MC_samples':0,
         'scope':'Three original-G essential floor targets excluded; fixed-I0 certificate expression exceeds128. No actualI/M/C7 lower bound, new globalg, newI, or whole-Green-model exclusion.',
         'source_sha256':sources}
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
