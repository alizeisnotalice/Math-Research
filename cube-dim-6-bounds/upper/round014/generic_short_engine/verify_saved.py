"""Saved endpoints/source verification only:0native or full reference calls."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,sys
from flint import arb,ctx,fmpq
sys.set_int_max_str_digits(0);ctx.prec=1024
HERE=Path(__file__).resolve().parent
valid=json.loads((HERE/'validation_results.json').read_text());native=json.loads((HERE/'native_result.json').read_text())
for path,sha in valid['input_sha256'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha

def fb(row):
    f=Q(2)**int(row['common_binary_exponent']);return[Q(int(row[k]))*f for k in ['lower_mantissa','upper_mantissa']]
def B(v):return arb(fmpq(v.numerator,v.denominator))
def fullball(row):
    l,u=fb(row);return arb(B((l+u)/2),B((u-l)/2))
anchors=json.loads((HERE/'anchor_validation.json').read_text());an=json.loads((HERE/'anchors_native.json').read_text())
for row,n in zip(anchors['cases'],an['anchors']):
    if row.get('zero_rejected'):assert n=={'unsupported_zero_node':True};continue
    l,u=fb(n['value']);rl,ru=fb(row['phi_DD5_reference_dyadic']);assert l<=rl<=ru<=u
    for nr,rr in zip(n['partials'],row['partial_DD6_reference_dyadic']):
        l,u=fb(nr);rl,ru=fb(rr);assert l<=rl<=ru<=u
    l,u=fb(n['p_cone']);p=row['p_DD5_equals_complexprod_exact'];exact=Q(int(p['numerator']),int(p['denominator']));assert l<=exact<=u
paths=valid['input_sha256'];GFILE=next(Path(v)for v in paths if v.endswith('exact_local_reference_results.json'));HFILE=next(Path(v)for v in paths if v.endswith('g_time_moments/results.json'));DGFILE=next(Path(v)for v in paths if v.endswith('gradient_exact_reference.json'));DHFILE=next(Path(v)for v in paths if v.endswith('short_minimum_diagnostic/gradient_results.json'))
g=json.loads(GFILE.read_text());h=json.loads(HFILE.read_text());dg=json.loads(DGFILE.read_text());dh=json.loads(DHFILE.read_text());hr=next(v for v in h['time_nodes']if v['tau']=='4/125');dhr=next(v for v in dh['gradient_passes']if v['tau']=='4/125')
G=arb(fmpq(int(g['Q_over_pi5_numerator']),int(g['Q_over_pi5_denominator'])))/arb.pi()**5;DG=arb(fmpq(int(dg['gradient_Q_over_pi5_numerator']),int(dg['gradient_Q_over_pi5_denominator'])))/arb.pi()**5
Kref=G-arb(hr['H_arb384']);Dref=DG-arb(dhr['gradient_x2_to_x6'][0]);l,u=fb(native['K']);assert B(l)>Kref.lower()and B(u)<Kref.upper();l,u=fb(native['gradient'][1]);assert B(l)>Dref.lower()and B(u)<Dref.upper()
e=5*fullball(native['K'])+arb(4)/125*fullball(native['p_T'])+sum((B(Q(v,g['denominator']))*fullball(d)for v,d in zip(g['point_numerators'],native['gradient'])),arb(0));assert e.contains(0)
assert native['full_group_terms']==967680 and valid['full_evaluator_executions']==1
print('PASS6anchors+K+grad_x2+pidentity+Euler;0kernel calls')
