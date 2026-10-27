#!/usr/bin/env python3
"""Future explicit mathematical check. NOT executed in the byte archive task."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,itertools,os
from io175 import read,put,checked,lines,canonical,sha
from view_recipe import build_views

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify-mathematics',action='store_true')
    ap.add_argument('--packet',type=Path,required=True);ap.add_argument('--restored179',type=Path,required=True)
    ap.add_argument('--repository',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    if not a.verify_mathematics:ap.error('Mathematical verification requires explicit --verify-mathematics;archive mode never executes it.')
    c=read(a.packet/'certificate.json')
    for name,h in c['portable_code_hashes'].items():assert sha(Path(__file__).resolve().parent/name)==h
    a.output.mkdir(exist_ok=False);put(a.output/'registration.json',{'pid':os.getpid(),'explicit_mathematical_verification':True,'certificate_sha256':sha(a.packet/'certificate.json')})
    v=build_views(a.restored179,a.repository,a.output/'views',c['view_bindings']);updates={};N=1<<160;B=1<<64
    for i in range(15):
        pair=c['pairs'][i];p=a.packet/pair['relative_path'];checked({'path':str(p),'bytes':pair['bytes'],'sha256':pair['sha256']})
        digest=hashlib.sha256();count=0
        for r,ends in itertools.zip_longest(lines(v['shards'][i]['path']),lines(p)):
            assert r is not None and ends is not None and len(ends)==2;digest.update(canonical(r)+b'\n');Y,Z=ends
            assert type(Y) is int and type(Z) is int and Y>=0 and Z>=0
            U,L=r['range_U_i'],r['range_L_i'];assert 0<=L<U and r['strict112_active']
            assert Y**4*N<=U*(1<<256)<(Y+1)**4*N and Z**4*U<=L*(1<<256)<(Z+1)**4*U
            z=F(Z,B);den=1+z+z*z+z*z*z;new=F(Y,B)*(den+z**4)/den
            slope=max(F(r['slope164_numerator'],r['slope164_denominator']),new);assert slope>=F(1,16)
            num=r['weight']*(r['mean112_i']-r['mean162_i'])*slope.numerator;divisor=slope.denominator<<(4*r['level'])
            total,rem=divmod(num,divisor);assert 0<=rem<divisor
            value=r['self112']-total;assert r['selfL']<=value<=r['self166'];assert r['key'] not in updates
            updates[r['key']]=value;count+=1
        assert count==pair['count'] and digest.hexdigest()==c['parameter_equivalence']['shards'][i]['canonical_parameters_sha256']
    t=read(v['topology']['path']);assert hashlib.sha256(canonical(t)).hexdigest()==c['parameter_equivalence']['canonical_topology_sha256']
    nodes={r['key']:r for r in t['nodes']};assert len(nodes)==97972 and len(updates)==97971
    assert set(nodes)-set(updates)=={'5:31,31,31,31'};roots=t['roots'];assert len(roots)==len(set(roots))==52360
    assert set(roots)=={k for k,n in nodes.items() if n['parent'] is None};eff={};visiting=set()
    def visit(k):
        if k in eff:return eff[k]
        assert k not in visiting;visiting.add(k);n=nodes[k];sv=updates.get(k,n['self112']);assert n['selfL']<=sv<=n['self166']
        children=n['children'];assert len(children)==len(set(children))
        for ch in children:assert nodes[ch]['parent']==k and nodes[ch]['level']==n['level']+1
        value=min(sv,sum(visit(ch) for ch in children)) if children else sv
        assert value<=n['eff166'];eff[k]=value;visiting.remove(k);return value
    units=sum(visit(k) for k in roots);assert set(eff)==set(nodes);UH=F(units,N);assert str(UH)==c['new_UH']
    put(a.output/'results.json',{'status':'explicit_portable_mathematical_check_pass','new_UH':str(UH),'root_pairs':97971,
        'nodes':97972,'roots':52360,'global_Cn_claim':False,'full_historical_output_JSON_byte_restoration_claim':False})
