#!/usr/bin/env python3
"""Portable R004 verifier. Fresh cheap checks; recorded complete C4, no C4 replay."""
import argparse
import hashlib
import json
import re
import sys
from fractions import Fraction as F
from pathlib import Path

if not __debug__:
    raise SystemExit('Run without -O: mathematical assertions must remain enabled.')
parser=argparse.ArgumentParser()
parser.add_argument('--out',type=Path,help='Write a new receipt outside the frozen package; default prints only.')
args=parser.parse_args()
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
VENDOR=PROJECT/'kernel-rescaling/vendor'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest_path=HERE/'package-manifest.json'
manifest=json.loads(manifest_path.read_text())
for rel,digest in manifest['files'].items():
    p=HERE/rel
    assert not Path(rel).is_absolute() and '..' not in Path(rel).parts
    assert sha(p)==digest
for rel,digest in manifest['reused_project_files'].items():
    p=PROJECT/rel
    assert not Path(rel).is_absolute() and '..' not in Path(rel).parts
    assert sha(p)==digest
if args.out is not None:
    protected={ (HERE/p).resolve() for p in manifest['files'] }
    protected|={ (PROJECT/p).resolve() for p in manifest['reused_project_files'] }
    protected|={manifest_path.resolve(),(HERE/'portable-receipt.json').resolve(),(HERE/'portable-execution.json').resolve()}
    if args.out.resolve() in protected:
        raise SystemExit('Choose a new receipt path; frozen input/output records must not be overwritten.')

sys.path.insert(0,str(VENDOR))
from flint import ctx
from hl_core import Spec, A, pw
from hl_check import check_c12, c4_cells, prov_digest, provenance
ctx.prec=192
base_path=PROJECT/'kernel-rescaling/certificates/scaled.json'
path=HERE/'candidate.json'
assert sha(base_path)=='d570ac3ffc1357e129183eda9579fde8a1bf3775e914efe6ef2d6be336348c99'
assert sha(path)=='88af41536e1117d953bbc4594ee99f1319b79b75d98d740c5d3c82caaf1f306f'
base_d=json.loads(base_path.read_text())
d=json.loads(path.read_text())
assert set(d)==set(base_d)
for key in set(d)-{'a','bound','note'}:
    assert d[key]==base_d[key]
assert set(d['a'])==set(base_d['a']) and len(d['a'])==1181
assert {k:F(d['a'][k])-F(base_d['a'][k]) for k in d['a'] if F(d['a'][k])!=F(base_d['a'][k])}=={'0,0':-F(1,200)}
claim=F(3590865312933,10**12)
assert F(d['bound'])==claim < F(base_d['bound'])
sp=Spec(path)
expected=c4_cells(sp)
assert len(expected)==len(set(expected))==1368
report_path=HERE/'records/full-c1234-report.json'
report=json.loads(report_path.read_text())
assert report['spec_sha256']==sha(path)
assert F(report['claim'])==claim
assert set(report['args']['only'].split(','))=={'c1','c2','c3','c4'}
assert report['args']['cells']==''
assert report['args']['procs']==10 and report['args']['maxdepth']==8 and report['args']['kdepth']==10
assert all(report[k]['ok'] for k in ['c1','c2','c3','c4'])
assert report['c1']['cells']==1368 and report['c1']['leaves']==2025 and report['c1']['fails']==[]
assert report['c2']['cells']==272 and report['c2']['leaves']==272 and report['c2']['fails']==[]
rp=report['provenance']
assert rp['hl_check.py']==sha(VENDOR/'hl_check.py')=='0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8'
assert rp['hl_core.py']==sha(VENDOR/'hl_core.py')=='c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf'
rows=report['c4']['results']
assert len(rows)==report['c4']['cells']==len(expected)
assert {(r['p'],r['q']) for r in rows}==set(expected)
assert report['c4']['prov_consistent']
assert all(r['ok'] and r['nfail']==0 and r['fails']==[] and r['spec_sha']==sha(path)
           and r['prov']==prov_digest(rp) for r in rows)
assert sum(r['leaves'] for r in rows)==report['c4']['leaves']==443757

# Independent stdout cross-check, with both length and set equality.
stdout=(HERE/'records/full-c1234.stdout.log').read_text()
pattern=r'^cell \((\d+),(\d+)\): (ok|FAIL) min .* leaves (\d+) taylor (\d+) depth (\d+) fails (\d+) '
logged=re.findall(pattern,stdout,re.MULTILINE)
assert len(logged)==len(rows)
assert len({(int(p),int(q)) for p,q,*_ in logged})==len(rows)
by_cell={(r['p'],r['q']):r for r in rows}
for p,q,status,leaves,taylor,depth,nfail in logged:
    r=by_cell[int(p),int(q)]
    assert status=='ok' and int(nfail)==0
    assert (int(leaves),int(taylor),int(depth))==(r['leaves'],r['taylor'],r['maxdepth'])
assert 'OVERALL (c1,c2,c3,c4): PASS' in stdout

# Exact integration of the two even spline pieces, independent of the WL log.
integral=lambda c,a,b:sum(v*(b**(k+1)-a**(k+1))/F(k+1) for k,v in enumerate(c))
int_B=2*(integral([F(2,3),0,-1,F(1,2)],F(0),F(1))
          +integral([F(4,3),-2,1,-F(1,6)],F(1),F(2)))
assert int_B==1 and sp.h**2*int_B**2/2==F(1,2048)

# Independently expand all signed/permuted coefficient indices, not orbit weights.
tot=sum(sp.coef(i,j) for i in range(-sp.N,sp.N+1) for j in range(-sp.N,sp.N+1))
old=Spec(base_path)
oldtot=sum(old.coef(i,j) for i in range(-old.N,old.N+1) for j in range(-old.N,old.N+1))
assert oldtot-tot==F(1,200)
assert F(report['c3']['sum_a'])==tot
mass=sp.cA*sp.a_*pw(sp.RA,2-sp.a_)/(2-sp.a_)+A(sp.h**2*tot)/2
lo,hi=F(str(mass.lower().fmpq())),F(str(mass.upper().fmpq()))
assert F('3.590865312932196537') < lo <= hi < F('3.590865312932196538') < claim
checks=check_c12(sp,maxdepth=10)
assert not checks['c1_fail'] and not checks['c2_fail']
assert checks['c1_cells']==1368 and checks['c2_cells']==272
assert checks['c1_leaves']==2025 and checks['c2_leaves']==272
assert checks['c1_min'].lower().fmpq()>=0
c2lo=F(str(checks['c2_min'].lower().fmpq()))
assert c2lo>1

# The large-step negative control excludes only this affine ray.
unit=Spec(HERE/'unit-origin-spline.json')
assert unit.c==0 and unit.acls=={(0,0):F(1)}
assert (unit.alpha,unit.R,unit.h)==(sp.alpha,sp.R,sp.h)
u,v=F(35,64),F(5,256)
x,y=old.AS_point(u,v)
ak=x+y+old.cA*old.APsi_point(u,v)
x,y=unit.AS_point(u,v)
ab=x+y
bad=ak-A(F(3,200))*ab
assert ab>0 and bad<0
assert F(str(bad.upper().fmpq()))<F('-0.0000102147')

bindings=json.loads((HERE/'source-bindings.json').read_text())
assert bindings['source_freeze']['sha256']=='ea1ff01d53e5d2ab3edd9078e802ef05939c7e9c88a057ef51ad9853ab01f36a'
assert bindings['adapted_verifier']['source_sha256']=='e4c91430dde998ce486dc168e0d2440ade51cbb3cf0e3a69f3c650a641f3f6aa'
assert bindings['adapted_verifier']['published_sha256']==sha(HERE/'verify_result.py')
assert bindings['mathematical_result']['claim']==str(claim)
assert bindings['mathematical_result']['direction_halfmass_decrease']=='1/409600'
for item in bindings['unchanged_copies']:
    assert sha(HERE/item['package_path'])==item['source_sha256']
current=provenance()
assert current['hl_check.py']==rp['hl_check.py'] and current['hl_core.py']==rp['hl_core.py']
# Deliberately no equality test of historical/current paths or FLINT binaries.
result={'status':'PASS','scope':'Fresh 192-bit C1/C2/direct mass/negative control; recorded full C4 coverage/provenance/stdout checked, no C4 replay.',
        'package_manifest_sha256':sha(manifest_path),'candidate_sha256':sha(path),'full_report_sha256':sha(report_path),
        'static_package_files_checked':len(manifest['files']),
        'reused_project_files_checked':len(manifest['reused_project_files']),
        'source_copies_checked':len(bindings['unchanged_copies']),
        'recorded_full_c4_precision_bits':128,'fresh_precision_bits':ctx.prec,
        'full_c4_cells':len(rows),'full_c4_leaves':sum(r['leaves'] for r in rows),'full_c4_rerun':False,
        'fresh_c2_lower_rational':str(c2lo),'mass':mass.str(45),
        'mass_lower_rational':str(lo),'mass_upper_rational':str(hi),'claim':str(claim),
        'integral_B_exact':str(int_B),'direction_halfmass_decrease':'1/409600',
        'excluded_ray_epsilon':'[3/200,infinity)','obstruction_upper_rational':str(bad.upper().fmpq()),
        'current_provenance':current,'historical_current_provenance_equal':current==rp}
if args.out is not None:
    args.out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
