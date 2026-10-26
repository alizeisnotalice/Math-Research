#!/usr/bin/env python3
"""Register schema/source and one fresh symbolic diagnostic;no real factor."""
import ast,gzip,hashlib,json,os
from pathlib import Path
from fractions import Fraction as F

HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p);return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def write(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')

def main():
    assert not (HERE/'registration.json').exists()
    names=['upper_time159/proof.md','upper_time159/symbolic159.py','upper_time159/symbolic_results.json',
      'upper_alpha1_029/report.md','upper_mean085/report.md','upper_shift091/report.md',
      'upper_pathmean111/report.md','upper_pathmean111/pathmean.py','upper_pathmean111/secant.py','upper_pathmean111/geometry_plans.json',
      'upper_path112/geometries.json.gz','upper_path112/pathmean.py','upper_path112/factor_catalog.json',
      'upper_integral150/report.md','upper_integral150/consumer.py','upper_integral150/preflight_results.json','upper_integral150/actual/factors.json']
    amendment=MAIN/'upper_time159/amendment159.md'
    if amendment.exists():names.append(str(amendment.relative_to(MAIN)))
    sources={str(MAIN/n):sha(MAIN/n) for n in names}
    pre=read(MAIN/'upper_integral150/preflight_results.json');times=pre['time_endpoints'];assert len(times)==73 and all(F(a)<F(b) for a,b in zip(times,times[1:])) and F(times[0])>0
    geos=read(MAIN/'upper_path112/geometries.json.gz');catalog=read(MAIN/'upper_path112/factor_catalog.json');rows=read(MAIN/'upper_integral150/actual/factors.json')
    assert len(geos)==5 and set(geos)==set(catalog) and len(rows)==73
    count=0
    for gid,g in geos.items():
        p=g['profile'];assert set(['h','mu_infinity','J0','M2','offsets','pieces'])<=set(p) and len(p['offsets'])==4
        j0,m2=F(p['J0']),F(p['M2']);assert j0>=0 and m2>=0 and F(p['mu_infinity'])>=0
        # Relational stored-field validation,not a profile/J0 recomputation.
        assert 2*j0.numerator*m2.denominator==m2.numerator*j0.denominator
        assert set(catalog[gid])==set(times)
        for t,row in catalog[gid].items():
            assert set(['J_upper','factor_upper','plan_sha256'])<=set(row)
            assert (1<<192)%F(row['factor_upper']).denominator==0;count+=1
    for i,row in enumerate(rows):
        assert row['a']==times[i] and row['b']==(times[i+1] if i+1<73 else None)
        assert (1<<192)%F(row['new_upper']).denominator==0;count+=1
    old=HERE.parent/'audit_final147';assert sha(old/'manifest.json')=='141eaf7bf1883d5fabeaa22b76eb07ede520648de5cdab56504cbd98f98e6ebb'
    protected={str(old/'manifest.json'):sha(old/'manifest.json'),**{str(old/n):h for n,h in read(old/'manifest.json')['files'].items()}}
    for p,h in protected.items():assert sha(p)==h
    for p in HERE.glob('*.py'):ast.parse(p.read_text())
    code={str(p):sha(p) for p in HERE.glob('*.py')}
    write(HERE/'schema_receipt.json',dict(status='real111_112_150_stored_schema_and_grid_linked_no_new_factor',pid=os.getpid(),profile_geometries=list(sorted(geos)),profile_count=5,stored_J_endpoint_rows=365,stored150_factor_rows=73,stored192_grid_fields_checked=count,
        source_time_key_semantics='positive physical time t;J_upper encloses J(1/t)',zero_time_key_allowed=False,infinite_tail_key=None,
        profile_semantics=dict(mu_infinity='nonnegative infinity norm of mean offset',J0='stored exact M2/2',M2='stored squared sup-offset second moment'),
        real_new_scalar_factor_mean_corner_H_R_J_calls=0,profile_evaluator_calls=0,old_baseline_replays=0,source153_imports=0))
    write(HERE/'registration.json',dict(task='audit_time160_pure_proof_interface',case='capped_endpoint_upper_polynomial_first_second_derivatives_only',source_locks=sources,code_locks=code,protected147_locks=protected,
        new_symbolic_diagnostic_authorized=True,old159_five_groups_must_not_replay=True,real_scalar_factor_mean_corner_H_R_J_authorized=False,finite_symbolic_cases=1,workers=0,wall_cap=None,
        future_factor_contract='pure moment_core functions only;explicit source/geometry/time/catalog/oldfactor binding and separate primary authorization required before real calls',
        future_time_coefficients_at_most=73,real_factor_cases_registered=0,real_mean_cases_registered=0))
    print(json.dumps(dict(registration_sha256=sha(HERE/'registration.json'),source_files=len(sources),code_files=len(code),schema_grid_fields=count,new_real_factor_mean_calls=0)),flush=True)

if __name__=='__main__':main()
