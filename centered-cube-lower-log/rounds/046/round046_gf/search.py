#!/usr/bin/env python3
"""Replay the frozen bounded search; verify again using the explicit winning basis."""
import pathlib,json,hashlib,time,subprocess,sys
from certify import field_certificate,call,forms,py_check
HERE=pathlib.Path(__file__).resolve().parent

def main():
    plan=json.loads((HERE/'FROZEN_PLAN.json').read_text())
    assert hashlib.sha256((HERE/'gf_check.c').read_bytes()).hexdigest()==plan['C_source_sha256']
    assert hashlib.sha256((HERE/'PREFLIGHT.json').read_bytes()).hexdigest()==plan['certificate_sha256']
    cert=field_certificate(plan['d'],plan['modulus'],plan['generator'])
    subprocess.run(['/usr/bin/clang','-O3','-Wall','-Wextra','-std=c11',str(HERE/'gf_check.c'),'-o',str(HERE/'gf_check')],check=True)
    start=time.perf_counter()
    runs=call(plan['d'],plan['modulus'],plan['b'],plan['generator'],plan['seed'],plan['max_attempts'])
    (HERE/'SEARCH_ATTEMPTS.json').write_text(json.dumps(runs,indent=2)+'\n')
    winners=[r for r in runs if r['passed']]
    if not winners:
        print('EXHAUSTED',len(runs),'without success');return
    result=winners[0]
    check=call(plan['d'],plan['modulus'],plan['b'],plan['generator'],0,1,result['basis'])[0]
    assert all(result[k]==v for k,v in check.items() if k not in ['attempt','seconds'])
    # Independent Python linear algebra and selected field checks. Full sorted checks are C.
    W=forms(result['basis']);assert len(set(W))==plan['B']
    pairs=[]
    for i in range(1,65):
        a=W[(i*31)%2047+1];b=W[(i*97)%2047+1]
        from certify import mul,power
        inv=power(b,plan['K']-1,plan['modulus']);assert mul(b,inv,plan['modulus'])==1
        pairs.append({'a':a,'b':b,'product':mul(a,b,plan['modulus']),'inverse_b':inv,'ratio':mul(a,inv,plan['modulus'])})
    subprocess.run(['/usr/bin/clang','-O3','-std=c11',str(HERE/'gf_spot.c'),'-o',str(HERE/'gf_spot')],check=True)
    bridge=subprocess.run([str(HERE/'gf_spot'),str(plan['d']),str(plan['modulus']),str(plan['generator'])],input=''.join(f"{r['a']} {r['b']}\n" for r in pairs),text=True,capture_output=True,check=True)
    for r,line in zip(pairs,bridge.stdout.splitlines(),strict=True):
        v=list(map(int,line.split()));assert v==[r['product'],r['inverse_b'],r['ratio'],mul(plan['generator'],r['a'],plan['modulus'])]
        r['C_match']=True
    manifest={'schema':'round046.exact_gf_source.v1','d':plan['d'],'b':plan['b'],'B':plan['B'],'q':plan['q'],'K':plan['K'],'modulus':plan['modulus'],'modulus_hex':hex(plan['modulus']),'generator':plan['generator'],'basis':result['basis'],'seed':plan['seed'],'max_attempts':plan['max_attempts'],'successful_attempt':result['attempt'],'attempts_used':len(runs),'parameters':plan['accepted_target_parameters'],'field_certificate':cert,'validation':result,'explicit_basis_recheck':check,'independent_Python_span_size':len(set(W)),'independent_Python_selected_arithmetic_checks':len(pairs),'construction_formula':{'x_j':'g^j for 0<=j<K','y_i(j)':'Tr(beta_i*x_j)','y_j':'sum_i 2^i*y_i(j)','U_j':'(y_j+1/2)/B','a_j':'(j+U_j)/K'},'search_and_recheck_seconds':time.perf_counter()-start,'claim_scope':'Exact GF source hypotheses only. This does not certify a new log factor for C_n.','source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['gf_check.c','gf_spot.c','certify.py','search.py','FROZEN_PLAN.json']}}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (HERE/'PYTHON_SPOT_CHECKS.json').write_text(json.dumps(pairs,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
