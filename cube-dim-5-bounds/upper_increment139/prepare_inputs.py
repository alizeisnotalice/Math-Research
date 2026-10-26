#!/usr/bin/env python3
"""Copy the one qualified fixed input by byte/JSON bindings only. No J math."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,os

HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_bytes())
def write(name,x):
    with (HERE/name).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def main():
    assert not (HERE/'inputs.json').exists(),'Preserve fixed input'
    original=MAIN/'round015/join135_inputs.json';prior=MAIN/'round015/join135_actual/results.json'
    assert sha(original)=='c6eb027efee798f140f6cb7362721b246de257a722b58cf65ba57ea1a6ad426c'
    assert sha(prior)=='8b0623db708b493cf7d196e5f8dc0f88724d401154a7c888da95987568502686'
    old=read(original);previous=read(prior);accepted_path=MAIN/'round015/primary_upper_acceptance.json';accepted=read(accepted_path)
    assert accepted['status']=='primary_accepts_complete_C5_upper' and accepted['g']==old['g']=='21/200000'
    locks={'original_fixed_input':{'path':str(original),'sha256':sha(original)},'saved135_result':{'path':str(prior),'sha256':sha(prior)},
        'accepted_round015':{'path':str(accepted_path),'sha256':sha(accepted_path)}}
    for quantity,key in [('Cbar','Cbar_upper'),('M','M_upper'),('IH','IH_upper')]:
        row=old['qualifications'][quantity];assert sha(row['source_path'])==row['source_sha256']
        value=read(row['source_path'])
        for part in row['field'].split('.'):value=value[part]
        assert value==old[key]==previous[key]
        locks[quantity]={'path':row['source_path'],'sha256':row['source_sha256']}
    pi_path=MAIN/'upper_shorttime076/results.json';assert sha(pi_path)=='caeda104af7973ccb153daf44cfd3711cac332ae262f5960415263a32f102ec3'
    pi=read(pi_path)['derived_constants']['pi_outward_384_dyadic'];assert previous['pi_interval']==[pi['lower'],pi['upper']]
    locks['stored_pi']={'path':str(pi_path),'sha256':sha(pi_path)}
    for name,path,expected in [('proof133',MAIN/'upper_rearrangement133/report.md','7e811dbecfc823aad8553c9bdf659242e2e09ca5f92f7b1fb0d2ff889bce4bde'),
        ('transfer029',MAIN/'upper_alpha1_029/report.md','298f0de7277ea6a7a516b85344c5189562bf7026f3f8de43d4d3f304ff6fed81')]:
        assert sha(path)==expected;locks[name]={'path':str(path),'sha256':expected}
    out={k:old[k] for k in ['epsilon','g','Cbar_upper','M_upper','IH_upper','qualifications']}
    out.update({'pi_interval':[pi['lower'],pi['upper']],'saved135_cap_upper':previous['cap_upper'],
        'accepted_I_upper':accepted['complete_I_upper'],'source_locks':locks,
        'auxiliary_expression':'5epsilon^2/(2pi)*integral_0^1 z^2[D+(M-D)z^4]^(1/4)dz;D=Cbar/(2epsilon^2)',
        'actual_J_R_I_lower_claim':False})
    write('inputs.json',out)
    write('source_binding.json',{'task':'upper_increment139_fixed_input_source_binding','pid':os.getpid(),
        'utc':datetime.now(timezone.utc).isoformat(),'prepare_code_sha256':sha(__file__),'inputs_sha256':sha(HERE/'inputs.json'),
        'source_locks':locks,'status':'PASS_exact_source_fields_copied_no_arithmetic','actual_J_cases':0,
        'actual_root_atoms':0,'actual_panels':0,'old_scalar_norm_Gamma_or_toy_runs':0})
    print(json.dumps({'status':'fixed_input_prepared_only','pid':os.getpid(),'inputs_sha256':sha(HERE/'inputs.json'),'actual_J_cases':0},sort_keys=True))
if __name__=='__main__':main()
