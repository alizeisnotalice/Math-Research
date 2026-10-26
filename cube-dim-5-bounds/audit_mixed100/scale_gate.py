#!/usr/bin/env python3
"""Exact finite binding of EPS128 values to EPS64 evaluation at doubled points."""
import sys
sys.dont_write_bytecode=True
import json,gzip,os,time
from pathlib import Path
from fractions import Fraction as F
from math import lcm
from interfaces import REG,HERE,common as C


def main():
    start=time.monotonic();rp=HERE/'scale_gate_registration.json';registration=json.loads(rp.read_text())
    for p,h in registration['input_locks'].items():assert C.file_sha(p)==h
    fan=json.loads(gzip.decompress(Path(REG['fan']).read_bytes()))
    tasks=json.loads(Path(registration['execution_tasks']).read_text());assert len(tasks)==194
    dilation=F(2);dimension=5;time_jacobian=dilation;density_multiplier=dilation**(-dimension)
    H_multiplier=1/(density_multiplier*time_jacobian);gradient_multiplier=H_multiplier*dilation
    assert H_multiplier==16 and gradient_multiplier==32
    assert F(1,64)/dilation==F(1,128)
    assert 1/(16*F(1,64)**2)==256 and 1/(16*F(1,128)**2)==1024
    proofs=[];seen=set()
    for row in tasks:
        point=list(map(F,row['point']));doubled=[dilation*x for x in point];D=lcm(*(x.denominator for x in point));E=lcm(*(x.denominator for x in doubled))
        integers=[int(x*D) for x in point];scaled=[int(x*E) for x in doubled];ratio=2*F(E,D)
        assert point[0]==1 and doubled[0]==2 and all(F(b)==ratio*a for a,b in zip(integers,scaled))
        original_projections=[sum(a*b for a,b in zip(v,integers)) for v in fan['vertices']]
        doubled_projections=[sum(a*b for a,b in zip(v,scaled)) for v in fan['vertices']]
        assert len(set(original_projections))==len(set(doubled_projections))==280 and 0 not in original_projections
        assert all(b==ratio*a for a,b in zip(original_projections,doubled_projections))
        assert F(144*E,64)==ratio*F(144*D,128)
        # Physical pole ratio and time parameter scale together. No H/log evaluated.
        pair=('1/128',','.join(map(str,point)));new_input=','.join(map(str,doubled));assert new_input not in seen;seen.add(new_input)
        assert row['epsilon']==pair[0] and row['pair_cache_key']==list(pair)
        proofs.append({'ordinal':row['ordinal'],'source_point':list(map(str,point)),'source_pair_key':list(pair),'independent_EPS64_input':list(map(str,doubled)),
                       'H_multiplier':16,'four_physical_gradient_multiplier':32,'projection_time_common_scale':str(ratio),'source_D':str(D),'independent_D':str(E)})
    for p,h in registration['input_locks'].items():assert C.file_sha(p)==h
    result={'task':'audit_mixed100_scale_gate','status':'EPS128_to_EPS64_doubled_input_exact_binding_qualified_no_values','pid':os.getpid(),'wall_seconds':time.monotonic()-start,
            'registration_sha256':C.file_sha(rp),'p_t_identity':'p_t(2x)=2^-5 p_(t/2)(x)','time_substitution':'u=t/2, dt=2du; lower1/64 becomes1/128',
            'H_identity':'H_(1/128)(x)=16H_(1/64)(2x)','gradient_identity':'partial_xj H_(1/128)(x)=32 partial_yj H_(1/64)(y) at y=2x, j=2..5',
            'same_epsilon_homogeneity_used':False,'Q128_Taylor_penalty':1024,'new_H_gradient_log_calls':0,'finite_inputs':194,'inputs':proofs}
    path=HERE/'scale_gate_result.json';path.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['pid','status','wall_seconds','finite_inputs','new_H_gradient_log_calls']}));print(json.dumps({'result_sha256':C.file_sha(path)}))


if __name__=='__main__':main()
