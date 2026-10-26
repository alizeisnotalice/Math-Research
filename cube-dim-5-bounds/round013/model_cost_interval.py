"""New cost-interval/S accounting from frozen complete certificates only."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib
import json
import os
import time

H = Path(__file__).resolve().parent
P = H.parent


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def fourth(q):
    D = 1 << 256
    a = isqrt(isqrt(q.numerator*D**4//q.denominator))
    assert a**4*q.denominator <= q.numerator*D**4 < (a+1)**4*q.denominator
    return F(a,D), F(a+1,D)


def main():
    start = time.monotonic()
    rp = H/'model_cost_interval_registration.json'
    reg = json.loads(rp.read_text())
    for name, expected in reg['locks'].items():
        assert sha(P/name) == expected, name
    acceptance = json.loads((H/'primary_I095_acceptance.json').read_text())
    gu_source = json.loads((P/'audit_zero041/result.json').read_text())
    assert gu_source['status'] == 'PASS_selected019_first_point_time_resolvent_fullGreen'
    assert gu_source['point'] == ['1','1','1','2/3','1/16']
    # This point has three distinct absolute coordinates. A combination of two
    # body directions has only |a+b| and |a-b|, so the point is outside S2.
    assert len(set(map(lambda x: abs(F(x)),gu_source['point']))) == 3
    gL, gU = F(acceptance['qualified_g']), F(gu_source['pure_Fraction_G_interval'][1])
    b = json.loads((P/'upper_adaptive095/run_receipt.json').read_text())['after']
    IL, IU = F(b['LH_truncated_integral']), F(acceptance['I_upper'])
    UH = F(b['UH_truncated_integral'])
    assert 0 < gL < gU and 0 < IL <= UH <= IU
    low4, up4 = 3125*IL**4/gU**5, 3125*IU**4/gL**5
    costlo, costhi = fourth(low4)[0], fourth(up4)[1]
    threshold = tuple(32*x for x in fourth((gL/5)**5))
    S_certificate = [threshold[0]-IU, threshold[1]-IU]
    S_for_fixed_gL_true_I = [threshold[0]-IU, threshold[1]-IL]
    assert S_certificate[1] < 0 and costlo < 32 < costhi
    exact_gate = 32**4*gL**5-3125*IU**4
    assert exact_gate < 0
    out = {'status':'complete_fixed_alpha1_model_cost_interval_crosses32',
           'model':'n5 alpha1 body16 unnormalized weight1 common time Fourier(2pi)^-5',
           'g_L':str(gL),'g_U':str(gU),'I_L':str(IL),'I_U':str(IU),
           'g_U_qualification':'041 complete rule point x=(1,1,1,2/3,1/16) lies outside S2 (three distinct absolute coordinates);029 compact small-time Ct and long-time Ct^-5 imply continuity, hence ess inf cube <= point value',
           'I_L_qualification':'G>=H512>=0, so full I >= certified full-face integral of H512^q >= LH; not a pointwise lower bound for G',
           'model_cost_lower_fourth_power':str(low4),
           'model_cost_upper_fourth_power':str(up4),
           'model_cost_outward_interval':[str(costlo),str(costhi)],
           'model_cost_display_only':[float(costlo),float(costhi)],
           'fixed_model_excluded':False,'not_a_C5_lower_bound':True,
           'I_threshold_for_qualified_g105_interval':list(map(str,threshold)),
           'certified_budget_S_interval':list(map(str,S_certificate)),
           'fixed_g105_true_I_S_interval':list(map(str,S_for_fixed_gL_true_I)),
           'exact_certificate_32_gate_fourth_power_margin':str(exact_gate),
           'certificate_S_display_only':list(map(float,S_certificate)),
           'combined_spatial_time_envelope_gap_UH_minus_LH':str(UH-IL),
           'shorttime_norm_payment_plus_final_cap_rounding_IU_minus_UH':str(IU-UH),
           'gap_to_threshold_outward_interval':[str(IU-threshold[1]),str(IU-threshold[0])],
           'loss_limits':'UH-LH includes positive-time spatial/range/time-envelope losses; the old data do not yet separately certify each true loss. IU-UH is an upper-budget difference, not an equality for the actual short-time error.',
           'selection':'prioritize new H512 mean-variance zero-heat consumer before further bulk spatial refinement; g110 prepared repair remains unqualified',
           'new_H_or_old_baseline_calls':0,'pid':os.getpid(),
           'seconds':time.monotonic()-start,'source_sha256':sha(__file__),
           'registration_sha256':sha(rp),'input_locks':reg['locks']}
    target = H/'model_cost_interval.json'
    assert not target.exists()
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','model_cost_display_only',
         'certificate_S_display_only','fixed_model_excluded','pid','seconds')}))


if __name__ == '__main__':
    main()
