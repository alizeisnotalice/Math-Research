"""Seven new exact rational implications, pending complete analytic acceptance.

No kernel, numerical log, optimization, old baseline, or floating arithmetic.
Also emit a function-local MATLAB Symbolic cross-check of the same finite chain.
"""
from fractions import Fraction as F
import argparse
import hashlib
import json
import os
from pathlib import Path
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--matlab-code', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.matlab_code.exists()
    started = time.monotonic()
    reg = json.loads(args.registration.read_bytes())
    assert sha(__file__) == reg['code_sha256'] and reg['cases'] == 7
    for item in reg['sources'].values():
        assert sha(item['path']) == item['sha256']
    def source(key):
        return json.loads(Path(reg['sources'][key]['path']).read_bytes())
    c = source('facepeak134')
    m = source('sup132')
    old = source('primary_I112')
    p = source('stored076')
    C = F(c['Cbar_upper_256dyadic'])
    M = F(m['objects']['analytic_stationary']['M_upper_portable_cap'])
    H = F(old['complete_H512_integral_upper'])
    pi_lower = F(p['derived_constants']['pi_outward_384_dyadic']['lower'])
    eps, g = F(1, 512), F(21, 200000)
    Ccap, Mcap, Hcap = F(1, 300000), F(2), F(44372, 10**9)
    pcap, rcap, Icap, Ucap = F(157, 50), F(33, 32), F(2271, 50000000), F(799, 25)
    D = Ccap / (2*eps**2)
    root_argument = (4*D+3*Mcap)/7
    Jcap = F(5, 6)*eps**2*rcap/pcap
    margins = [Ccap-C, Mcap-M, Hcap-H, pi_lower-pcap,
               rcap**4-root_argument, Icap-(Hcap+Jcap),
               Ucap**4*g**5-3125*Icap**4]
    labels = ['Cbar<=1/300000', 'M<=2', 'IH<=44372/1e9',
              'pi>=157/50', 'fourthroot<=33/32',
              'complete_I<=2271/50000000', 'conditional_C5<=799/25']
    values = [str(v) for v in (C, M, H, pi_lower)]
    expr = (
        "feval(@(C,M,H,P,e,g) disp(["
        "sign(sym(1)/300000-C),sign(sym(2)-M),"
        "sign(sym(44372)/sym(10)^9-H),sign(P-sym(157)/50),"
        "sign((sym(33)/32)^4-(4*((sym(1)/300000)/(2*e^2))+6)/7),"
        "sign(sym(2271)/50000000-(sym(44372)/sym(10)^9+"
        "sym(5)/6*e^2*(sym(33)/32)/(sym(157)/50))),"
        "sign((sym(799)/25)^4*g^5-3125*(sym(2271)/50000000)^4)]),"
        + ','.join('str2sym("'+v+'")' for v in values)
        + ',sym(1)/512,sym(21)/200000);\n'
    )
    args.matlab_code.write_text(expr)
    result = {'status': 'seven_rational_implications_complete_analytic_qualification_separate',
              'pid': os.getpid(), 'seconds': time.monotonic()-started,
              'registration_sha256': sha(args.registration), 'cases': 7,
              'checks': [{'claim': name, 'strict_margin': str(value), 'passed': value > 0}
                         for name, value in zip(labels, margins)],
              'all_passed': all(value > 0 for value in margins),
              'input_Cbar': str(C), 'input_M': str(M), 'input_IH': str(H),
              'input_pi_lower': str(pi_lower), 'epsilon': str(eps), 'qualified_g': str(g),
              'simple_root_argument': str(root_argument), 'simple_J_upper': str(Jcap),
              'simple_I_cap': str(Icap), 'simple_C_cap': str(Ucap),
              'conditional_gain_against32': str(32-Ucap),
              'matlab_code_sha256': sha(args.matlab_code),
              'formal_C5_claim': False,
              'needs': 'complete134/133/135 theory, same-model IH/M/g105 and029, independent136 and primary',
              'new_H_gradient_log_or_old_baseline_calls': 0, 'Astra_calls': 0}
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'pid', 'cases', 'all_passed', 'simple_C_cap')}))


if __name__ == '__main__':
    main()
