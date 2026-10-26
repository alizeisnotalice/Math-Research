"""Three new exact logical/linear-algebra checks; no body/source integration."""
import signal
import time
START = time.monotonic()
CPU_START = time.process_time()
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.alarm(170)
import resource
hard = resource.getrlimit(resource.RLIMIT_CPU)[1]
CPU_LIMIT = min(25, hard) if hard != resource.RLIM_INFINITY else 25
resource.setrlimit(resource.RLIMIT_CPU, (CPU_LIMIT, CPU_LIMIT))
import hashlib
import json
import os
from pathlib import Path
import sys
import traceback
sys.dont_write_bytecode = True
P = Path(__file__).resolve().parent
LOCK = json.loads((P/'input_manifest.json').read_text())
SPEC = json.loads((P/'fixture_spec.json').read_text())
ROWS = []


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def save(path, value):
    with path.open('x') as stream:
        json.dump(encode(value), stream, indent=2, sort_keys=True)
        stream.write('\n')


def pins():
    assert sha(P/'registration.json') == LOCK['registration_sha256']
    author = Path(LOCK['author_root'])
    assert sha(author/'manifest.json') == LOCK['author_manifest_sha256']
    for name, digest in LOCK['author_artifact_sha256'].items():
        assert sha(author/name) == digest, name
    for name, digest in json.loads((P/'inherited_input_pins.json').read_text())['readonly_sha256'].items():
        assert sha(Path(name)) == digest, name
    for name, digest in json.loads((P/'pre_run_manifest.json').read_text())['sha256'].items():
        assert sha(P/name) == digest, name


def main():
    pins()
    save(P/'run_started.json', dict(pid=os.getpid(), python=sys.executable,
        pre_run_manifest_sha256=sha(P/'pre_run_manifest.json'), CPU_hard_limit_seconds=CPU_LIMIT,
        total_alarm_seconds=170, automatic_retries=0))
    code, failure, current = 0, None, None
    try:
        import sympy as s

        current = 'C1_closed_set_and_singular_source'
        data = SPEC['configurations'][0]
        n = s.Symbol('n', positive=True, integer=True)
        continuous_test = (1-1/n)/3
        limit_test = s.limit(continuous_test, n, s.oo)
        closed_restricted_limit = s.Rational(0)
        limiting_closed_mass = s.Rational(1, 3)
        assert limit_test == limiting_closed_mass
        assert closed_restricted_limit < limiting_closed_mass
        mu = list(map(s.Rational, data['three_component_measure']))
        a = list(map(s.Rational, data['signed_correction']))
        eps = s.Rational(data['epsilon'])
        invalid = s.Rational(data['invalid_epsilon'])
        changed = [x-eps*y for x,y in zip(mu,a)]
        invalid_changed = [x-invalid*y for x,y in zip(mu,a)]
        assert sum(a)==0 and sum(mu)==sum(changed)==1
        assert min(changed)>=0 and changed[0]==mu[0]
        assert all(max(y,0)<=3*x for x,y in zip(mu,a))
        assert min(invalid_changed)<0
        ROWS.append(dict(id=current,status='PASS',closed_set_limsup=closed_restricted_limit,
            limiting_closed_mass=limiting_closed_mass,reverse_Portmanteau_is_false=True,
            valid_corrected_masses=changed,invalid_amplitude_masses=invalid_changed,
            singular_component_unchanged=True,source_integrals=0))
        print(current,'PASS',flush=True)

        current = 'C2_cylinder_transversality'
        data = SPEC['configurations'][1]
        v,w,r,z = [s.Matrix(list(map(s.Rational,data[key]))) for key in ('v','w','r','z')]
        param = s.Symbol('s', real=True)
        def projection(columns):
            b=s.Matrix.hstack(*columns)
            return s.eye(6)-b*(b.T*b).inv()*b.T
        containing=projection((v,w))
        assert containing*v==s.zeros(6,1)
        assert s.expand(((z+param*v).T*containing*(z+param*v))[0]-(z.T*containing*z)[0])==0
        transverse=projection((w,r)); enlarged=projection((v,w,r))
        det2=(s.Matrix.hstack(w,r).T*s.Matrix.hstack(w,r)).det()
        det3=(s.Matrix.hstack(v,w,r).T*s.Matrix.hstack(v,w,r)).det()
        slope=(v.T*transverse*v)[0]
        minimum=(z.T*enlarged*z)[0]
        shift=(z.T*transverse*v)[0]/slope
        distance_squared=((z+param*v).T*transverse*(z+param*v))[0]
        assert slope==det3/det2 and slope>=s.Rational(1,36)
        assert minimum>0 and s.expand(distance_squared-minimum-slope*(param+shift)**2)==0
        # Off-S2 alone is insufficient: three distinct absolute coordinate
        # values rule out every body pair plane, but this center lies in S3.
        bad=2*v+w+2*r
        assert len(set(abs(x) for x in bad))==3
        assert transverse*(bad-2*v)==s.zeros(6,1)
        ROWS.append(dict(id=current,status='PASS',pair_Gram_det=det2,triple_Gram_det=det3,
            transverse_direction_distance_squared=slope,full_line_minimum_squared=minimum,
            minimizing_parameter=-shift,bad_off_S2_center=list(bad),
            bad_cylinder_hits_pair_plane_at=-2,body_frames_enumerated=0,source_integrals=0))
        print(current,'PASS',flush=True)

        current = 'C3_general_source_transfer_model'
        data = SPEC['configurations'][2]; size=data['group_size']
        vector=lambda key:s.Matrix(list(map(s.Rational,data[key])))
        k0,mu,b,bad,u,sigma=[vector(key) for key in ('K0','mu0','good_b','bad_b','u','sigma')]
        eps=s.Rational(data['epsilon']);kappa=s.Rational(data['kappa'])
        origin=s.Matrix([1]+[0]*(size-1))
        def generator(x):
            return s.Matrix([x[(i-1)%size]+x[(i+1)%size]-2*x[i] for i in range(size)])
        def convolution(x,y):
            return s.Matrix([sum(x[(i-j)%size]*y[j] for j in range(size)) for i in range(size)])
        assert generator(k0)==mu-origin and min(k0)>0 and sum(mu)==1
        correction=generator(b);k=k0-eps*b;changed_mu=mu-eps*correction
        assert min(k)>=0 and min(changed_mu)>=0 and sum(changed_mu)==1
        assert correction[0]==0 and generator(k)==changed_mu-origin
        assert sum(k)==sum(k0)-eps and sum(correction)==0
        assert all(max(correction[i],0)<=2*mu[i] for i in range(size))
        bad_kernel=k0-eps*bad;bad_source=mu-eps*generator(bad)
        assert min(bad_kernel)>0 and bad_source[0]<0
        forcing=sigma-generator(u)
        assert min(forcing)>=0 and min(sigma)>=0 and max(sigma)<=kappa
        active=[i for i in range(size) if u[i]>0]
        assert all(sigma[i]==kappa for i in active) and sum(forcing)==sum(sigma)
        assert convolution(k,generator(u))==convolution(changed_mu,u)-u
        outside=[i for i in range(size) if i not in active]
        conv=convolution(k,forcing)
        assert all(conv[i]<=kappa*sum(k) for i in outside)
        ROWS.append(dict(id=current,status='PASS',corrected_kernel=list(k),corrected_probability=list(changed_mu),
            mass_decrement=eps,invalid_bump_still_positive_kernel=list(bad_kernel),
            invalid_bump_negative_source=list(bad_source),obstacle_set=active,
            forcing=list(forcing),outside_convolution_values=[conv[i] for i in outside],
            comparison_upper=kappa*sum(k),continuum_all_scale_claim_from_toy=False,source_integrals=0))
        print(current,'PASS',flush=True)
        pins()
    except BaseException as exc:
        code=1;failure=dict(configuration=current,error=repr(exc),traceback=traceback.format_exc())
    result=dict(schema='A095-independent-new-check-results-v1',
        status='PASS_THREE_NEW_EXACT_CONFIGURATIONS' if code==0 else 'FAIL_STOP_NO_RETRY',
        actual_configurations=len(ROWS)+(1 if failure else 0),passed_configurations=len(ROWS),
        results=ROWS,failure=failure,wall_seconds=time.monotonic()-START,CPU_seconds=time.process_time()-CPU_START,
        body_kernel_calls=0,source_integrals=0,old_sanity_reruns=0,automatic_retries=0,
        new_C6_bound=False,continuum_proof_replaced_by_tests=False)
    save(P/'results.json',result);print(json.dumps(encode(result),sort_keys=True),flush=True)
    return code


if __name__=='__main__':sys.exit(main())
