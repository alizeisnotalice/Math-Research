"""Three new moment caps only; default hashes, explicit --run prints PID."""
import argparse,hashlib,json,os,sys,time
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx

HERE=Path(__file__).resolve().parent;RUN='audit_I028_fixed_001'
MANIFEST='27f69351a186acfe18cf3071943933350a88de912d482e60a3dd41e60b96c867'
M24='96f65eaa792afc290d3b83aa2284f64170386de045e54f04070b96335f93b764'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def project(p):
    p=p.resolve();return p/'cube-dim-5-bounds' if (p/'cube-dim-5-bounds').is_dir() else p
def ar(x):
    x=F(x);return arb(x.numerator)/x.denominator
def bound(x):return F(str(x.lower().fmpq())),F(str(x.upper().fmpq()))
def encoded(row):return F(row['lower_fraction']),F(row['upper_fraction'])
def parse(s):
    c,r=map(F,s[1:-1].split(' +/- '));return c-r,c+r
def decimal(x,upper=False):
    scale=10**50;n=(x*scale).__ceil__() if upper else (x*scale).__floor__()
    whole,rem=divmod(n,scale);return f'{whole}.{rem:050d}'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',required=True,type=Path);parser.add_argument('--model-root',required=True,type=Path)
    parser.add_argument('--run',action='store_true');args=parser.parse_args()
    source,model=project(args.source_root),project(args.model_root);author=source/'upper_I027'
    manifest=author/'manifest.json';assert sha(manifest)==MANIFEST
    locked={'upper_I027/manifest.json':MANIFEST};info=json.loads(manifest.read_text())
    for row in info['white_list']:
        p=author/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes'];locked['upper_I027/'+row['path']]=sha(p)
    resolved={}
    for row in info['source_inputs']:
        suffix=Path(row['source']).parts[-2:];choices=[r.joinpath(*suffix) for r in (source,model)]
        p=next((p for p in choices if p.exists() and sha(p)==row['sha256']),None);assert p is not None
        resolved[suffix]=p;locked['/'.join(suffix)]=sha(p)
    assert sha(source/'upper_I024/manifest.json')==M24;locked['upper_I024/manifest.json']=M24
    saved=json.loads((author/'results.json').read_text())
    opt=json.loads(resolved[('upper_I024','results.json')].read_text())
    Zlo,Zhi=encoded(opt['Z_exact_interval']);assert 0<Zlo<Zhi and (Zlo,Zhi)==encoded(saved['Zflag_source_interval'])
    old=json.loads(resolved[('upper','round011_shape_integral_results.json')].read_text())
    oldQ=parse(old['arb_256bit']['negative_determinant_moment_upper_Q'])
    if not args.run:
        print(json.dumps({'run_handle':RUN,'status':'locked_inputs_only','calculation_executed':False}));return
    pid=os.getpid();print(json.dumps({'run_handle':RUN,'pid':pid,'status':'started_three_new_points'}),flush=True)
    began=time.monotonic();ctx.prec=512;beta=F(3,4)
    tilt_inverse=16*ar(beta)*ar(beta).gamma();rows=[]
    for gamma,expected_k in ((F(1,5),F(1,4)),(F(3,14),F(9,28)),(F(1,4),F(1,2))):
        k=5*gamma-beta;assert k==expected_k and 0<k<beta
        exponent=1+k/beta-5*gamma;assert exponent==5*gamma/3
        m=ar(gamma/beta).gamma()/(ar(beta)*ar(gamma).gamma())
        positive_unit=ar(-k/beta).gamma()/(ar(beta)*ar(-k).gamma())
        assert m>0 and positive_unit>0
        positive_T=ar(16)**ar(k/beta)*positive_unit
        raw_negative_cap=ar(16)**ar(-5*gamma)*m**5
        prefactor=tilt_inverse*positive_T*raw_negative_cap
        lower=bound(prefactor*ar(Zhi)**(-ar(gamma)))[0]
        upper=bound(prefactor*ar(Zlo)**(-ar(gamma)))[1]
        row=next(row for row in saved['rows'] if F(row['gamma'])==gamma)
        original=encoded(row['cap_interval']);assert original[0]<=lower<upper<=original[1]
        cap=F(row['strict_rational_R_negative_moment_cap']);assert cap==F(info['caps'][str(gamma)]) and upper<cap
        item={'gamma':str(gamma),'k':str(k),'k_over_beta':str(k/beta),'net_16_exponent':str(exponent),
              'new_formula_interval':[str(lower),str(upper)],'decimal_outward':[decimal(lower),decimal(upper,True)],
              'strict_rational_cap':str(cap),'inside_locked_full_formula_interval':True}
        if gamma==F(3,14):
            assert upper<oldQ[0] and cap<oldQ[0]
            assert oldQ==encoded(row['old_locked_Q3_over14_printed_interval'])
            item['old_Q_full_printed_interval']=[str(x) for x in oldQ]
            item['new_cap_below_old_entire_ball']=True
        rows.append(item)
    out={'run_handle':RUN,'pid':pid,'status':'pass_three_true_R_negative_moment_caps',
         'rows':rows,'Harris_domain':'3/20 <= gamma < 3/10; endpoint k=0 is exact, only three interior points computed',
         'Zflag_geometry_or_entropy_recomputed':False,'old_Q_formula_recomputed':False,
         'new_I_bound':False,'new_C5_bound':False,
         'numerical_path':'Uncombined 16 scale factors; Gamma(beta), negative-argument positive Mellin Gamma ratio, unshifted negative Mellin; Arb512.',
         'shared_library_limit':'Gamma/Arb backend is shared; new exact parameters, algebraic scaling and endpoint containment independently checked. No cross-library Gamma certification claimed.',
         'locked_inputs_sha256':locked,'source_sha256':sha(Path(__file__)),
         'registration_sha256':sha(HERE/'registration.md'),
         'runtime':{'python':sys.executable,'Arb_bits':ctx.prec,'seconds':time.monotonic()-began,'new_points':3}}
    (HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('run_handle','pid','status','rows','runtime')},indent=2))

if __name__=='__main__':main()
