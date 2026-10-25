"""Selected018 only: barycentric exact log coefficients and rational atanh bounds."""
import argparse,gzip,hashlib,json,math,os,sys,time
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent;RUN='audit_green036_fixed_001';N=96;GRID=2**384
MANIFEST='e387fd41e3bebc8ba94427200fb13fc6693f15c92118796b34d704d738d018bc'
TABLE='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
GATE={'alpha1_gate.py':'a6deb402cdd9ea8d1e212486364d7cc2a4f1573c8262afb1df576e1cdacc7293',
      'alpha1_gate.json':'2e8601968b6ed99582ab9885efda11dcd892e1ffb70ae0995350a868aaaf3142',
      'alpha1_gate_proof.md':'e348a5480c6570da3a9f9481b567e017ca481213fd97e5b393823d654f7df04b'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def project(p):
    p=p.resolve();return p/'cube-dim-5-bounds' if (p/'cube-dim-5-bounds').is_dir() else p
def parse(s):
    a,b=map(F,s[1:-1].split(' +/- '));return a-b,a+b
def log_bound(u):
    assert 0<=u<=F(1,3)
    if not u:return 0,0
    power=u;s=F(0);u2=u*u
    for j in range(N):s+=2*power/(2*j+1);power*=u2
    hi=s+2*power/F(2*N+1)/(1-u2)
    return (s*GRID).__floor__(),(hi*GRID).__ceil__()
def pi_bounds():
    def atan(q):
        lo=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(80)),F(0))
        return lo,lo+F(1,161*q**161)
    a,b=atan(5),atan(239);return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',required=True,type=Path);parser.add_argument('--model-root',required=True,type=Path)
    parser.add_argument('--run',action='store_true');args=parser.parse_args()
    source,model=project(args.source_root),project(args.model_root);author=source/'astra_green018'
    manifest=author/'manifest.json';assert sha(manifest)==MANIFEST;locked={'astra_green018/manifest.json':MANIFEST}
    for row in json.loads(manifest.read_text())['white_list']:
        p=author/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes'];locked['astra_green018/'+row['path']]=sha(p)
    tablepath=source/'astra_vertex017/boundary_triangulation.json.gz';assert sha(tablepath)==TABLE;locked['astra_vertex017/boundary_triangulation.json.gz']=TABLE
    for name,expected in GATE.items():
        p=model/'round009'/name;assert sha(p)==expected;locked['round009/'+name]=expected
    gate=json.loads((model/'round009/alpha1_gate.json').read_text());resolved={}
    for legacy,expected in gate['source_sha256'].items():
        suffix=Path(legacy).parts[-2:];options=[r.joinpath(*suffix) for r in (source,model)]
        p=next((p for p in options if p.exists() and sha(p)==expected),None);assert p is not None
        locked['/'.join(suffix)]=sha(p);resolved[suffix]=p
    batch=json.loads((author/'batch_result.json').read_text());selected=batch['selected_value']
    assert batch['status']=='complete' and selected==batch['records'][batch['best_record']] and batch['best_record']==49
    x=tuple(map(F,selected['point']));assert x==(F(1),F(104115,104729),F(103354,104729),F(71094,104729),F(6327,104729))
    if not args.run:
        print(json.dumps({'run_handle':RUN,'status':'locked_inputs_only','calculation_executed':False}));return
    pid=os.getpid();began=time.monotonic();print(json.dumps({'run_handle':RUN,'pid':pid,'status':'started_selected_new_point_only'}),flush=True)
    table=json.loads(gzip.decompress(tablepath.read_bytes()));Z=table['vertices'];simplices=table['simplices']
    assert table['denominator']==144 and len(Z)==280 and len(simplices)==3504
    L=math.lcm(*(v.denominator for v in x));n=[int(v*L) for v in x]
    k=[sum(a*b for a,b in zip(n,z)) for z in Z]
    assert all(k) and len(set(k))==280 and set(k)=={-v for v in k}
    A=[F(0)]*280
    for s in simplices:
        ids=s['vertices'];D=s['det_numerator'];local=[]
        for i in ids:
            value=F(D,math.prod(k[i]-k[j] for j in ids if j!=i));local.append(value);A[i]+=value
        assert sum(local,F(0))==0
    assert sum(A,F(0))==0
    raw={}
    for kv,av in zip(k,A):raw[abs(kv)]=raw.get(abs(kv),F(0))+av
    prefactor=-F(L**4,144*32)
    C={key:prefactor*value for key,value in raw.items() if value}
    assert len(raw)==140 and len(C)<=140 and sum(C.values(),F(0))==0
    reference=max(C);terms={};eweight=F(0)
    for integer,c in C.items():
        if integer==reference:continue
        y=F(integer,reference);e=0
        while y<1:y*=2;e-=1
        while y>=2:y/=2;e+=1
        u=(y-1)/(y+1);assert 0<=u<F(1,3)
        if u:terms[u]=terms.get(u,F(0))+c
        eweight+=e*c
    terms[F(1,3)]=terms.get(F(1,3),F(0))+eweight
    terms={u:c for u,c in terms.items() if c};assert len(terms)<=140
    coefficients={'point':[str(v) for v in x],'L':L,'n':n,'integer_projections':k,
                  'raw_vertex_A_coefficients':[str(a) for a in A],
                  'signed_prefactor_for_G_times_pi5':str(prefactor),
                  'merged_abs_C_coefficients':{str(key):str(value) for key,value in sorted(C.items())},
                  'sum_A':'0','sum_C':'0','log_reference':reference,
                  'rescaled_atanh_terms':{str(u):str(c) for u,c in sorted(terms.items())}}
    coefficientpath=HERE/'coefficients.json';coefficientpath.write_text(json.dumps(coefficients,indent=2)+'\n')
    common=math.lcm(*(c.denominator for c in terms.values()));lo_num=hi_num=0
    for i,(u,c) in enumerate(sorted(terms.items())):
        low,high=log_bound(u);weight=c.numerator*(common//c.denominator)
        lo_num+=weight*(low if weight>=0 else high);hi_num+=weight*(high if weight>=0 else low)
        if time.monotonic()-began>120:raise TimeoutError('finite log-batch checkpoint; exact coefficients saved')
    Slo,Shi=F(lo_num,common*GRID),F(hi_num,common*GRID);assert 0<Slo<Shi
    pil,pih=pi_bounds();glo,ghi=Slo/pih**5,Shi/pil**5
    scale=10**85;glo,ghi=F((glo*scale).__floor__(),scale),F((ghi*scale).__ceil__(),scale)
    printed=parse(selected['G_ball']);dyadic=F(selected['G_lower_dyadic']),F(selected['G_upper_dyadic'])
    assert 0<dyadic[0]<=glo<ghi<=dyadic[1] and printed[0]<=glo<ghi<=printed[1]
    projections=[F(v,144*L) for v in k];delta=min(map(abs,projections));sorted_b=sorted(projections)
    gap=min(b-a for a,b in zip(sorted_b,sorted_b[1:]));M=max(sum(map(abs,z)) for z in Z)/F(144)
    assert delta==F(selected['minimum_abs_projection']) and gap==F(selected['minimum_distinct_projection_gap'])
    radius=min(delta/(2*M),gap/(4*M));assert radius>0
    assert len(set(map(abs,x)))==5 # any two-body-direction span has at most two absolute coordinate values
    assert gate['point']==selected['point'] and list(map(F,gate['G_point_full_dyadic_interval']))==list(dyadic)
    assert F(gate['nonzero_projection_neighbourhood_infinity_radius'])==delta/(2*M)
    gu=dyadic[1];budget=json.loads(resolved[('upper_alpha1_029','check_results.json')].read_text())
    oldbounds=budget['independent_512bit_checks']['I']['independent_interval'];il,iu=F(oldbounds['lower_fraction']),F(oldbounds['upper_fraction'])
    assert [il,iu]==list(map(F,gate['old_complete_I_budget_interval']))
    necessary4=F(32**4,3125)*gu**5;assert necessary4==F(gate['necessary_I_face_for_cost_below32_fourth_power'])
    gatebounds=gate['necessary_I_face_threshold_interval'];a,b=F(gatebounds['lower_fraction']),F(gatebounds['upper_fraction'])
    assert a**4<=necessary4<b**4
    cost4=3125*il**4/gu**5;costbounds=gate['best_possible_cost_using_old_formula_lower_interval'];a,b=F(costbounds['lower_fraction']),F(costbounds['upper_fraction'])
    assert a**4<=cost4<b**4 and a>F(41951,1000)
    residual=3125*il**4-32**4*gu**5;assert residual>0 and residual==F(gate['strict_incompatibility_residual'])
    out={'run_handle':RUN,'pid':pid,'status':'pass_selected018_full_Green_pure_rational_log_and_gate_qualification',
         'point':[str(v) for v in x],'model':'alpha1 Cauchy, psi=sum16 |v.xi|, complete G integral_0^infinity p_t dt',
         'exact_coefficient_checks':{'nonzero_distinct_projections':280,'raw_vertex_coefficients':len(A),
                                     'absolute_projection_count':len(raw),'nonzero_merged_coefficients':len(C),
                                     'per_simplex_and_global_constant_sum_zero':True,'series_count':len(terms),'N':N,
                                     'log_quantization_bits':384,'coefficient_common_denominator_bits':common.bit_length()},
         'pure_Fraction_G_interval':[str(glo),str(ghi)],'inside_source_G_ball_and_dyadic':True,
         'continuity':{'minimum_abs_projection':str(delta),'minimum_gap':str(gap),
                       'own_distinct_nonzero_open_infinity_radius':str(radius),
                       'outside_S2_by_five_distinct_absolute_coordinates':True},
         'gate':{'same_I_face_and_model':True,'essential_cube_g_upper':str(gu),
                 'saved_root_brackets_checked_by_fourth_powers_only':True,'exact_residual_positive':True,
                 'old_budget_best_cost_exceeds_41_951':True,'not_a_lower_bound_for_actual_I_or_C5':True},
         'coefficients_sha256':sha(coefficientpath),'locked_inputs_sha256':locked,
         'source_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.md'),
         'scope':'Only selected new point. No old point,032 anchor,vertex generation,DD/Arb kernel,root search,spatial search or Gamma baseline executed. Geometry and prior complete I budget remain locked inputs.',
         'runtime':{'python':sys.executable,'seconds':time.monotonic()-began,'point_evaluations':1,'transcendental_library_calls':0}}
    (HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({key:out[key] for key in ('run_handle','pid','status','exact_coefficient_checks','pure_Fraction_G_interval','runtime')},indent=2))

if __name__=='__main__':main()
