#!/usr/bin/env python3
"""Low-resource certificate audit; one side4 profile, no global union rebuild."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib,json,os,resource,sys,threading,time

ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
D=2**80; G=2**40
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def root2():
    lo,hi=0,2*D
    while hi-lo>1:
        m=(lo+hi)//2
        if m**5<2*D**5:lo=m
        else:hi=m
    assert lo**5<2*D**5<hi**5
    return Q(lo,D),Q(hi,D)
def floor_grid(q):return (q*G).numerator//(q*G).denominator
def ceil_grid(q):return -floor_grid(-q)
def main():
    ap=ArgumentParser();ap.add_argument('--project-root',required=True,type=Path)
    args=ap.parse_args();assert os.getpriority(os.PRIO_PROCESS,0)>=10
    assert all(os.environ.get(k)=='1' for k in ENV) and threading.active_count()==1
    start=time.monotonic();p=args.project_root.resolve();b=p/'lower_round016/delivery007';out=Path(__file__).resolve().parent
    def load(name):return json.loads((b/name).read_text())
    manifest=load('manifest.json');bound_inputs={}
    for row in manifest['files']:
        assert sha(b/row['path'])==row['sha256'];bound_inputs[row['path']]=row['sha256']
    ref=manifest['existing_reference_dependency'];assert sha(b/ref['path'])==ref['sha256']
    registration={'scope':'saved arithmetic/certificate checks; ONE independent side4 profile; no full geometry union rebuild',
                  'producer_imports':False,'nice_min':10,'threads':1,'limits_seconds':120,
                  'locked_manifest_sha256':sha(b/'manifest.json'),'source_sha256':sha(Path(__file__))}
    (out/'registration.json').write_text(json.dumps(registration,indent=2)+'\n')
    setup=load('profile_setup_v4.json');cert=load('radius_certificate.json');bundle=load('geometry_bundle.json')
    hist=load('aggregate_histograms.json');result=load('result.json');evidence=load('case_evidence.json')
    for name,h in result['input_hashes'].items():assert sha(b/name)==h
    uL,uH=root2();aa=(1+uL**3/2,1+uH**3/2);ss=((1+uL)/2,(1+uH)/2)
    assert list(map(Q,setup['spacings_intervals']['a']))==list(aa)
    assert list(map(Q,setup['spacings_intervals']['s']))==list(ss)
    assert aa[0]>Q(7,4) and ss[0]>1
    F8=(1+aa[1]/8)*(1+ss[1]/8)**2*(1+Q(14,5)*ss[1]/8+ss[1]**2/64)/(aa[0]*ss[0]**4)
    assert F8==Q(setup['F8_rigorous_upper']) and F8<1
    def iv(coeff):return tuple(coeff[0]+coeff[1]*a+coeff[2]*s for a,s in zip(aa,ss))
    crit={(Q(4),Q(0),Q(0)),(Q(8),Q(0),Q(0))}
    for axis in (1,2):
        for k in range(1,101):
            c=[Q(0)]*3;c[axis]=Q(k,5) if axis==2 else Q(k)
            interval=iv(c)
            if 4<interval[0] and interval[1]<8:crit.add(tuple(c))
    ordered=sorted(crit,key=lambda c:iv(c)[0]);points=setup['critical_points']
    assert len(points)==len(ordered)==23
    for i,(c,point) in enumerate(zip(ordered,points)):
        assert tuple(map(Q,point['expression_coefficients']))==c
        assert tuple(map(Q,point['interval']))==iv(c)
        if i:assert Q(points[i-1]['interval'][1])<Q(point['interval'][0])
    assert len(cert['interval_results'])==setup['open_intervals']==22
    assert cert['setup_sha256']==sha(b/'profile_setup_v4.json')
    interval_summary=[]
    for i,d in enumerate(cert['interval_results']):
        assert d['interval_index']==i and d['left']==points[i] and d['right']==points[i+1]
        lb=tuple(map(Q,points[i]['interval']));rb=tuple(map(Q,points[i+1]['interval']))
        side=tuple((x+y)/2 for x,y in zip(lb,rb));assert tuple(map(Q,d['interior_side_interval']))==side
        assert lb[1]<side[0]<=side[1]<rb[0]
        ms=[]
        for spacing in (aa,ss,ss):
            valid=[k for k in range(10) if k*spacing[1]<side[0] and side[1]<(k+1)*spacing[0]]
            assert len(valid)==1;ms.append(valid[0]+1)
        assert ms==d['max_unrotated_counts']
        count=d['capture']['exact_maximum']
        assert d['capture']['inner_outer_maxima']==[count,count] and d['capture']['dual_cell_and_event_algorithms_agree']
        mass=ms[0]*ms[1]*ms[2]*count
        assert mass==d['max_total_mass'] and Q(mass)<lb[0]**5
        assert d['excluded_for_all_sides_in_open_interval'] and not d['large_side_cover']
        interval_summary.append({'index':i,'planar_maximum_saved_production':count,'total_mass':mass,'strict_lower_endpoint_fifth_power_margin':str(lb[0]**5-mass)})
    # The only independently rebuilt planar arrangement in this audit: T=4.
    direct={}
    for m,n in product(range(-20,21),repeat=2):
        xy=(Q(4*m-3*n,5),Q(3*m+4*n,5))
        if all(-5<=c<10 for c in xy):direct[m,n]=xy
    atoms=load('rotation_atoms.json');assert len(direct)==225
    assert {(a['m'],a['n']) for a in atoms}==set(direct)
    for a in atoms:assert direct[a['m'],a['n']]==(Q(a['x_coefficient_of_s']),Q(a['y_coefficient_of_s']))
    maxima=[];cells=[]
    for mode in (0,1):
        boxes=[]
        for xy in direct.values():
            box=[]
            for c in xy:
                low,high=(c*ss[0],c*ss[1]) if c>=0 else (c*ss[1],c*ss[0])
                if mode==0:lo=ceil_grid(max(0,high-2));hi=floor_grid(min(5*ss[0],low+2))
                else:lo=floor_grid(max(0,low-2));hi=ceil_grid(min(5*ss[1],high+2))
                box.extend((lo,hi))
            if box[0]<box[1] and box[2]<box[3]:boxes.append(box)
        xs=sorted({x for a in boxes for x in a[:2]});ys=sorted({y for a in boxes for y in a[2:]})
        maximum=0
        for x0,x1 in zip(xs,xs[1:]):
            for y0,y1 in zip(ys,ys[1:]):
                count=sum(2*a[0]<x0+x1<2*a[1] and 2*a[2]<y0+y1<2*a[3] for a in boxes)
                maximum=max(maximum,count)
        maxima.append(maximum);cells.append((len(xs)-1)*(len(ys)-1))
    assert maxima==setup['side4_overlap']['inner_outer_maxima']==[15,15]
    Ks=sorted({a*b*c for a,b,c in product(range(1,5),repeat=3)})
    assert {(d['K'],d['n']) for d in evidence['cases']}==set(product(Ks,range(1,17)))
    assert len(evidence['cases'])==256 and len(bundle['families'])==26
    positive={(d['K'],d['n']) for d in evidence['cases'] if Q(d['area_interval'][1])>0}
    assert positive=={(f['K'],f['n']) for f in bundle['families']}
    # Exact re-evaluation of all four SAVED production histogram measures.
    # The histograms themselves are NOT reconstructed from geometry here.
    assert hist['grid']==bundle['grid']==G
    assert {(h['kind'],h['mode']) for h in hist['histograms']}==set(product(('plane','unrot'),(0,1)))
    answer=[];prior2=[];prior3=[];hist_details=[]
    for mode in (0,1):
        def get(kind):return next(h for h in hist['histograms'] if h['kind']==kind and h['mode']==mode)
        ph,bh=get('plane'),get('unrot');assert ph['power']==2 and bh['power']==3
        joint=defaultdict(int)
        for pm,a in ph['integer_units'].items():
            assert a>0 and 0<int(pm)<1<<26
            for bm,v in bh['integer_units'].items():
                assert v>0 and 0<int(bm)<1<<26
                mask=int(pm)&int(bm)
                if mask:joint[mask]+=a*v
        factor=Q(8,25*G**5);answer.append(factor*sum(joint.values()))
        for included,target in (({1,8},prior2),({1,2,8},prior3)):
            mask=sum(1<<i for i,f in enumerate(bundle['families']) if f['K'] in included)
            target.append(factor*sum(v for m,v in joint.items() if m&mask))
        for kind,h in (('plane',ph),('unrot',bh)):
            axes=sorted({b['box'][k] for f in bundle['families'] for b in (f['inner_planar'] if mode==0 else f['outer_planar']) for k in (0,1)}) if kind=='plane' else sorted({b[k] for f in bundle['families'] for b in (f['inner_boxes'] if mode==0 else f['outer_boxes']) for k in (0,1)})
            blocks=[d for d in hist['checkpoint_evidence'] if d['kind']==kind and d['mode']==mode]
            last=0
            for d in blocks:
                assert d['x_indices'][0]==last;last=d['x_indices'][1]
                assert d['x_grid_bounds']==[axes[d['x_indices'][0]],axes[last]]
            assert last==len(axes)-1
            hist_details.append({'kind':kind,'mode':mode,'mask_bins':len(h['integer_units']),'production_block_certificates':len(blocks),'block_partition_indices_verified':True,'independently_rebuilt':False})
    assert answer==list(map(Q,result['ratio_interval'])) and answer[0]<answer[1]
    assert prior2==list(map(Q,result['prior_signature_ratios']['K1_K8']))
    assert prior3==list(map(Q,result['prior_signature_ratios']['K1_K2_K8']))
    reference=json.loads((b/ref['path']).read_text());assert answer[1]<Q(reference['ratio_interval'][0])
    assert time.monotonic()-start<120
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    output={'status':'PASS_limited_independent_radius_logic_and_saved_rollup_arithmetic',
        'scope':'one fixed rotation/spacings at level1 all scales; NOT a general C5 upper bound or other-threshold exclusion',
        'global_histograms_independently_rebuilt':False,'all_22_profiles_independently_rebuilt':False,
        'independent_profile':{'side':4,'inner_outer_maxima':maxima,'elementary_cells_each_mode':cells},
        'F8_outward_upper':str(F8),'F8_display':float(F8),'critical_points':len(points),'strictly_excluded_intervals':interval_summary,
        'small_side_actual_count_bounds':[3,4,4,15],'small_side_mass_bound':720,
        'histograms':hist_details,'ratio_interval':[str(q) for q in answer],
        'ratio_display':[float(q) for q in answer],'strictly_below_complete_C':True,
        'production_inputs':bound_inputs,'manifest_sha256':sha(b/'manifest.json'),
        'source_sha256':sha(Path(__file__)),'runtime':{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),
        'thread_environment':{k:os.environ[k] for k in ENV},'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,'wall_seconds':time.monotonic()-start}}
    (out/'results.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps({k:v for k,v in output.items() if k not in ('production_inputs','strictly_excluded_intervals')},indent=2))
if __name__=='__main__':main()
