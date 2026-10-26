#!/usr/bin/env python3
"""Registered period3 mass profiles at the exact C lattice; complete raw inventory."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from collections import Counter
import hashlib,json,os,resource,sys,threading,time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def save(p,d):
    t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(d,separators=(',',':'))+'\n');t.replace(p)
def rt(start):
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss*(1 if sys.platform=='darwin' else 1024),'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
def enc(x):return {f'{i//5},{i%5}':str(c) for i,c in enumerate(x.a) if c}
class Limit(Exception):pass

def main(direction):
    start=time.monotonic();reg=load(HERE/'registration.json')
    if direction not in reg['directions']:reg=load(HERE/'fan_registration.json')
    assert direction in reg['directions']
    assert sha(Path(__file__))==reg['source_sha256']
    for row in reg['inputs']:assert sha(Path(row['path']))==row['sha256']
    pattern=tuple(Q(z) for z in reg['profiles'][direction]);F,u,v,bounds,roots,den=field();nodes=0;masknodes=0
    def budget(which):
        nonlocal nodes,masknodes
        if which=='sweep':nodes+=1
        else:masknodes+=1
        if nodes+masknodes>200000:raise Limit('200000_finite_union_nodes')
    class J:
        __slots__=('c','d','g')
        def __init__(self,c=0,d=0,g=None):
            self.c=c if isinstance(c,F) else F(c);self.d=d if isinstance(d,F) else F(d)
            self.g=tuple(F(0) for _ in range(5)) if g is None else g
        def __hash__(self):return hash((self.c,self.d))
        def __eq__(self,o):
            o=o if isinstance(o,J) else J(o);return self.c==o.c and self.d==o.d
        def __lt__(self,o):
            o=o if isinstance(o,J) else J(o);return self.c<o.c if self.c!=o.c else self.d<o.d
        def __le__(self,o):return self==o or self<o
        def __add__(self,o):
            o=o if isinstance(o,J) else J(o);return J(self.c+o.c,self.d+o.d,tuple(x+y for x,y in zip(self.g,o.g)))
        __radd__=__add__
        def __neg__(self):return J(-self.c,-self.d,tuple(-x for x in self.g))
        def __sub__(self,o):return self+-o
        def __rsub__(self,o):return J(o)+-self
        def __mul__(self,o):
            o=o if isinstance(o,J) else J(o)
            return J(self.c*o.c,self.c*o.d+self.d*o.c,tuple(self.c*y+x*o.c for x,y in zip(self.g,o.g)))
        __rmul__=__mul__
        def __truediv__(self,q):return J(self.c/q,self.d/q,tuple(x/q for x in self.g))
    def encodej(x):return {'value':enc(x.c),'speed':enc(x.d),'faces':[enc(z) for z in x.g]}
    def decodef(d):
        a=[Q(0)]*25
        for ij,c in d.items():i,j=map(int,ij.split(','));a[5*i+j]=Q(c)
        return F(a=a)
    def decodej(d):return J(decodef(d['value']),decodef(d['speed']),tuple(decodef(z) for z in d['faces']))
    a=(1+u**3/2,)+((1+u)/2,)*4
    @lru_cache(None)
    def radius(n,dm):
        n0=n;e2=e3=0
        while n0%2==0:n0//=2;e2+=1
        while n0%3==0:n0//=3;e3+=1
        assert n0==1
        r=u**e2*v**e3/2;assert 32*r**5==F(n)
        return J(r,r*dm/(5*n))
    @lru_cache(None)
    def receiver_interval(k,ell,j,n,dm):
        r=radius(n,dm);aa=J(a[k]);domain=3*aa if k==0 else aa/2
        lo=max(J(0),(j+ell-1)*aa-r);hi=min(domain,j*aa+r)
        def mark(x):return J(x.c,x.d,tuple(x.d if t==k else F(0) for t in range(5)))
        return mark(lo),mark(hi)
    catalog=[];raw=[];counts=Counter();opening=[];zero_discard=[]
    for ls in product(range(1,5),repeat=5):
        n=prod(ls);tail=prod(ls[1:])
        for phase in range(3):
            blocksum=sum(pattern[(phase+t)%3] for t in range(ls[0]));dm=blocksum*tail
            assert n==ls[0]*tail and abs(dm)<=2*n
            catalog.append({'lengths':ls,'phase':phase,'N':n,'mass_derivative':str(dm),'first_axis_block_sum':str(blocksum)})
        for j in range(-1,5):
            counts['raw_candidates']+=1;dm=sum(pattern[(j+t)%3] for t in range(ls[0]))*tail
            endpoints=[];starts=[j]+[-(ell-1)//2 if ell%2 else 1-ell//2 for ell in ls[1:]]
            flags=[]
            for k,ell in enumerate(ls):
                lo,hi=receiver_interval(k,ell,starts[k],n,dm);endpoints.extend((lo,hi))
                flags.append('negative_base' if hi.c<lo.c else 'zero_nonopening' if hi.c==lo.c and hi.d<=lo.d else 'zero_opening' if hi.c==lo.c else 'positive_base')
            if 'negative_base' in flags:counts['discard_negative_base']+=1;continue
            if 'zero_nonopening' in flags:
                counts['discard_zero_nonopening']+=1
                if len(zero_discard)<12:zero_discard.append({'lengths':ls,'start':j,'flags':flags})
                continue
            counts['active_raw']+=1
            if 'zero_opening' in flags:
                counts['active_zero_base_boxes']+=1
                opening.append({'lengths':ls,'start':j,'flags':flags,'mass_derivative':str(dm)})
            # Cartesian indices are distinct; direct first-axis weight count
            # times tail cardinality is the exact common-radius mass slope.
            assert len({(j+t)%3 for t in range(ls[0])})<=3
            raw.append((tuple(endpoints),{'lengths':ls,'first_axis_start':j,'phase':j%3,'starts':starts,'distinct_atom_count':n,'mass_derivative':str(dm)}))
    assert len(catalog)==3072 and counts['raw_candidates']==6144
    unique={b:m for b,m in raw};kept=[]
    for b in unique:
        if not any(o!=b and all(o[2*k]<=b[2*k] and b[2*k+1]<=o[2*k+1] for k in range(5)) for o in unique):kept.append(b)
    kept.sort();counts['unique_active_boxes']=len(unique);counts['kept_jet_boxes']=len(kept)
    events=[]
    for k in range(5):
        groups={}
        for b in kept:
            for x in b[2*k:2*k+2]:groups.setdefault(x.c,set()).add(x.d)
        events.append({'coordinate':k,'distinct_base_endpoint_values':len(groups),'base_values_with_different_speeds':sum(len(s)>1 for s in groups.values())})
    geometry={'direction':direction,'counts':dict(counts),'mass_phase_catalog':catalog,'new_zero_base_opening_boxes':opening,'sample_zero_nonopening_diagnostics':zero_discard,'kept_jets':[dict(unique[b],endpoints=[encodej(x) for x in b]) for b in kept],'endpoint_ties':events,'source_sha256':reg['source_sha256'],'runtime':rt(start)}
    gp=HERE/f'{direction}_geometry.json';save(gp,geometry)
    print(json.dumps({'stage':'complete_raw_jet_geometry','direction':direction,'counts':dict(counts),'endpoint_ties':events,'runtime':rt(start)}),flush=True)
    @lru_cache(None)
    def sweep(boxes):
        budget('sweep')
        if not boxes:return J(0)
        if len(boxes[0])==2:
            answer=J(0);left=right=None
            for lo,hi in sorted(boxes):
                if left is None:left,right=lo,hi
                elif right<lo:answer+=right-left;left,right=lo,hi
                else:right=max(right,hi)
            return answer+right-left
        cuts=sorted({x for b in boxes for x in b[:2]});out=J(0)
        for lo,hi in zip(cuts,cuts[1:]):
            sub=tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]}))
            out+=(hi-lo)*sweep(sub)
        return out
    cuts=sorted({x for b in kept for x in b[:2]});prefix=0;volume=J(0);cp=HERE/f'{direction}_checkpoint.json'
    cut_hash=hashlib.sha256(json.dumps([encodej(x) for x in cuts],sort_keys=True).encode()).hexdigest()
    if cp.exists():
        ck=load(cp);assert ck['source_sha256']==reg['source_sha256'] and ck['first_axis_cut_hash']==cut_hash
        prefix=ck['completed_first_axis_slabs'];volume=decodej(ck['partial_volume_jet'])
    def checkpoint(status):
        save(cp,{'status':status,'direction':direction,'completed_first_axis_slabs':prefix,'total_first_axis_slabs':len(cuts)-1,'partial_volume_jet':encodej(volume),'first_axis_cut_hash':cut_hash,'source_sha256':reg['source_sha256'],'sweep_nodes':nodes,'mask_nodes':masknodes,'runtime':rt(start),'no_partial_derivative_claim':status!='COMPLETE_two_union_engines'})
    try:
        for index in range(prefix,len(cuts)-1):
            lo,hi=cuts[index:index+2];sub=tuple(sorted({b[2:] for b in kept if b[0]<=lo and hi<=b[1]}))
            volume+=(hi-lo)*sweep(sub);prefix=index+1;checkpoint('FIRST_AXIS_exact_jet_prefix_saved')
        segments=[]
        for k in range(5):
            axis=sorted({x for b in kept for x in b[2*k:2*k+2]});ss=[]
            for lo,hi in zip(axis,axis[1:]):
                mask=sum(1<<i for i,b in enumerate(kept) if b[2*k]<=lo and hi<=b[2*k+1])
                if mask:ss.append((hi-lo,mask))
            segments.append(ss)
        @lru_cache(None)
        def masks(k,mask):
            budget('mask')
            if k==5:return J(1)
            return sum((d*masks(k+1,mask&nextmask) for d,nextmask in segments[k] if mask&nextmask),J(0))
        second=masks(0,(1<<len(kept))-1)
        assert volume==second and volume.g==second.g
        ratio=volume*Q(16,3)
        C=load(ROOT/'lower_round003/a2_equals_a3_results.json');assert ratio.c==decodef(C['coefficients'])
        assert sum(ratio.g,F(0))==ratio.d
        dlo,dhi=bounds(ratio.d.a)
        status='STRICT_POSITIVE_DIRECTIONAL_DERIVATIVE' if dlo>0 else 'STRICT_NEGATIVE_DIRECTIONAL_DERIVATIVE' if dhi<0 else 'EXACT_ZERO_DIRECTIONAL_DERIVATIVE' if ratio.d==F(0) else 'UNRESOLVED_ROOT_ISOLATION'
        result={'status':status,'direction':direction,'pattern':list(map(str,pattern)),'parameter':'epsilon->0+;registered finite profile','period_mass':'3','period_mass_derivative':'0','receiver_domain':'[0,3a1] x [0,s/2]^4','ratio_normalization':'16*physical_union_volume/3;no cellvolume denominator','baseline_inherited_from_full_C_by_translation_and_reflection':True,'jet_constant_consistency_equals_full_C_coefficients':True,'raw_length_phase_count':len(catalog),'counts':dict(counts),'directional_derivative_coefficients':enc(ratio.d),'derivative_interval':[str(dlo),str(dhi)],'derivative_decimal_display_only':[float(dlo),float(dhi)],'synchronous_coordinate_face_contributions':[{'coordinate':k,'coefficients':enc(x),'interval':list(map(str,bounds(x.a)))} for k,x in enumerate(ratio.g)],'endpoint_ties':events,'union_engines_agree_constant_derivative_and_faces':True,'first_axis_slabs':len(cuts)-1,'coordinate_segments':list(map(len,segments)),'sweep_nodes':nodes,'mask_nodes':masknodes,'source_sha256':reg['source_sha256'],'registration_sha256':sha(HERE/('registration.json' if direction=='middle' else 'fan_registration.json')),'geometry_sha256':sha(gp),'C_source_sha256':sha(ROOT/'lower_round003/a2_equals_a3_results.json'),'root_denominator':str(den),'root_lower_numerators':roots,'runtime':rt(start),'scope':'registered period3 profile or fan cell;finite1024length family only;no finite epsilon or C5 improvement certified;main owns proof/L1/collar'}
        save(HERE/f'{direction}_derivative.json',result);checkpoint('COMPLETE_two_union_engines')
        print(json.dumps({k:v for k,v in result.items() if k not in ('directional_derivative_coefficients','derivative_interval','synchronous_coordinate_face_contributions')},indent=2),flush=True)
    except Limit as e:
        checkpoint('CHECKPOINT_FINITE_NODE_LIMIT')
        save(HERE/f'{direction}_limit.json',{'status':'CHECKPOINT_FINITE_NODE_LIMIT','reason':str(e),'no_partial_derivative_claim':True,'runtime':rt(start)})
    except Exception as e:
        checkpoint('ACTUAL_FAULT_no_automatic_restart')
        save(HERE/f'{direction}_fault.json',{'status':'ACTUAL_FAULT_no_automatic_restart','type':type(e).__name__,'message':str(e),'runtime':rt(start)});raise

if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    if len(sys.argv)!=2:raise SystemExit('usage: profile_variation.py registered_profile_id')
    main(sys.argv[1])
