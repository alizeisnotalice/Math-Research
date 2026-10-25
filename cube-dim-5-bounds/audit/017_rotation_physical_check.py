#!/usr/bin/env python3
"""Independent one-block K8/n4 physical-capture audit, standard library only."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib,json,time

D=2**80; G=2**40

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def root_bracket(n):
    target=n*D**5; low,high=0,1
    while high**5<=target: high*=2
    while high-low>1:
        mid=(low+high)//2
        if mid**5<=target: low=mid
        else: high=mid
    if low**5==target: low-=1
    assert low**5<target<=(low+1)**5
    return Q(low,D),Q(low+1,D)
def scale(c,b): return (c*b[0],c*b[1]) if c>=0 else (c*b[1],c*b[0])
def floor_grid(q): return Q((q*G).numerator//(q*G).denominator,G)
def ceil_grid(q): return -floor_grid(-q)

def main():
    parser=ArgumentParser(); parser.add_argument('--project-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args(); project=args.project_root.resolve(); out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    source=project/'lower_round016'; pp=source/'K8_n4_partition.json'; rp=source/'K8_n4_results.json'
    result=json.loads(rp.read_text()); packet=json.loads(pp.read_text())
    registration={'status':'REGISTERED_one_independent_small_block_before_computation',
        'scope':'K=8,n=4 only; no global 26-family rerun; original integer m,n enumeration',
        'limits':{'seconds':120,'elementary_cells':10000},'producer_packet_sha256':sha(pp),
        'producer_result_sha256':sha(rp),'audit_source_sha256':sha(Path(__file__))}
    (out/'017_rotation_physical_registration.json').write_text(json.dumps(registration,indent=2)+'\n')
    started=time.monotonic()
    def cap():
        if time.monotonic()-started>120: raise TimeoutError('audit 017 120-second limit')
    # Direct original-lattice enumeration. Inverse bounds imply all points
    # in the selected physical coefficient square have m,n in [-20,20].
    direct={}
    for m,n in product(range(-20,21),repeat=2):
        xy=(Q(4*m-3*n,5),Q(3*m+4*n,5))
        if all(-5<=c<10 for c in xy): direct[m,n]=xy
    reps={mn:xy for mn,xy in direct.items() if all(0<=c<5 for c in xy)}
    assert len(reps)==25 and len(direct)==225 and len(set(direct.values()))==225
    # Exact inverse for any xcoef,ycoef in [-5,10) gives m in [-7,14),
    # n in (-10,11), hence enumeration [-20,20] is complete.
    byid={}
    for record in packet['atom_manifest']:
        mn=(record['m'],record['n']); ident=record['id']
        assert mn in direct and ident not in byid
        assert direct[mn]==(Q(record['x_coefficient_of_s']),Q(record['y_coefficient_of_s']))
        byid[ident]=(mn,direct[mn])
    assert len(byid)==225 and {p[0] for p in byid.values()}==set(direct)
    ul,uh=root_bracket(2); ss=((1+ul)/2,(1+uh)/2); aa=(1+ul**3/2,1+uh**3/2)
    rl,rh=root_bracket(32); rl/=2; rh/=2
    assert rh==1 and (2*rl)**5<32 and str(rl)==result['strict_inner_radius']
    all_windows=[]; unrot=[]; areas=[]; diagnostics=[]
    for mode in (0,1):
        windows={}
        for ident,(_,xy) in byid.items():
            box=[]
            for c in xy:
                lower,upper=scale(c,ss)
                if mode==0: lo=ceil_grid(max(0,upper-rl)); hi=floor_grid(min(5*ss[0],lower+rl))
                else: lo=floor_grid(max(0,lower-rh)); hi=ceil_grid(min(5*ss[1],upper+rh))
                box.extend((lo,hi))
            if box[0]<box[1] and box[2]<box[3]: windows[ident]=tuple(box)
        all_windows.append(windows)
        ub=[]
        for ab in (aa,ss,ss):
            if mode==0: ub.extend((ceil_grid(max(0,ab[1]-rl)),floor_grid(min(ab[0]/2,rl))))
            else: ub.extend((floor_grid(max(0,ab[0]-rh)),ceil_grid(min(ab[1]/2,rh))))
        ub=tuple(ub);unrot.append(ub)
        assert ub==tuple(map(Q,packet['unrotated_inner_outer_boxes'][mode]))
        xs=sorted({v for b in windows.values() for v in b[:2]})
        ys=sorted({v for b in windows.values() for v in b[2:]})
        assert (len(xs)-1)*(len(ys)-1)<=10000
        counts=Counter();area=Q(0);coverage={}; qualifying=0
        for i,(xl,xh) in enumerate(zip(xs,xs[1:])):
            cap(); xc=(xl+xh)/2
            for j,(yl,yh) in enumerate(zip(ys,ys[1:])):
                yc=(yl+yh)/2
                selected=tuple(ident for ident,b in windows.items() if b[0]<xc<b[1] and b[2]<yc<b[3])
                counts[len(selected)]+=1; coverage[i,j]=len(selected)>=4
                if len(selected)>=4: area+=(xh-xl)*(yh-yl);qualifying+=1
        # Independent disjoint-cell audit of the compressed packet, including
        # checks that every rectangle keeps the SAME four IDs throughout.
        key='inner_fixed_atom_rectangles' if mode==0 else 'outer_count_rectangles'
        idkey='atom_ids' if mode==0 else 'covering_window_ids'
        occupied=Counter();packet_area=Q(0);guard_count=0;nonrect=0
        for record in packet[key]:
            b=tuple(map(Q,record['rectangle'])); ids=record[idkey]
            assert len(ids)==len(set(ids))==4
            for ident in ids:
                w=windows[ident]
                assert w[0]<=b[0]<b[1]<=w[1] and w[2]<=b[2]<b[3]<=w[3]
            if mode==0:
                # Rebuild 32 distinct full five-integer atom indices.
                pts={(j1,*byid[ident][0],j4,j5) for j1,j4,j5 in product((0,1),repeat=3) for ident in ids}
                assert len(pts)==32
                extremes=[((Q(0),Q(0)),ab) for ab in (aa,ss,ss)]
                for axis in (0,1):
                    cs=[byid[ident][1][axis] for ident in ids]
                    extremes.append((scale(min(cs),ss),scale(max(cs),ss)))
                full=ub+b
                for axis,(mn,mx) in enumerate(extremes):
                    assert mx[1]-rl<=full[2*axis] and full[2*axis+1]<=mn[0]+rl
                    guard_count+=2
                ms={byid[k][0][0] for k in ids}; ns={byid[k][0][1] for k in ids}
                if len(ms)*len(ns)!=4: nonrect+=1
            x0,x1=xs.index(b[0]),xs.index(b[1]); y0,y1=ys.index(b[2]),ys.index(b[3])
            for i,j in product(range(x0,x1),range(y0,y1)): occupied[i,j]+=1
            packet_area+=(b[1]-b[0])*(b[3]-b[2])
        assert all(occupied[cell]==int(good) for cell,good in coverage.items())
        assert packet_area==area==Q(result['area_interval'][mode])
        vol=(ub[1]-ub[0])*(ub[3]-ub[2])*(ub[5]-ub[4])
        assert Q(8,25)*vol*area==Q(result['subproblem_ratio_interval'][mode])
        areas.append(area)
        diagnostics.append({'active_windows':len(windows),'segments':[len(xs)-1,len(ys)-1],
            'qualifying_cells':qualifying,'coverage_histogram':dict(counts),'compressed_rectangles':len(packet[key]),
            'geometry_endpoint_guards':guard_count,'nonrectangular_witnesses':nonrect,
            'packet_is_disjoint_exact_cover_of_threshold_cells':True})
    assert areas[0]<=areas[1] and sha(pp)==result['partition_sha256']
    output={'status':'PASS_independent_rotation_physical_K8_n4_mechanism_and_packet',
        'scope':'one independently rebuilt small block; no certification of all 26 families or global new C',
        'producer_imports':False,'direct_original_integer_range':[-20,20],
        'supercell_atoms':len(reps),'neighbor_atoms':len(direct),'K':8,'n':4,'distinct_mass_per_witness':32,
        'exact_common_radius':'1','strict_inner_radius':str(rl),'ratio_normalization':'8*V/25',
        'area_interval':[str(t) for t in areas],'ratio_interval':result['subproblem_ratio_interval'],
        'inner':diagnostics[0],'outer':diagnostics[1],
        'root_spacing_denominator':D,'coordinate_grid_denominator':G,
        'elapsed_seconds':time.monotonic()-started,'limits':registration['limits'],
        'producer_packet_sha256':sha(pp),'producer_result_sha256':sha(rp),
        'producer_source_sha256':sha(source/'physical_capture.py'),'audit_source_sha256':sha(Path(__file__))}
    (out/'017_rotation_physical_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))

if __name__=='__main__':main()
