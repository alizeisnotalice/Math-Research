"""Independent lightweight audit, no union batch and no source mutation.
Masses are counted from the infinite periodic lattice directly, not a finite
window or the author's prefix sum. Only positive stored receivers are audited.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib,json,math,os,platform,time,argparse

D=Path(__file__).resolve().parent;P=D.parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,default=D/'generated');S=parser.parse_args().output_dir.resolve()
start=time.time()
z=json.loads((S/'decomposition.json').read_text())
old=json.loads((P/'checks/lower/exactinputs/c2lower_1p70272.json').read_text())
new=json.loads((P/'round002/verification/candidate002.json').read_text())
assert old['pts']==new['pts'][:-1] and len(new['pts'])==51
assert all(old[k]==new[k] for k in ('Px','Py','reach','lam')) and F(new['lam'])==1
dx,dm,dr,L=(int(z[k]) for k in ('Dx','Dm','Dr','L'))
px,py,reach=(F(new[k]) for k in ('Px','Py','reach'))
atoms=[tuple(map(F,p)) for p in new['pts']]
assert dx==math.lcm(*(v.denominator for p in atoms for v in p[:2]),px.denominator,py.denominator,reach.denominator)
assert dm==math.lcm(*(p[2].denominator for p in atoms)) and L==math.lcm(dx,dr) and dr==2**40
atoms=[(int(x*dx),int(y*dx),int(m*dm)) for x,y,m in atoms]
pxi,pyi,reachi=int(px*dx),int(py*dx),int(reach*dx)
oldx={x%pxi for x,y,m in atoms[:-1]};oldy={y%pyi for x,y,m in atoms[:-1]}
assert atoms[-1][0]%pxi not in oldx and atoms[-1][1]%pyi not in oldy
ceildiv=lambda a,b:-((-a)//b)
def count(a,b,x,p):return max(0,(b-x)//p-ceildiv(a-x,p)+1)
counts=Counter();max_radius=0
for row in z['records']:
    a,b,c,d,m0,m1=map(int,row[:6]);assert 0<=b-a<=reachi and 0<=d-c<=reachi
    masses=[count(a,b,x,pxi)*count(c,d,y,pyi)*m for x,y,m in atoms]
    assert sum(masses[:-1])==m0 and sum(masses)==m1
    oldbbox=all(v%pxi in oldx for v in (a,b)) and all(v%pyi in oldy for v in (c,d))
    def receiver(m):
        if m<=0:return None
        k=math.isqrt(m*dr*dr//(4*dm))-1
        if k<=0:return None
        assert F(k*k,dr*dr)<F(m,4*dm)
        r=k*(L//dr)
        R=(max(b*(L//dx)-r,0),min(a*(L//dx)+r,pxi*(L//dx)),max(d*(L//dx)-r,0),min(c*(L//dx)+r,pyi*(L//dx)))
        return R if R[0]<R[1] and R[2]<R[3] else None
    r0=receiver(m0) if oldbbox else None;r1=receiver(m1)
    cls='old_positive' if oldbbox and r0 else 'old_empty_new_positive' if oldbbox else 'new_boundary_positive'
    assert r1==tuple(map(int,row[7])) and r0==(tuple(map(int,row[8])) if row[8] else None) and cls==row[6]
    if r0:assert r1[0]<=r0[0]<r0[1]<=r1[1] and r1[2]<=r0[2]<r0[3]<=r1[3]
    counts[cls]+=1;max_radius=max(max_radius,math.isqrt(m1*dr*dr//(4*dm))-1)
assert all(counts[k]==v for k,v in z['counts'].items())
A={k:F(v) for k,v in z['areas'].items()};m0=sum(F(p[2]) for p in old['pts']);m1=sum(F(p[2]) for p in new['pts'])
assert m1-m0==F(1,1000)
for key,stem,m in [('A0','checks/lower/exactinputs/c2lower_1p70272',m0),('Afull','round002/verification/candidate002',m1)]:
    r=json.loads((P/(stem+'.check2.json')).read_text());assert A[key]/m==F(int(r['ratio_num']),int(r['ratio_den']))
assert A['A0']<=A['A1']==A['A12']<A['Afull']
assert A['A1']/m1<A['A0']/m0<A['Afull']/m1
assert F(z['increments']['class1'])==A['A1']-A['A0']
assert F(z['increments']['class2_given_class1'])==0
assert F(z['increments']['class3_given_1_and_2'])==A['Afull']-A['A12']
# A narrowly targeted extra check: a pure new-point rectangle is insufficient
# even after adding it to all updated old-coordinate boxes. No global union run.
zero_old=list({tuple(map(int,row[7])) for row in z['records'] if int(row[4])==0})
target=max(zero_old,key=lambda r:(r[1]-r[0])*(r[3]-r[2]))
assert all(target[0]<=r[0]<=r[1]<=target[1] and target[2]<=r[2]<=r[3]<=target[3] for r in zero_old)
old_rectangles=[tuple(map(int,row[7])) for row in z['records'] if row[6]=='old_positive']
clips=[(max(target[0],r[0]),min(target[1],r[1]),max(target[2],r[2]),min(target[3],r[3])) for r in old_rectangles]
clips=[r for r in clips if r[0]<r[1] and r[2]<r[3]]
xs=sorted({v for r in clips for v in r[:2]});covered=0
for a,b in zip(xs,xs[1:]):
    ys=sorted((r[2],r[3]) for r in clips if r[0]<=a and b<=r[1]);length=0;end=None
    for c,d in ys:
        length+=d-c if end is None else max(0,d-max(c,end))
        end=d if end is None else max(end,d)
    covered+=(b-a)*length
pure_gain=F((target[1]-target[0])*(target[3]-target[2])-covered,L*L)
assert (A['A1']+pure_gain)/m1<A['A0']/m0
old_base_rectangles=[tuple(map(int,row[8])) for row in z['records'] if row[8] is not None]
no_new_mass=[tuple(map(int,row[7])) for row in z['records'] if row[6]=='new_boundary_positive' and int(row[5])==int(row[4])]
assert all(any(r[0]<=q[0]<=q[1]<=r[1] and r[2]<=q[2]<=q[3]<=r[3] for r in old_base_rectangles) for q in no_new_mass)
files=list(S.glob('*'))+[P/'checks/lower/exactinputs/c2lower_1p70272.json',P/'round002/verification/candidate002.json',P/'checks/lower/exact_check.py',P/'checks/lower/lower_check.py']
bindings={(str(f.relative_to(P)) if f.is_relative_to(P) else 'generated/'+f.name):hashlib.sha256(f.read_bytes()).hexdigest() for f in files if f.is_file()}
out={'pid':os.getpid(),'python':platform.python_version(),'seconds':time.time()-start,'positive_receivers':len(z['records']),'counts':dict(counts),'mass0':str(m0),'mass1':str(m1),'ratios':{k:str(A[k]/(m0 if k=='A0' else m1)) for k in A},'areas':{k:str(v) for k,v in A.items()},'no_area_batch_rerun':True,'zero_activation_exhaustiveness':'source loop audited, not independently re-enumerated','infinite_periodic_direct_mass_and_all_receiver_checks':True,'sha256':bindings}
out['pure_new_point_check']={'zero_old_records':sum(int(r[4])==0 for r in z['records']),'distinct_zero_old_rectangles':len(zero_old),'positive_intersections_with_A1_records':len(clips),'net_gain':str(pure_gain),'A1_plus_pure_point_ratio':str((A['A1']+pure_gain)/m1),'below_old_ratio':True}
out['pure_new_point_check']['new_boundary_no_new_mass_records']=len(no_new_mass)
out['pure_new_point_check']['all_no_new_mass_rectangles_contained_in_single_A0_record']=True
out['pure_new_point_check']['mixed_new_boundary_records']=sum(row[6]=='new_boundary_positive' and 0<int(row[4])<int(row[5]) for row in z['records'])
(S/'mechanism_checks.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='sha256'},indent=2))
