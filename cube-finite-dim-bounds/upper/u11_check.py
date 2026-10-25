from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import hashlib,json,time
import sympy as s
import mpmath as mp
import u11_green as g

start=time.time()
R=g.BASE_RAYS
rates={k:sum(abs(g.dot(v,r)) for v in g.V) for k,r in R.items()}
assert rates==dict(A=8,B=8,C=12,D=12,E=16,F=24)
base=[]
for labs,det in g.BASE_CONES:
 mat=s.Matrix.hstack(*[s.Matrix(R[z]) for z in labs]);assert abs(mat.det())==det
 norm=[sum(R[z]) for z in labs]
 base.append({'labels':labs,'det':det,'sum_section_det':str(Q(det,1)/g.prod(norm)),
              'spectral_simplex_volume':str(Q(det,24)/g.prod(rates[z] for z in labs))})
assert sum(Q(z['sum_section_det']) for z in base)==Q(1,24)
assert len(g.RAYS)==240 and len(g.CONES)==1152
assert len({tuple(sorted(ids)) for ids,d in g.CONES})==1152
qr=g.rates();specvol=sum(Q(det,24)/g.prod(qr[i] for i in ids) for ids,det in g.CONES)
p0coeff=sum(Q(det)/g.prod(qr[i] for i in ids) for ids,det in g.CONES)/16
assert specvol==Q(11,4608) and p0coeff==Q(11,3072)
rankzero={}
for name,r in R.items():
 rows=[list(v) for v in g.V if g.dot(v,r)==0]
 rankzero[name]=s.Matrix(rows).rank() if rows else 0
assert rankzero==dict(A=0,B=3,C=0,D=3,E=2,F=1)
vertices={tuple(Q(a)/q for a in r) for r,q in zip(g.RAYS,qr) if sorted(map(abs,r)) in ([0,0,1,1],[1,1,1,1])}
assert len(vertices)==40
# Independent scalar time-integral normalizations, fixed node lists only.
mp.mp.dps=60
scalar=[]
for nodes,expected in [((-1,-2,1,2),mp.pi/12),((-1,-1,1,1),mp.pi/4)]:
 val=mp.quad(lambda t:mp.re(1/g.prod(t-1j*b for b in nodes)),[0,1,mp.inf])
 dd=-mp.pi/2*g._mpq(g.dd_sign(nodes))
 assert abs(val-expected)<mp.mpf('1e-50') and abs(dd-expected)<mp.mpf('1e-50')
 scalar.append({'nodes':nodes,'quadrature':mp.nstr(val,28),'DD_value':mp.nstr(dd,28)})
eps=mp.mpf(1)/8
v=mp.re(-1j*g.dd_cutoff_log((0,0,1,-1),eps))
expected=1/eps-mp.pi/2+mp.atan(eps)
assert abs(v-expected)<mp.mpf('1e-50')
scalar.append({'nodes':[0,0,1,-1],'epsilon':'1/8','cutoff_formula':mp.nstr(v,28),'expected':mp.nstr(expected,28)})
bases=[(1,2,3,5),(1,3,4,7),(2,3,5,11),(1,2,4,9),(2,5,7,13),(3,4,8,15),
       (1,4,6,10),(2,3,7,12),(1,5,8,14),(3,5,9,17),(2,6,11,19),(4,7,13,23)]
primes=(151,157,163,167);points=[tuple(Q(a)+Q(1,p) for a,p in zip(base,primes)) for base in bases]
rows=[]
for j,x in enumerate(points):
 xx,qq,b=g.point_nodes(x)
 assert all(z for z in b) and all(len({b[i] for i in ids})==4 for ids,det in g.CONES)
 co=g.green_pi3_coefficient(x);assert co>0
 assert g.green_allnodes_pi3_coefficient(x)==co
 res=g.generic_residues(x)
 reco=sum(c/bj for bj,c in res.items())/32
 assert reco==co
 exact=g._mpq(co)/mp.pi**3
 cut=g.green_cutoff_mp(x,mp.mpf(2)**-40,dps=80)
 # Arbitrary-precision diagnostics are not roundoff intervals.
 assert abs(cut-exact)<mp.mpf('1e-35')
 t=mp.mpf(1)/3
 pm=g.density_mp_from_residues(t,res);pn=g.density_float(float(t),list(map(float,x)))
 assert pm>0 and abs(pn-float(pm))<=1e-8*float(pm)
 fl=g.green_float_generic(list(map(float,x)))
 item={'id':j+1,'x':list(map(str,x)),'G_pi3_exact':str(co),'G_display':mp.nstr(exact,28),
       'cutoff_2^-40_display':mp.nstr(cut,28),'cutoff_diagnostic_absdiff':mp.nstr(abs(cut-exact),8),
       'float_generic_relative_error':float(abs(mp.mpf(fl)-exact)/exact),
       'density_t1/3_relative_error':float(abs(mp.mpf(pn)-pm)/pm),'nonzero_residue_pairs':len(res)}
 if j<2:
  qv=g.green_time_quad(x,dps=50)
  rel=abs(qv-exact)/exact
  assert rel<mp.mpf('1e-12')
  item.update({'full_time_quad':mp.nstr(qv,28),'full_time_quad_relative_error':mp.nstr(rel,8)})
 rows.append(item)
 print('fixed point',j+1,'G',item['G_display'],'floatrel',item['float_generic_relative_error'],flush=True)
xdeg=(1,0,0,0);rrdeg=g.confluent_residues(xdeg);codeg=g.green_allnodes_pi3_coefficient(xdeg)
assert codeg==Q(3,32) and not any(b==0 for b,k in rrdeg)
exdeg=g._mpq(codeg)/mp.pi**3
vdeg=mp.quad(lambda t:g.density_mp_from_confluent(t,rrdeg),[0,mp.mpf(1)/12,mp.mpf(1)/8,1,mp.inf])
assert abs(vdeg-exdeg)<mp.mpf('1e-45')
degenerate={'x':['1','0','0','0'],'G_pi3_exact':str(codeg),'G_display':mp.nstr(exdeg,28),
            'full_confluent_time_quad':mp.nstr(vdeg,28),'relative_error':mp.nstr(abs(vdeg-exdeg)/exdeg,8),
            'zero_poles_remaining':0,'nonzero_residues':[{'b':str(b),'order':k,'R':str(v)} for (b,k),v in sorted(rrdeg.items())]}
assert abs(g.dd_cutoff_log((0,0,1,-1),Q(1,8))-g.dd_cutoff_log((0,0,1,-1),mp.mpf(1)/8))<mp.mpf('1e-50')
geom={'base_rays':R,'body_rates':rates,'base_cones':base,'zero_hyperplane_ranks':rankzero,
      'full_cone_count':len(g.CONES),'unique_ray_count':len(g.RAYS),'global_spectral_vertex_count':40,
      'spectral_unit_volume':str(specvol),'p_t0_pi4_coefficient':str(p0coeff),
      'all_geometry_sha256':hashlib.sha256(json.dumps({'rays':g.RAYS,'cones':g.CONES},sort_keys=True).encode()).hexdigest()}
data={'batch':'U11','status':'fixed geometry and numerical interface diagnostics PASS; no cost or new bound',
      'geometry':geom,'scalar_time_integrals':scalar,'fixed_points':rows,'degenerate_symbolic_point':degenerate,'new_upper_bound':False,
      'run_notes':'Initial same12 run; epsilon Fraction coercion repair same12 rerun; allnode exact implementation final same12 equivalence + one preregistered degenerate e1, no extra generic points or costs',
      'elapsed_seconds':time.time()-start,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'interface_sha256':hashlib.sha256(Path(g.__file__).read_bytes()).hexdigest(),
      'scope':'exact Fraction geometry/G coefficients, mp numerical comparisons not interval roundoff certificates'}
Path(__file__).with_name('u11_results.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
print('PASS',data['elapsed_seconds'],flush=True)
