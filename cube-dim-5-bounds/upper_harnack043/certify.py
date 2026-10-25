#!/usr/bin/env python3
"""One H_(1/64) center value+four gradients; global log-paraboloid tests."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from itertools import product
from math import lcm
from pathlib import Path
import gzip,hashlib,json,os,resource,sys,tempfile,time
from flint import arb,fmpq,ctx
import flint
import sympy as sp

if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rq(x):
    x=F(str(x));return fmpq(x.numerator,x.denominator)
def ab(x):return arb(rq(x))
def end(x,up=False):return F(str((x.upper() if up else x.lower()).fmpq()))
def atomic(path,data):
    temp=path.with_suffix(path.suffix+'.tmp')
    with temp.open('w') as f:
        json.dump(data,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(temp,path)
def show(x):return x.str(80)
def floor_dyadic(x,bits=24):return F((x.numerator*(1<<bits))//x.denominator,1<<bits)


def first_jet(ids,vertices,k,D,A,B):
    ids=sorted(ids,key=lambda i:k[i]);nodes=[k[i] for i in ids]
    assert len(set(nodes))==5
    slopes=[tuple(D*z for z in vertices[i][1:]) for i in ids]
    w=[B*B*n*n+A*A for n in nodes];table=[]
    for right,node in enumerate(nodes):
        jet=[[fmpq(0) for _ in range(6)] for _ in range(5)]
        jet[0][right+1]=fmpq(1,2)
        for j in range(4):jet[j+1][0]=fmpq(slopes[right][j]*B*B*node,w[right])
        row=[jet]
        for order in range(1,right+1):
            left=right-order;d=node-nodes[left];dl=[x-y for x,y in zip(slopes[right],slopes[left])]
            P=row[-1];Q=table[right-1][order-1];out=[]
            for ai in range(5):
                coef=[]
                for symbol in range(6):
                    value=P[ai][symbol]-Q[ai][symbol]
                    if ai:value-=dl[ai-1]*out[0][symbol]
                    coef.append(value/d)
                out.append(coef)
            row.append(out)
        table.append(row)
    return table[-1][-1],w


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream-root',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args();began=time.monotonic();root=Path(__file__).resolve().parent
    regpath=root/'registration.json';reg=json.loads(regpath.read_text());ctx.prec=reg['precision_bits']
    for rel,digest in reg['read_only_inputs'].items():assert sha(args.upstream_root/rel)==digest
    source=sha(Path(__file__));regsha=sha(regpath)
    fingerprint=hashlib.sha256(json.dumps({'source':source,'reg':regsha,'inputs':reg['read_only_inputs']},sort_keys=True).encode()).hexdigest()
    work=args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix='upper-harnack043-'))
    work.mkdir(parents=True,exist_ok=True)
    assert not any((p/'.git').exists() for p in (work,*work.parents)),'cache must be outside Git'
    checkpoint=work/'checkpoint.json'
    if checkpoint.exists() and not args.resume:raise ValueError('existing checkpoint requires --resume')
    data=json.loads(gzip.decompress((args.upstream_root/'astra_vertex017/boundary_triangulation.json.gz').read_bytes()))
    assert data['denominator']==144 and len(data['vertices'])==280 and len(data['simplices'])==3504
    vertices=data['vertices'];simplices=data['simplices'];point=tuple(map(F,reg['center']))
    D=lcm(*(p.denominator for p in point));n=[int(D*p) for p in point]
    k=[sum(x*z for x,z in zip(n,v)) for v in vertices];assert len(set(k))==280
    epsilon=F(reg['epsilon']);scaled=144*D*epsilon;A,B=scaled.numerator,scaled.denominator
    stencil=[{} for _ in range(5)];start=0;segments=[]
    if args.resume:
        old=json.loads(checkpoint.read_text());assert old['fingerprint']==fingerprint
        start=old['next_simplex'];segments=old['segments']
        stencil=[{int(key):rq(value) for key,value in row.items()} for row in old['stencil']]
    def save(done,status):
        atomic(checkpoint,{'fingerprint':fingerprint,'status':status,'next_simplex':done,
          'stencil':[{str(key):str(value) for key,value in row.items()} for row in stencil],
          'segments':segments+[{'pid':os.getpid(),'seconds':time.monotonic()-began,'end_simplex':done}]})
    for ci in range(start,3504):
        simplex=simplices[ci];jet,keys=first_jet(simplex['vertices'],vertices,k,D,A,B);det=simplex['det_numerator']
        for ai,row in enumerate(jet):
            for symbol,value in enumerate(row):
                if value==0:continue
                key=0 if symbol==0 else keys[symbol-1]
                combined=stencil[ai].get(key,fmpq(0))-det*value
                if combined==0:stencil[ai].pop(key,None)
                else:stencil[ai][key]=combined
        done=ci+1
        if done%reg['checkpoint_stride_simplices']==0 or done==3504:
            save(done,'complete' if done==3504 else 'checkpoint')
            print(f'completed {done}/3504; {time.monotonic()-began:.3f}s',flush=True)
        if time.monotonic()-began>=reg['checkpoint_seconds'] and done<3504:
            save(done,'checkpoint');raise SystemExit(75)
    keys={key for row in stencil for key in row if key};logs={key:arb(key).log() for key in keys}
    factor=arb(D)**4/(144*(2*arb.pi())**5);quantities=[]
    for row in stencil:
        value=arb(0)
        for key,c in row.items():value+=arb(c)*(logs[key] if key else 1)
        quantities.append(factor*value)
    H=quantities[0];gradH=quantities[1:];assert H>0
    saved=json.loads((args.upstream_root/'astra_green018/batch_result.json').read_text())['selected_value']
    assert tuple(map(F,saved['point']))==point and H.upper()<arb(saved['G_ball']).upper()
    g=[v/H for v in gradH]
    N=sp.Matrix.hstack(*(sp.Matrix(v) for v in reg['fixed_body_basis_columns']))
    assert abs(N.det())==48
    Ninv=N.inv();Q5=Ninv.T*Ninv;Qsp=Q5[1:5,1:5]
    assert Qsp==sp.eye(4)/4-sp.ones(4)/36
    Q=[[F(str(Qsp[i,j])) for j in range(4)] for i in range(4)]
    def quadratic(s):return sum((s[i]*Q[i][j]*s[j] for i in range(4) for j in range(4)),F())
    corners=[];h0=F(reg['fixed_closed_box_tail_halfwidth'])
    for sign in product((-1,1),repeat=4):
        b=sum((v*s for v,s in zip(g,sign)),arb(0));quad=quadratic(sign)
        corners.append({'sign':sign,'b_lower':end(b),'b_upper':end(b,True),'Q':quad})
    assert max(c['Q'] for c in corners)==1
    Hlo,Hhi=end(H),end(H,True)
    rlo=min(c['b_lower']*h0-c['Q']*h0*h0/(epsilon*epsilon) for c in corners)
    rhi=max(c['b_upper']*h0+F(21,8)*c['Q']*h0*h0/(epsilon*epsilon) for c in corners)
    box_lower=ab(Hlo)*ab(rlo).exp();box_upper=ab(Hhi)*ab(rhi).exp()
    target=F(reg['primary_target']);Delta_ball=(ab(Hlo)/ab(target)).log();Delta=end(Delta_ball)
    region=None
    if Delta>0:
        root_rows=[]
        for c in corners:
            b=ab(c['b_lower']);curv=c['Q']/(epsilon*epsilon)
            positive_root=(b+(b*b+4*ab(curv)*ab(Delta)).sqrt())/(2*ab(curv))
            root_rows.append({'sign':list(c['sign']),'Q_sign':str(c['Q']),
                'log_gradient_dot_lower':str(c['b_lower']),
                'positive_root_ball':show(positive_root),'lower':end(positive_root),'upper':end(positive_root,True)})
        hstar_lo=min(c['lower'] for c in root_rows);hstar_hi=min(c['upper'] for c in root_rows)
        hsafe=floor_dyadic(hstar_lo)
        assert hsafe>0
        cube_residuals=[Delta+c['b_lower']*hsafe-c['Q']*hsafe*hsafe/(epsilon*epsilon) for c in corners]
        assert min(cube_residuals)>0
        face_margin=min(min(p,1-p) for p in point[1:]);hface=min(hsafe,face_margin)
        # Q^-1=4I+(4/5)J; use a directed upper norm, not midpoint gradients.
        B2=4*sum((v*v for v in g),arb(0))+ab(F(4,5))*sum(g,arb(0))**2
        B2up=end(B2,True);assert B2up>0
        Bup=end(ab(B2up).sqrt(),True)
        rho_ball=2*ab(Delta)/((ab(Bup)**2+4*ab(Delta)/(ab(epsilon)**2)).sqrt()+ab(Bup))
        rho=floor_dyadic(end(rho_ball));assert rho>0
        ellipse_residual=Delta-Bup*rho-rho*rho/(epsilon*epsilon);assert ellipse_residual>0
        region={'target':str(target),'Delta_lower':str(Delta),
          'maximal_safe_h_for_conservative_corner_test_interval':[str(hstar_lo),str(hstar_hi)],
          'maximal_h_root_rows':[{key:(str(value) if isinstance(value,F) else value) for key,value in row.items()} for row in root_rows],
          'safe_centered_cube_halfwidth':str(hsafe),'safe_halfwidth_clipped_to_unit_face':str(hface),
          'minimum_cube_log_residual':str(min(cube_residuals)),
          'centered_ellipsoid':{'definition':'v^T Q v <= rho^2, x=(1,x0_tail+v)',
            'Q':[[str(v) for v in row] for row in Q],'rho':str(rho),
            'dual_loggradient_norm_upper':str(Bup),'log_residual_lower':str(ellipse_residual)},
          'warning':'Maximum is for this fixed conservative log-paraboloid corner test, not the actual H superlevel set. The ellipsoid may be intersected with the unit face.'}
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result={'status':'center_H_gradient_and_global_logparaboloids_certified',
      'model':'alpha1','epsilon':str(epsilon),'point':reg['center'],'precision_bits':ctx.prec,
      'global_g_or_C5_claim':False,'new_H_objects':1,'all_simplices_processed':3504,
      'H_center_ball':show(H),'H_center_directed_interval':[str(Hlo),str(Hhi)],
      'four_tail_gradient_H_balls':list(map(show,gradH)),'four_tail_loggradient_balls':list(map(show,g)),
      'tail_loggradient_directed_intervals':[[str(end(v)),str(end(v,True))] for v in g],
      'N_columns':reg['fixed_body_basis_columns'],'N_determinant':str(N.det()),
      'N_inverse':[[str(v) for v in row] for row in Ninv.tolist()],
      'Q_face':[[str(v) for v in row] for row in Q],
      'theory_global_Hessian_bounds':'-2||N^-1u||2^2/epsilon^2 <= u^T Hess(log H)u <=(21/4)||N^-1u||2^2/epsilon^2',
      'theory_global_logparaboloid_remainders':['-||N^-1u||2^2/epsilon^2','+(21/8)||N^-1u||2^2/epsilon^2'],
      'fixed_closed_box':{'tail_halfwidth':str(h0),'H_uniform_lower_bound_ball':show(box_lower),
          'H_uniform_lower_bound_dyadic':str(end(box_lower)),
          'H_uniform_upper_bound_ball':show(box_upper),'primary_target_passed':bool(box_lower>ab(target)),
          'secondary_comparison_passed':bool(box_lower>ab(F(reg['secondary_comparison_only']))),
          'all16_corners_checked':True},
      'paraboloid_regions':region,'exact_rational_stencil_entries':sum(len(row) for row in stencil),
      'checkpoint':{'path':str(checkpoint),'sha256':sha(checkpoint),'bytes':checkpoint.stat().st_size,
         'fingerprint':fingerprint,'segments':segments+[{'pid':os.getpid(),'seconds':time.monotonic()-began,'end_simplex':3504}]},
      'runtime':{'PID':os.getpid(),'wall_seconds_this_process':time.monotonic()-began,
         'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,'flint_version':flint.__version__,
         'sympy_version':sp.__version__,'basis_optimizations':0,'epsilon_optimizations':0,'mesh_cells':0},
      'source_sha256':{'registration.json':regsha,'certify.py':source,**reg['read_only_inputs']}}
    args.output_dir.mkdir(parents=True,exist_ok=True);out=args.output_dir/'results.json';atomic(out,result)
    print(json.dumps({'output':str(out),'sha256':sha(out),'H_center':show(H),
      'box_lower':show(box_lower),'box_passed':result['fixed_closed_box']['primary_target_passed'],
      'safe_h':None if region is None else region['safe_centered_cube_halfwidth'],
      'ellipsoid_rho':None if region is None else region['centered_ellipsoid']['rho'],
      'runtime':result['runtime']},indent=2))


if __name__=='__main__':main()
