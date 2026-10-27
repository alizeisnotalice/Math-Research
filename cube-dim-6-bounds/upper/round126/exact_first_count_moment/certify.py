#!/usr/bin/env python3
"""UNRUN. One gated Arb512 run of the eight true first-count moments.
No MC, native, time grid, old peak/fan producer, or125 optimizer is called.
"""
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
import traceback
from collections import Counter
from fractions import Fraction as Q

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
SOURCE=Path(__file__).resolve()

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,data):
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');temporary.replace(path)
def cpu():
    r=resource.getrusage(resource.RUSAGE_SELF);return r.ru_utime+r.ru_stime
def display(iv):return [format(float(q),'.17g') for q in iv]

def read_nodes(spec,old,sp):
    """Read the old source, validate the new spectral nodes, never rebuild a fan."""
    paths=[(HERE/spec['fan_file'],spec['fan_sha256']),
           (PROJECT/spec['audited_cells_relative'],spec['audited_cells_sha256']),
           (PROJECT/spec['independent_audit_relative'],spec['independent_audit_sha256'])]
    for path,pin in paths: assert sha(path)==pin,str(path)
    certificate=json.loads(paths[1][0].read_text());audit=json.loads(paths[2][0].read_text())
    R0=Q(spec['old_full_R0'])
    assert Q(certificate['full_exp_minus_symbol_integral'])==Q(audit['full_exp_minus_symbol_integral'])==R0==Q(old['R0'])
    assert audit['status']=='PASS'
    cells=certificate['cells'];assert len(cells)==21
    tokens=list(map(int,paths[0][0].read_text().split()));assert tokens[0]==21
    cursor=1;V=[(1,)+tail for tail in itertools.product((-1,1),repeat=5)]
    frame=sp.Matrix(V).T;assert frame*frame.T==32*sp.eye(6)
    grouped={};pairs={};base=Q(0);zero_variables=0
    for cid,cell in enumerate(cells):
        determinant=tokens[cursor];cursor+=1;assert 0<determinant<=10**9
        rays=[];rates=[]
        for j in range(6):
            rate=tokens[cursor];ray=tokens[cursor+1:cursor+7];cursor+=7
            assert 0<rate<=10**6 and len(ray)==6 and max(map(abs,ray))<=10**6
            audited=list(map(Q,cell['rays_a0_equals_one'][j]));scale=Q(ray[0])/audited[0]
            assert scale>0 and all(Q(a)==scale*b for a,b in zip(ray,audited))
            assert Q(rate)==scale*Q(cell['rates'][j])
            rays.append(ray);rates.append(rate)
        assert abs(int(sp.Matrix([[rays[j][i] for j in range(6)] for i in range(6)]).det()))==determinant
        scales=[Q(rays[j][0])/Q(cell['rays_a0_equals_one'][j][0]) for j in range(6)]
        assert Q(determinant)==Q(cell['determinant'])*math.prod(scales)
        assert all(rate==sum(abs(sum(v[i]*ray[i] for i in range(6))) for v in V) for rate,ray in zip(rates,rays))
        alpha=Q(determinant,math.prod(rates));base+=alpha
        for vid,v in enumerate(V):
            dots=[sum(v[i]*ray[i] for i in range(6)) for ray in rays]
            assert min(dots)>=0 or max(dots)<=0,('mixed spectral signs',cid,vid,dots)
            nodes=tuple(sorted(Q(abs(dot),rate) for dot,rate in zip(dots,rates) if dot))
            assert 1<=len(nodes)<=6 and all(0<b<Q(1,2) for b in nodes)
            zero_variables+=6-len(nodes)
            grouped[nodes]=grouped.get(nodes,Q(0))+46080*alpha
            pairs.setdefault(nodes,[]).append([cid,vid])
    assert cursor==len(tokens) and 46080*base==R0  # Old normalization consistency only.
    assert sum(len(v) for v in pairs.values())==672 and len(grouped)<=672
    assert sum(grouped.values())==32*R0
    table=[]
    for tid,nodes in enumerate(sorted(grouped)):
        table.append({'id':tid,'nodes':[str(b) for b in nodes],'weight':str(grouped[nodes]),'source_cone_label_pairs':pairs[nodes]})
    save(HERE/'node_table.json',{'scope':'new spectral functional tuples derived from pinned audited source','ungrouped_pairs':672,
         'unique_tuples':len(table),'zero_Exp_variables_removed':zero_variables,'old_R0_normalization_consistency':True,
         'max_multiplicity':max(max(Counter(row['nodes']).values()) for row in table),'tuples':table})
    return table

def rational_jets(sp):
    b,h=sp.symbols('b h');c=(1-b)/b;compiled={};source={}
    def polynomial(expr):
        return [(powers,Q(int(coef.p),int(coef.q))) for powers,coef in sp.Poly(expr,b,h).terms()]
    def rational_form(expr):
        numerator,denominator=sp.fraction(sp.cancel(expr))
        return polynomial(numerator),polynomial(denominator)
    for m in range(1,7):
        A=b**m/(1-2*b);coeff=[A,b**(m-1)*(1-b)/(1-2*b),sp.Integer(0),-A]
        for k in range(m):
            coeff=[sp.cancel(v) for v in coeff]
            assert sp.cancel(coeff[0].subs(h,0)+coeff[3].subs(h,0))==0
            assert sp.cancel(coeff[1]+coeff[3]-sp.diff(b**(m-1),b,k))==0
            compiled[(m,k)]=[rational_form(v) for v in coeff]
            source[str((m,k))]=[sp.sstr(v) for v in coeff]
            A,B,C,D=coeff
            coeff=[sp.diff(A,b)-2*c*h*A/b**2,sp.diff(B,b),sp.diff(C,b)+2*A/b**2,sp.diff(D,b)]
    save(HERE/'rational_jet_coefficients.json',{'initial_and_recurrence':'proof(6)','h0_and_hInfinity_symbolic_checks':True,'coefficients':source})
    return compiled

def eval_rational(compiled,b,h):
    def poly(terms):return sum((coefficient*b**powers[0]*h**powers[1] for powers,coefficient in terms),Q(0))
    numerator,denominator=compiled;den=poly(denominator);assert den!=0
    return poly(numerator)/den

def hermite(nodes,jet,denominator=lambda value:value):
    values=[jet(node,0) for node in nodes]
    for order in range(1,len(nodes)):
        for j in range(len(nodes)-order):
            if nodes[j]==nodes[j+order]:values[j]=jet(nodes[j],order)
            else:values[j]=(values[j+1]-values[j])/denominator(nodes[j+order]-nodes[j])
    return values[0]

def main():
    if (HERE/'state.json').exists():print('Existing state: no overwrite/restart.',file=sys.stderr);return 2
    spec=json.loads((HERE/'inputs.json').read_text());release=json.loads((HERE/'primary_release.json').read_text())
    assert release['authorized'] is True and release['source_sha256']==sha(SOURCE) and release['input_sha256']==sha(HERE/'inputs.json')
    started={'pid':os.getpid(),'unix_time':time.time(),'source_sha256':sha(SOURCE),'input_sha256':sha(HERE/'inputs.json'),'release_sha256':sha(HERE/'primary_release.json')}
    with (HERE/'started.json').open('x') as handle:
        json.dump(started,handle,indent=2);handle.write('\n');handle.flush();os.fsync(handle.fileno())
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
    start=time.monotonic();c0=cpu();current_h=None;current_tuple=None;code=0
    state={'status':'RUNNING','pid':os.getpid(),'completed_h':[],'source_sha256':sha(SOURCE),'input_sha256':sha(HERE/'inputs.json')};save(HERE/'state.json',state)
    try:
        import sympy as sp
        from flint import arb,fmpq,ctx
        ctx.prec=512
        assert spec['precision_bits']==512 and spec['workers']==1
        assert spec['h_values']==[str(Q(1,2**j)) for j in range(12,4,-1)]
        assert sha(HERE/spec['old_caps_file'])==spec['old_caps_sha256']
        old=json.loads((HERE/spec['old_caps_file']).read_text())
        assert Q(old['pi_cubed_B1'])==64*Q(old['R1'])
        for row in old['rows']:
            path=Path(row['source_path']);assert sha(path)==row['source_sha256']
            original=json.loads(path.read_text())
            assert row['old_B1z_expression_interval']==[original['caps']['b_lower'][1],original['caps']['beta_upper'][1]]
        table=read_nodes(spec,old,sp);compiled=rational_jets(sp)
        def A(q):q=Q(q);return arb(fmpq(q.numerator,q.denominator))
        def exact_ends(value):
            out=[]
            for endpoint in (value.lower(),value.upper()):
                man,exp=endpoint.man_exp();out.append(Q(int(man))*Q(2)**int(exp))
            return out
        def grid(value,bits):
            lo,hi=exact_ends(value);scale=1<<bits
            return Q((lo*scale).numerator//(lo*scale).denominator,scale),Q(-((-hi*scale).numerator//(-hi*scale).denominator),scale)
        # The exact polynomial hInfinity check uses the SAME confluent table.
        for row in table:
            nodes=list(map(Q,row['nodes']));m=len(nodes)
            def jet_poly(node,k):
                return Q(math.factorial(m-1),math.factorial(m-1-k)*math.factorial(k))*node**(m-1-k)
            assert hermite(nodes,jet_poly)==1
        pi=arb.pi();sqrtpi=pi.sqrt();rows=[];previous=None
        for index,hs in enumerate(spec['h_values']):
            current_h=hs;h=Q(hs);ah=A(h);sqrt_h=ah.sqrt();exp_h=(-ah).exp();f=sqrt_h.erf();E=exp_h*sqrt_h/sqrtpi
            cache={};Tcache={};functions=[];JL=Q(0);JU=Q(0)
            def jets(m,node):
                key=(m,node)
                if key not in cache:
                    if node not in Tcache:
                        c=1/A(node)-1;T=exp_h*(c*c*ah).hypgeom_u(arb(1)/2,arb(1)/2)/sqrtpi
                        assert T.is_finite() and 0<T<1;Tcache[node]=T
                    values=[]
                    for k in range(m):
                        a,b,c,d=[eval_rational(r,node,h) for r in compiled[(m,k)]]
                        values.append((A(a)*Tcache[node]+A(b)*f+A(c)*E+A(d))/math.factorial(k))
                    cache[key]=values
                return cache[key]
            for row in table:
                current_tuple=row['id'];nodes=list(map(Q,row['nodes']));m=len(nodes)
                value=hermite(nodes,lambda node,k:jets(m,node)[k],A)
                assert value.is_finite()
                lo,hi=grid(value,spec['functional_endpoint_grid_bits'])
                assert 0<lo<=hi<1 and hi-lo<=lo/Q(2**spec['relative_functional_interval_bits_required'])
                weight=Q(row['weight']);JL+=weight*lo;JU+=weight*hi
                functions.append({'tuple_id':row['id'],'interval':[str(lo),str(hi)]})
            assert 0<JL<=JU<32*Q(spec['old_full_R0'])
            if previous is not None:assert previous<JL
            previous=JU
            muL=grid(A(JL)/pi**3,spec['moment_endpoint_grid_bits'])[0]
            muU=grid(A(JU)/pi**3,spec['moment_endpoint_grid_bits'])[1]
            saved=old['rows'][index];assert saved['h']==hs
            oldL,oldU=map(Q,saved['old_B1z_expression_interval'])
            checkL,checkU=grid(A(old['pi_cubed_B1'])*(1-exp_h)/pi**3,384)
            assert oldL<=checkL<=checkU<=oldU
            assert muL<=oldU,'contradiction with the old valid first-count cap'
            result={'index':index,'h':hs,'status':'TRUE_FIRST_COUNT_MOMENT_INTERVAL','pi_cubed_EdN_interval':[str(JL),str(JU)],
                    'EdN_interval':[str(muL),str(muU)],'old_B1z_expression_interval':[str(oldL),str(oldU)],
                    'strict_cap_improvement':muU<oldL,'old_lower_minus_new_upper':str(oldL-muU),
                    'old_cap_over_new_true_interval':[str(oldL/muU),str(oldU/muL)],
                    'EdN_display':display([muL,muU]),'cap_improvement_factor_display':display([oldL/muU,oldU/muL]),
                    'functional_values':functions,'functional_count':len(functions),'precision_bits':512,
                    'node_table_sha256':sha(HERE/'node_table.json'),'coefficient_sha256':sha(HERE/'rational_jet_coefficients.json')}
            save(HERE/('h_%d.json'%index),result);rows.append({k:v for k,v in result.items() if k!='functional_values'})
            state['completed_h'].append(hs);save(HERE/'state.json',state)
            print(json.dumps({'h':hs,'EdN_display':result['EdN_display'],'strict_improvement':result['strict_cap_improvement'],'factor_display':result['cap_improvement_factor_display']}),flush=True)
        state['status']='COMPLETE';save(HERE/'state.json',state)
        result={'status':'COMPLETE_EIGHT_TRUE_AUXILIARY_MOMENT_INTERVALS','rows':rows,'unique_functional_tuples':len(table),
                'ungrouped_pairs_per_h':672,'precision_bits':512,'new_C6_bound':False,'M_or_global_g0_certificate':False,
                'old_R0_used_only_for_source_normalization':True,'MC_samples':0,'old_peak_producer_calls':0,'physical_point_calls':0,
                'dependencies':{'sympy':sp.__version__},'node_table_sha256':sha(HERE/'node_table.json'),'coefficient_sha256':sha(HERE/'rational_jet_coefficients.json')}
    except BaseException as exc:
        code=1;result={'status':'FAILED_NO_RETRY','h':current_h,'tuple_id':current_tuple,'exception':repr(exc),'traceback':traceback.format_exc()}
        state['status']='FAILED';save(HERE/'state.json',state)
    result.update({'source_sha256':sha(SOURCE),'input_sha256':sha(HERE/'inputs.json'),'pid':os.getpid(),
                   'wall_seconds':time.monotonic()-start,'cpu_seconds':cpu()-c0,'exit_code':code,'child_processes':0,
                   'time_or_CPU_ceiling':None,'automatic_precision_increase':False})
    save(HERE/'results.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','traceback')},indent=2),flush=True)
    return code

if __name__=='__main__':sys.exit(main())
