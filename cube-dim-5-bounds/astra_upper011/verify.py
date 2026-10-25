"""Three rational Young shapes, with positive series and rigorous external tails."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, factorial, floor
import ast
import hashlib
import json
import sys
import time
from flint import arb, arb_series, ctx

ROOT=Path(__file__).resolve().parent
BETA=F(3,4)
sys.set_int_max_str_digits(100000)


def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator


def rational_endpoint(q,upper=False):
    return F(str((q.upper() if upper else q.lower()).fmpq()))


def rising(q,n):
    out=F(1)
    for j in range(n):out*=q+j
    return out


def atomic_json(path,obj):
    tmp=path.with_suffix(path.suffix+'.partial');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(path)


def encode_ball(q):
    text=str(q)
    assert arb(text).contains(q)
    return text


def profile_terms(ones,twos):
    from itertools import product
    degree=sum(ones)+2*sum(twos);out={}
    for ell in product(*(range(k+1) for k in twos)):
        js=tuple(a+b+e for a,b,e in zip(ones,twos,ell));j=sum(js)
        coeff=F(factorial(j),16**j)/rising(BETA,degree)
        for k,e in zip(twos,ell):coeff*=comb(k,e)
        coeff*=BETA**(sum(ones)+sum(twos)+sum(ell))*(1-BETA)**(sum(twos)-sum(ell))
        exponents=tuple(3*a-4*(b+2*c) for a,b,c in zip(js,ones,twos))
        assert all(e<=0 for e in exponents) and coeff>0
        key=(j+1,exponents[0],exponents[1]);out[key]=out.get(key,F(0))+coeff
    return out


def make_terms(geometry):
    families={'d':{},'d2':{},'b':{},'one':{(1,0,0):F(1)}}
    def add(name,ones,twos,weight):
        for key,q in profile_terms(ones,twos).items():families[name][key]=families[name].get(key,F(0))+weight*q
    for key,weight in geometry['basis_group_coefficients'].items():add('d',ast.literal_eval(key),(0,0,0),weight)
    for key,weight in geometry['pair_pattern_weights'].items():
        ones,twos=ast.literal_eval(key);add('d2',ones,twos,weight)
    for p in range(5):add('b',(0,p,4-p),(0,0,0),5*comb(4,p)*4**p*12**(4-p))
    return {name:sorted(terms.items()) for name,terms in families.items()}


def generate_moments(N,terms,fingerprint,old,rawW):
    started=time.monotonic();ctx.cap=N+1
    checkpoint=Path('/tmp')/f'astra_upper011_series_{N}.json'
    if checkpoint.exists():
        state=json.loads(checkpoint.read_text());checksum=state.pop('payload_sha256')
        assert hashlib.sha256(json.dumps(state,sort_keys=True).encode()).hexdigest()==checksum
        assert state['fingerprint']==fingerprint and state['N']==N
    else:state={'fingerprint':fingerprint,'N':N,'completed':{},'family':None,'next_term':0,'partial':None}
    inv=((arb_series([1,-1])**ar(BETA)+5*arb_series([1,-ar(F(1,3))])**ar(BETA)+10)/16).inv()
    denominator_cache={};linear_cache={}
    degrees={'d':5,'d2':10,'b':4,'one':0}
    def save():
        payload=dict(state);payload['payload_sha256']=hashlib.sha256(json.dumps(state,sort_keys=True).encode()).hexdigest()
        atomic_json(checkpoint,payload)
    for name in ['d','d2','b','one']:
        if name in state['completed']:continue
        if state['family']==name:
            total=arb_series([arb(q) for q in state['partial']]);start=state['next_term']
        else:total=arb_series([0]);start=0
        for index in range(start,len(terms[name])):
            (j,e1,e2),coefficient=terms[name][index]
            if j not in denominator_cache:denominator_cache[j]=inv**j
            for e,w in [(e1,F(1)),(e2,F(1,3))]:
                if (e,w) not in linear_cache:linear_cache[e,w]=arb_series([1,-ar(w)])**ar(F(e,4))
            total+=ar(coefficient)*denominator_cache[j]*linear_cache[e1,F(1)]*linear_cache[e2,F(1,3)]
            if time.monotonic()-started>=120:
                state.update({'family':name,'next_term':index+1,'partial':[encode_ball(total[n]) for n in range(N+1)]});save()
                print(json.dumps({'checkpoint':'120 second batch reached','N':N,'family':name,'next_term':index+1}),flush=True)
                return None
        nu=BETA+degrees[name];factor=arb(1);moments=[]
        for n in range(N+1):
            if n:factor*=arb(n)/ar(nu+n-1)
            q=total[n]*factor;assert q>0;moments.append(q)
        W=rawW if name=='one' else list(map(F,old['moments_to_order_22'][{'d':'E_d_W_power','d2':'E_d_squared_W_power','b':'E_b_corner_W_power'}[name]]))
        for n in range(12):
            expected=sum((F(comb(n,k))*2**(n-k)*W[k] for k in range(n+1)),F(0))/12**n
            assert moments[n].contains(ar(expected))
        state['completed'][name]=[encode_ball(q) for q in moments];state.update({'family':None,'next_term':0,'partial':None});save()
        print(json.dumps({'completed_moment_family':name,'N':N,'positive_kernel_terms':len(terms[name]),'seconds':time.monotonic()-started}),flush=True)
    return {name:[arb(q) for q in state['completed'][name]] for name in degrees}


def expectation_bounds(moments,p,theta,N):
    rho=12*theta/(1+2*theta);r=rho*F(p+N+1,N+2);assert 0<rho<1 and r<1
    total=arb(0);weight=arb(1)
    for n in range(N+1):
        total+=weight*moments[n]
        weight*=ar(rho)*arb(p+n)/(n+1)
    total/=ar((1+2*theta)**p)
    lower=rational_endpoint(total);partial_upper=rational_endpoint(total,True)
    qN=rational_endpoint(moments[N],True)
    tail=qN*comb(p+N,N+1)*rho**(N+1)/((1-r)*(1+2*theta)**p)
    assert lower>0 and tail>=0
    return lower,partial_upper+tail,{'qN_upper':str(qN),'rho':str(rho),'ratio_bound':str(r),'tail_upper':str(tail),'partial_lower':str(lower),'partial_upper':str(partial_upper)}


def main():
    if len(sys.argv)!=2 or int(sys.argv[1]) not in [256,512,1024]:raise SystemExit('Pass one registered N: 256, 512, or 1024.')
    N=int(sys.argv[1]);began=time.monotonic();ctx.prec=256
    geometry_path=ROOT.parent/'astra_upper004'/'corner_checkpoint.json';geometry=json.loads(geometry_path.read_text())['geometry']
    old_path=ROOT.parent/'astra_upper005'/'checkpoint.json';old=json.loads(old_path.read_text())
    rawW_path=ROOT.parent/'astra_upper008'/'verification.json';rawW=list(map(F,json.loads(rawW_path.read_text())['raw_W_moments']))
    sources=[ROOT/'registration.md',Path(__file__),geometry_path,old_path,rawW_path]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources};fingerprint=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()
    terms=make_terms(geometry);moments=generate_moments(N,terms,fingerprint,old,rawW)
    if moments is None:return
    AL,BU,cL,cU=(F(old['constants'][key]) for key in ['A_lower','B_upper','c_lower','c_upper'])
    rows=[]
    for theta in [F(1,32),F(1,16),F(3,32)]:
        bounds={};tails={}
        for name,p in [('d',7),('d2',7),('one',7),('b',11)]:
            lower,upper,evidence=expectation_bounds(moments[name],p,theta,N);bounds[name]=(lower,upper);tails[name]=evidence
        Ad,A2,B,V=(bounds[name] for name in ['d','d2','one','b'])
        Ulo,Uhi=AL*Ad[0]-BU*A2[1],AL*Ad[1]-BU*A2[0]
        assert Uhi>0
        row={'theta':str(theta),'N':N,'expectation_bounds':{k:[str(a),str(b)] for k,(a,b) in bounds.items()},'tail_evidence':tails,'U_bounds':[str(Ulo),str(Uhi)]}
        if Ulo>0:
            a=F(floor(float((ar(Ulo)/ar(V[1]))**(arb(1)/4))*1000),1000)
            assert a**4<=Ulo/V[1]<(a+F(1,1000))**4
            lower=cL*(F(11,4)*Ulo*a**7-F(7,4)*V[1]*a**11)
            row.update({'a_selected':str(a),'complete_corner_g_lower':str(lower),'corner_g_lower_display':str(ar(lower)),
                        'passes_existing_g_screen':lower>F(2821,50000000),'passes_128e_6_corner_screen':lower>F(16,125000)})
        family4=cU**4*Uhi**11/V[0]**7
        ideal4=cU**4*Ad[1]**5*B[1]**6/V[0]**7
        ideal4lo=cL**4*Ad[0]**5*B[0]**6/V[1]**7
        row.update({'all_a_old_ell_ceiling_fourth_upper_factorization':'c_upper^4 * U_upper^11 / V_lower^7','old_ell_ceiling_display':str(ar(family4)**(arb(1)/4)),
                    'all_a_old_ell_below_existing_g':family4<F(2821,50000000)**4,
                    'ideal_Jensen_ceiling_fourth_upper_factorization':'c_upper^4 * A_upper^5 * B_upper^6 / V_lower^7','ideal_Jensen_ceiling_display':str(ar(ideal4)**(arb(1)/4)),
                    'ideal_Jensen_ceiling_below_128e_6':ideal4<F(16,125000)**4,
                    'ideal_Jensen_ceiling_below_1e_4':ideal4<F(1,10000)**4,
                    'ideal_Jensen_ceiling_strictly_above_128e_6':ideal4lo>F(16,125000)**4})
        row['classification_stable']=(row.get('passes_existing_g_screen',False) or row['all_a_old_ell_below_existing_g']) and (row['ideal_Jensen_ceiling_below_128e_6'] or row['ideal_Jensen_ceiling_strictly_above_128e_6'])
        rows.append(row)
    out={'status':'three_rational_shapes_necessary_corner_with_rigorous_tail','N':N,'new_global_g':False,'new_C5_upper_bound':False,
         'all_48_prefix_moments_crosschecked':True,'c_lower':str(cL),'c_upper':str(cU),'positive_Stieltjes_kernel_term_counts':{k:len(v) for k,v in terms.items()},
         'candidates':rows,'source_sha256':hashes,'checkpoint_fingerprint':fingerprint,
         'runtime':{'python':sys.executable,'arb_bits':ctx.prec,'seconds':time.monotonic()-began}}
    atomic_json(ROOT/f'verification_{N}.json',out)
    print(json.dumps({'N':N,'candidates':[{k:v for k,v in row.items() if k in ['theta','a_selected','corner_g_lower_display','passes_existing_g_screen','passes_128e_6_corner_screen','old_ell_ceiling_display','all_a_old_ell_below_existing_g','ideal_Jensen_ceiling_display','ideal_Jensen_ceiling_below_128e_6','ideal_Jensen_ceiling_strictly_above_128e_6','classification_stable']} for row in rows],'runtime':out['runtime']},indent=2))


if __name__=='__main__':main()
