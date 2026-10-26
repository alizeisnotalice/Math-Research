#!/usr/bin/env python3
"""Prepared portable checker. No absolute paths; not run in this session.

--decode-only only checks exact byte projection/schema. Arithmetic mode uses
accepted endpoint intervals,091 source theorem and saved qualified factors;
it does not certify those analytic sources by trusting a PASS string.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,ceil
from itertools import combinations_with_replacement,product
from collections import Counter
from math import factorial,prod
import argparse,gzip,hashlib,json,sys
import codec as C
N=1<<160
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def lines(p):
    with gzip.open(p,'rt') as f:
        for s in f:yield json.loads(s)
def fourth_upper(q):
    d=1<<224;n=q.numerator*d**4//q.denominator;v=isqrt(isqrt(n));assert v**4<=n<(v+1)**4
    return F(v if F(v,d)**4==q else v+1,d)
def qualify_factors(p,times,geos,factors):
    # Reader-only source arithmetic. Never invoked during167 preparation/export.
    sys.path.insert(0,str(p/'math'))
    import pathmean as P,secant as Q
    cat=read(p/'ratios112.json');old=read(p/'factors112.json')
    for g,obj in geos.items():
        prof=P.weight_profile(obj['offsets'])
        for n in ('h','mu_infinity','J0','M2'):assert prof[n]==F(obj['profile'][n])
        for j,t in enumerate(times):
            e=Q.secant_envelope(P.darboux_plan(obj['offsets'],F(t),8),cat)
            assert e['J_upper']==F(old[g][t]['J_upper'])
            a,b=1/F(t),F(0) if j==72 else 1/F(times[j+1])
            w=(F(4,5)*(a**5-b**5)/(a**4-b**4)-b)/(a-b)
            def end(z,key):
                if key is None:return F(1),F(1)
                return F(cat[str(prof['h']*z)]['upper']),1+prof['mu_infinity']*z+F(25,24)*z*z*F(old[g][key]['J_upper'])
            ra,ta=end(a,t);rb,tb=end(b,None if j==72 else times[j+1]);raw=min(F(old[g][t]['factor_upper']),w*ra+(1-w)*min(ra,rb),w*ta+(1-w)*min(ta,tb))
            assert F(factors[g][j]['new_upper'])==F(ceil(raw*(1<<192)),1<<192)
    # Certify stored091 interval uppers by a separately scaled integer enclosure.
    for xs,r in cat.items():
        x=F(xs);roots=[]
        for den in (16,256):
            q=1+x*x/den;D=1<<512;n=isqrt(q.numerator*D*D//q.denominator);roots.append(F(n+(n*n*q.denominator<q.numerator*D*D),D))
        assert (roots[0]+x/4)**2*(roots[1]+x/16)**8<=F(r['upper'])
def chord(li,ui,mi,level,weight):
    l,u,m=F(li,N),F(ui,N),F(mi,N);assert 0<l<=m<=u
    a,b=fourth_upper(l**5),fourth_upper(u**5)
    return weight*ceil((a if l==u else ((u-m)*a+(m-l)*b)/(u-l))*F(N,1<<(4*level)))
def arithmetic(row,vals,times,catalog,factors):
    masses=[];lower_grid=0;upper=F(0);mean=F(0)
    for j,(l,u) in enumerate(vals):
        nl,nu=vals[j+1] if j<72 else (F(0),F(0));ml,mh=l-nu,u-nl;assert 0<ml<=mh
        R=F(catalog[str(F(row['h_infinity'])/F(times[j]))]['upper']);assert R>=1
        term=ml/R;lower_grid+=(term.numerator<<320)//term.denominator
        masses.append(term);upper+=mh*R;mean+=mh*F(factors[row['geometry']][j]['new_upper'])
    li,ui=row['L_i'],row['U_i'];assert li*(1<<160)<=lower_grid or F(li,N)<=sum(masses,F(0))
    assert ui==3*N or upper<=F(ui,N)<3
    assert row['new_mean_i']==ceil(mean*N)<row['m112_i'] and li<=row['new_mean_i']<=ui
    assert chord(li,ui,row['m112_i'],row['level'],row['weight'])<=row['self112U']
    assert chord(li,ui,row['new_mean_i'],row['level'],row['weight'])<=row['new_selfU']<=row['self112U']
def main():
    ap=argparse.ArgumentParser();ap.add_argument('packet',type=Path);ap.add_argument('--decode-only',action='store_true');a=ap.parse_args();p=a.packet
    done=read(p/'completion.json')
    for n,h in done['files'].items():assert sha(p/n)==h
    info=read(p/'certificate.json');times=info['times'];groups=[[] for _ in range(15)]
    seen=set()
    for row in lines(p/'node_index.jsonl.gz'):
        assert row['key'] not in seen;seen.add(row['key']);groups[row['shard']].append(row)
    pm={r['id']:r for r in lines(p/'projection_map.jsonl.gz')};assert len(pm)==97971
    if not a.decode_only:
        cats={g:read(p/f'range_{g}.json') for g in ('092','095','103')};fac=read(p/'mean_factors.json')['geometry_factors'];geos=read(p/'geometries.json.gz')
        qualify_factors(p,times,geos,fac)
        for cat in cats.values():
            for xs,r in cat.items():
                x=F(xs);a0,a1,b0,b1=(F(r[n]) for n in ('sqrt1_lower','sqrt1_upper','sqrt2_lower','sqrt2_upper'))
                assert a0*a0<=1+x*x/16<=a1*a1 and b0*b0<=1+x*x/256<=b1*b1
                assert (a1+x/4)**2*(b1+x/16)**8<=F(r['upper'])
        # Same closed positive receiver face only;never a cap on outside anchor lines.
        assert F(info['M132_cap'])<3 and info['cap_domain']=='receiver_complete_F_H512_only'
    total=0;selfs={}
    for sid,rows in enumerate(groups):
        with gzip.open(p/f'endpoints{sid:02d}.bin.gz','rb') as f:
            assert C.take(f,8)==C.MAGIC and int.from_bytes(C.take(f,4),'big')==sid and int.from_bytes(C.take(f,4),'big')==len(rows)
            for row in rows:
                vals,proj,original=C.decode(f,row,times);assert pm[row['id']]=={'id':row['id'],'source_gzip_sha256':row['packet_sha256'],'original_raw_sha256':original,'normalized_projection_sha256':proj}
                if not a.decode_only:
                    off=[(F(i,1<<row['level'])-F(c),F(i+1,1<<row['level'])-F(c)) for i,c in zip(row['indices'],row['point'][1:])]
                    assert off==[tuple(map(F,s)) for s in geos[row['geometry']]['offsets']]
                    assert row['point'][0]=='1' and row['weight']==factorial(4)//prod(factorial(c) for c in Counter(row['indices']).values())
                    assert all(0<=i<(1<<row['level']) for i in row['indices'])
                    arithmetic(row,vals,times,cats[row['range_source']],fac)
                selfs[row['key']]=row['new_selfU'];total+=1
            assert f.read(1)==b''
    assert total==97971
    if not a.decode_only:
        tree=read(p/'tree.json.gz');eff={};special=info['unused_legacy_self_key'];ns=tree['nodes']
        assert len(tree['roots'])==52360 and len(set(tree['roots']))==52360 and len(ns)==97972
        assert {tuple(ns[k]['indices']) for k in tree['roots']}=={tuple(reversed(i)) for i in combinations_with_replacement(range(32),4)}
        assert set(ns)-set(selfs)=={special}
        for k,n in ns.items():
            assert n['weight']==factorial(4)//prod(factorial(c) for c in Counter(n['indices']).values())
            if n['parent'] is None:assert k in tree['roots'] and n['level']==5
        for k in sorted(tree['nodes'],key=lambda k:-tree['nodes'][k]['level']):
            n=tree['nodes'][k]
            if n['children']:
                expected={tuple(sorted((2*i+b for i,b in zip(n['indices'],bits)),reverse=True)) for bits in product((0,1),repeat=4)}
                assert expected=={tuple(ns[c]['indices']) for c in n['children']}
                assert len(n['children'])==len(expected) and sum(ns[c]['weight'] for c in n['children'])==16*n['weight']
                assert all(ns[c]['parent']==k and ns[c]['level']==n['level']+1 for c in n['children'])
                child=sum(eff[c] for c in n['children']);eff[k]=child if k==special else min(selfs[k],child)
            else:assert k!=special;eff[k]=selfs[k]
        uh=F(sum(eff[k] for k in tree['roots']),N);assert uh==F(info['candidate_UH'])
        I=uh+F(info['qualified_J139_upper']);assert I>=uh
    print(json.dumps({'status':'decode_only_pass' if a.decode_only else 'accepted_interval_input_arithmetic_pass','packets':total,'new_H_R_calls':0,
      'numeric_H_interval_truth_shared_accepted_dependency':True,'source_theorems_require_analytic_review':True,'global_C5_claim':False}))
if __name__=='__main__':main()
