#!/usr/bin/env python3
"""Prepared64-root zero-H time-spread consumer. NOT YET RUN.

73 new rational coefficients,4672 cached mass products, integer chord gains.
No H/gradient/R/log/exp/power-root/old-baseline function is imported.
"""
from datetime import datetime,timezone
from fractions import Fraction as F
from pathlib import Path
import argparse,gzip,hashlib,json,os,time

HERE=Path(__file__).resolve().parent
N=1<<160;D192=1<<192
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p);return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def put(p,x,atomic=False):
    p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp') if atomic else p
    with tmp.open('w' if atomic else 'x') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    if atomic:os.replace(tmp,p)
def ceilq(x):return -((-x.numerator)//x.denominator)
def theta(a,b):
    za,zb=1/a,F(0) if b is None else 1/b
    mean=F(4,5)*(za**5-zb**5)/(za**4-zb**4)
    t=(mean-zb)/(za-zb);assert 0<t<1
    return t
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reviewed-plan-sha',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--resume',action='store_true');args=ap.parse_args()
    assert sha(HERE/'plan.json')==args.reviewed_plan_sha;plan=read(HERE/'plan.json')
    assert sha(__file__)==plan['consumer_sha256'] and sha(HERE/'selected64.json')==plan['selected64_sha256']
    pre=read(HERE/'preflight_results.json');assert sha(HERE/'preflight_results.json')==plan['preflight_results_sha256']
    for src in pre['source_locks'].values():assert sha(src['path'])==src['sha256']
    for src in pre['selected_packet_locks'].values():assert sha(src['path'])==src['sha256']
    assert pre['selected']==64 and pre['geometry_ids']==['g002'] and pre['band_count']==73
    identity={'plan_sha256':args.reviewed_plan_sha,'selected64_sha256':plan['selected64_sha256'],'consumer_sha256':sha(__file__)}
    start=time.monotonic();sel=read(HERE/'selected64.json');g='g002';times=pre['time_endpoints']
    src=pre['source_locks'];geos=read(src['geometries.json.gz']['path']);old=read(src['factor_catalog.json']['path']);cat=read(src['ratio_catalog.json']['path'])
    tree=read(src['numeric_hierarchy.json.gz']['path']);prof=geos[g]['profile'];h,mu=F(prof['h']),F(prof['mu_infinity'])
    before=F(pre['saved_complete_UH'])*N;assert before.denominator==1;before=int(before)
    if args.resume:
        st=read(args.output/'state.json');assert st['identity']==identity and st['status']!='complete'
        assert 0<=len(st['factors'])<=73 and 0<=len(st['nodes'])<=64
        assert [x['a'] for x in st['factors']]==times[:len(st['factors'])]
        assert [x['key'] for x in st['nodes']]==[x['key'] for x in sel[:len(st['nodes'])]]
        assert st['gain_units']==sum(x['gain_units'] for x in st['nodes'])
        st['invocations'].append(os.getpid())
    else:
        args.output.mkdir(exist_ok=False)
        put(args.output/'registration.json',{'task':'upper_integral150_one64root_time_spread_batch','pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),**identity,
            'roots':64,'factors':73,'cached_band_products':4672,'workers':1,'pool':None,'wall_cap':None,'source_locks':src,'packet_locks':pre['selected_packet_locks'],
            'new_H_R_gradient_log_exp_powerroot_calls':0,'old_baseline_recomputations':0})
        st={'identity':identity,'status':'factors','factors':[],'nodes':[],'gain_units':0,'invocations':[os.getpid()]};put(args.output/'state.json',st,True)
    def endpoints(t):
        if t is None:return F(1),F(1)
        a=F(t);R=F(cat[str(h/a)]['upper']);T=1+mu/a+F(25,24)*F(old[g][t]['J_upper'])/a**2
        return R,T
    try:
        for j in range(len(st['factors']),73):
            a=F(times[j]);bt=times[j+1] if j+1<73 else None;b=None if bt is None else F(bt);w=theta(a,b)
            ra,ta=endpoints(times[j]);rb,tb=endpoints(bt);rb,tb=min(rb,ra),min(tb,ta)
            range_time=w*ra+(1-w)*rb;taylor_time=w*ta+(1-w)*tb;oldC=F(old[g][times[j]]['factor_upper'])
            exact=min(oldC,range_time,taylor_time);rounded=F(ceilq(exact*D192),D192)
            assert 1<=rounded<=oldC
            st['factors'].append({'a':times[j],'b':bt,'theta':str(w),'range_time_upper':str(range_time),'taylor_time_upper':str(taylor_time),'old_upper':str(oldC),'new_upper':str(rounded)})
            put(args.output/'state.json',st,True)
        st['status']='roots';put(args.output/'state.json',st,True)
        for i in range(len(st['nodes']),64):
            item=sel[i];row=item['source_node'];cert=item['source112_certificate'];key=item['key'];n=tree['nodes'][key]
            assert row['level']==5 and row['geometry']==g and key in tree['roots'] and n['parent'] is None and not n['children']
            assert n['effU']==n['selfU']==item['saved_selfU']==cert['new_selfU']<cert['old_selfU']
            assert row['old_U_i']>=1<<144 and row['old_L_i']<row['old_U_i']
            assert sha(row['packet_path'])==row['packet_sha256'];raw=read(row['packet_path']);ends=raw['endpoints']
            assert raw['point']==row['point'] and [e['t'] for e in ends]==times
            total=F(0)
            for j,e in enumerate(ends):
                nxt=ends[j+1]['H'] if j+1<73 else {'lower':'0','upper':'0'}
                ml=F(e['H']['lower'])-F(nxt['upper']);up=F(e['H']['upper'])-F(nxt['lower']);assert 0<ml<=up
                total+=up*F(st['factors'][j]['new_upper'])
            oldmi=cert['new_M_i'];newmi=min(oldmi,ceilq(total*N));assert row['old_L_i']<=newmi<=oldmi<=row['old_U_i']
            gain=row['weight']*(oldmi-newmi)//(16*(1<<(4*row['level'])))
            upper=item['saved_selfU']-gain;assert item['saved_selfL']<=upper<=item['saved_selfU']
            st['gain_units']+=gain;st['nodes'].append({'rank':i,'key':key,'old_mean_i':oldmi,'new_mean_i':newmi,'safe_secant_slope_lower':'1/16',
                'old_selfU':item['saved_selfU'],'new_selfU':upper,'gain_units':gain,'packet_sha256':row['packet_sha256'],'prefix_complete_UH_units':before-st['gain_units']})
            put(args.output/'state.json',st,True)
        for obj in list(src.values())+list(pre['selected_packet_locks'].values()):assert sha(obj['path'])==obj['sha256']
        result={'task':'upper_integral150','status':'finite64root_candidate_complete_pending_independent_primary','pid':os.getpid(),'factors':73,'roots':64,'cached_band_products':4672,
            'old_UH':str(F(before,N)),'new_UH':str(F(before-st['gain_units'],N)),'exact_UH_gain':str(F(st['gain_units'],N)),
            'strict_global_UH_improvement':st['gain_units']>0,'mean_improved_roots':sum(x['new_mean_i']<x['old_mean_i'] for x in st['nodes']),
            'no_powerroot_replay':True,'new_H_R_gradient_log_exp_powerroot_calls':0,'old_baseline_recomputations':0,'global_Cn_claim':False,'seconds':time.monotonic()-start,
            'registration_sha256':sha(args.output/'registration.json'),'all_sources_unchanged':True}
        put(args.output/'factors.json',st['factors']);put(args.output/'nodes.json',st['nodes']);put(args.output/'results.json',result)
        st['status']='complete';put(args.output/'state.json',st,True);print(json.dumps(result,sort_keys=True))
    except BaseException as ex:
        put(args.output/f"fault{len(st['invocations']):04d}.json",{'pid':os.getpid(),'exception':repr(ex),'saved_factors':len(st['factors']),'saved_roots':len(st['nodes']),'prefix_preserved':True,'action':'stop,no automatic retry/new H/refinement'})
        raise
if __name__=='__main__':main()
