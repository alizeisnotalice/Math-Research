#!/usr/bin/env python3
"""Exact closed-event maxima and relevance pruning of saved affine fan cells."""
from pathlib import Path
from fractions import Fraction as Q
import sys,os,time,json
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path[:0]=[str(ROOT/'lower'),str(P)]
from certify_rectangular import field
from profile_variation import ENV,load,save,sha,rt,enc

def main():
    start=time.monotonic();reg=load(P/'finalization_registration.json')
    assert sha(Path(__file__))==reg['source_sha256']
    for x in reg['inputs']:assert sha(Path(x['path']))==x['sha256']
    F,u,v,bounds,_,_=field();zero=F(0)
    def dec(d):
        x=[Q(0)]*25
        for ij,c in d.items():i,j=map(int,ij.split(','));x[5*i+j]=Q(c)
        return F(a=x)
    ev=load(P/'events.json')['events'];events=[(dec(x['numerator']),dec(x['positive_denominator'])) for x in ev]
    for n,d in events:assert zero<d
    for (n,d),(nn,dd) in zip(events,events[1:]):assert n*dd<nn*d
    assert events[0][0]==events[0][1]*Q(-1,3) and events[-1][0]==events[-1][1]*Q(1,3)
    rows=[];AB=[]
    for i in range(6):
        r=load(P/f'cell{i}_derivative.json');ck=load(P/f'cell{i}_checkpoint.json')
        assert r['union_engines_agree_constant_derivative_and_faces'] and ck['status']=='COMPLETE_two_union_engines'
        t=Q(r['sample_t']);assert events[i][0]<events[i][1]*t and events[i+1][1]*t<events[i+1][0]
        A=dec(r['fan_cell_linear_coefficients']['constant']);B=dec(r['fan_cell_linear_coefficients']['slope'])
        assert A+B*t==dec(r['directional_derivative_coefficients']);AB.append((A,B))
        rows.append({'cell':i,'left_event':i,'right_event':i+1,'constant':enc(A),'slope':enc(B),'sample_t':str(t),'result_sha256':sha(P/f'cell{i}_derivative.json'),'runtime':r['runtime']})
    def ratio_bounds(n,d):
        nl,nh=bounds(n.a);dl,dh=bounds(d.a);assert dl>0
        v=[nl/dl,nl/dh,nh/dl,nh/dh];return min(v),max(v)
    endpoint_values=[];true_events=[];maximum=None
    for i,(n,d) in enumerate(events):
        A,B=AB[min(i,5)];value=A*d+B*n
        if 0<i<6:
            AA,BB=AB[i-1];assert AA*d+BB*n==value
            if (AA,BB)!=(A,B):true_events.append(i)
        lo,hi=ratio_bounds(value,d);assert hi<0
        if maximum is None or maximum[0]*d<value*maximum[1]:maximum=(value,d,i)
        endpoint_values.append({'event':i,'value_numerator':enc(value),'positive_denominator':enc(d),'value_interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)],'strict_negative':True})
    n,d=events[0];oldplus=dec(load(ROOT/'lower_round023/plus_derivative.json')['directional_derivative_coefficients'])
    assert AB[0][0]*d+AB[0][1]*n==oldplus*d/3
    n,d=events[-1];oldminus=dec(load(ROOT/'lower_round023/minus_derivative.json')['directional_derivative_coefficients'])
    assert AB[-1][0]*d+AB[-1][1]*n==oldminus*d/3
    mid=dec(load(P/'middle_derivative.json')['directional_derivative_coefficients'])
    n,d=events[3];assert n==zero and AB[2][0]==mid/2 and AB[3][0]==mid/2
    A=mid/2;B=AB[1][1];E=oldplus/4
    compact=None
    if AB==[(E,-E),(A,B),(A,B),(A,-B),(A,-B),(E,E)]:
        an=E-A;ad=-B-E
        assert zero<an and zero<ad and zero<B and E<zero
        assert an*events[5][1]==events[5][0]*ad
        compact={'formula':'D(t)=A-B*abs(t) forabs(t)<=alpha;D(t)=E*(1+abs(t)) forabs(t)>=alpha','A':enc(A),'B':enc(B),'E':enc(E),'alpha_positive_numerator':enc(an),'alpha_positive_denominator':enc(ad),'A_interval':list(map(str,bounds(A.a))),'B_interval':list(map(str,bounds(B.a))),'E_interval':list(map(str,bounds(E.a))),'matches_all_six_cell_coefficients_exactly':True,'maximum_location':'t=0'}
    merged=[]
    for row in rows:
        if merged and row['constant']==merged[-1]['constant'] and row['slope']==merged[-1]['slope']:merged[-1]['right_event']=row['right_event']
        else:merged.append({k:row[k] for k in ('left_event','right_event','constant','slope')})
    lo,hi=ratio_bounds(maximum[0],maximum[1]);assert hi<0
    # Old files are a strict read-only dependency, including the previous rounds.
    for oldname in ('lower_round021','lower_round022','lower_round023'):
        old=ROOT/oldname
        for x in load(old/'archive_whitelist.json')['files']:assert sha(old/x['path'])==x['sha256']
    out={'status':'COMPLETE_STRICT_NEGATIVE_entire_period3_zero_mean_directional_fan','normalization':'range=max(z)-min(z)=1,sum(z)=0,z0>=z1>=z2','sorted_direction':'((1-t)/2,t,-(1+t)/2),-1/3<=t<=1/3','S3_mechanism':'phasecycles j->j+1 andreflection j->-j onfull3a1period generateS3;safe6144rawlistcompletebyexteriornonintersection;doesnotrequireliteralpermutationofstartlist','zero_direction_derivative':'0','scaling_rule':'D(lambda*z)=lambda*D(z) forlambda>=0','candidate_event_count':7,'true_breakpoint_event_indices':true_events,'compact_formula':compact,'merged_affine_pieces':merged,'all_six_cells':rows,'closed_event_values':endpoint_values,'normalized_fan_maximum':{'event_index':maximum[2],'value_numerator':enc(maximum[0]),'positive_denominator':enc(maximum[1]),'interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)],'strict_negative':True},'closed_continuity_verified_exactly':True,'old_endpoints_and_middle_verified_after_correct_range_scaling':True,'all021022023_archive_files_unchanged':True,'new_union_jobs_in_finalizer':0,'runtime':rt(start),'source_sha256':reg['source_sha256'],'registration_sha256':sha(P/'finalization_registration.json'),'scope':'strict negative firstvariation for every nonzero firstaxisperiod3 zero mean mass direction ofthisfinite1024lengthfunctional;noCn orgeneralweightedinputupper,noquantifiedepsilon,mainownsproofandL1/collar'}
    save(P/'fan_results.json',out);print(json.dumps({k:v for k,v in out.items() if k not in ('all_six_cells','closed_event_values','merged_affine_pieces')},indent=2))
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
