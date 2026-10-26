#!/usr/bin/env python3
"""Prepared, NOT run by the105 author.

Verify saved finite rational/root objects only; no root search, old budget,
kernel, source program or external-theorem proof is executed.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import json

H=Path(__file__).resolve().parent
load=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
n=0
def ck(ok,message):
    global n
    n+=1
    if not ok:raise AssertionError(message)
def pair(obj):return Q(obj['lower']),Q(obj['upper'])

m=load(H/'manifest.json');r=load(H/'results.json');inp=load(H/'coefficient_inputs.json')
for ent in m['artifacts']:
    ck(sha(H/ent['path'])==ent['sha256'],'artifact hash '+ent['path'])
sources=load(H/'source_manifest.json')
for ent in sources['local_source_files']:
    ck(sha(ent['path'])==ent['sha256'],'read-only source hash '+ent['path'])
ck(r['status']=='PASS_105_BACKUP_TRANSFER_CONDITIONAL108_DOMINANCE_NO_FULL104_RECALC','saved status')
ck(len(r['groups'])==3 and r['assertions']==91,'saved group/check ledger')
cap=r['cap_certificate'];delta=Q(cap['delta'])
ck(delta==Q(1,8) and delta**2/4==Q(cap['single_cap_fraction_upper'])==Q(1,256),'cap fraction')
ck(60*Q(1,256)==Q(cap['union_fraction_upper'])==Q(15,64),'union fraction')
eta=Q(15,64)+Q(49,64)*Q(114,115)**3
ck(eta==Q(cap['eta'])==Q(95408781,97336000) and Q(31,32)<eta<1,'eta')
ck(248**2>2*175**2 and Q(20)-14*Q(248,175)==Q(4,25),'cstar residuals')
ck(2*114**2>160**2 and 93**2>8*32**2,'relative gain/conditional comparison residuals')
expected={'sharp_norm_factor':(2,Q(8,27)),
          'sharp_I_factor':(5,Q(8,27)**3),
          '105_cap_I_factor':(5,eta**6*Q(8,27)**3),
          'sqrt3_for_conditional_BL_factor':(2,Q(3))}
roots={}
for cert in r['root_certificates']:
    label=cert['label'];degree=cert['degree'];rad=Q(cert['radicand']);a,b=pair(cert)
    ck(label not in roots and expected[label]==(degree,rad),'saved root binding '+label)
    ck(0<a<=b and a**degree<=rad<=b**degree,'saved integer-power root certificate '+label)
    roots[label]=(a,b)
ck(set(roots)==set(expected),'complete saved root list')
s=r['old104_formal_budget_scaling'];rk=roots['sharp_norm_factor'];ra=(eta*rk[0],eta*rk[1])
ck(pair(s['sharp_norm_factor'])==rk and pair(s['105_cap_norm_factor'])==ra,'norm factors')
ck(pair(s['sharp_I_factor'])==roots['sharp_I_factor'] and pair(s['105_cap_I_factor'])==roots['105_cap_I_factor'],'I factors')
ck(ra[1]<rk[0]<rk[1]<1,'strict105 norm improvement')
pl,pu=map(Q,inp['pi_interval']);a,b=roots['sqrt3_for_conditional_BL_factor']
ck(3<pl<pu,'inherited pi enclosure')
rb=(8/(3*b*pu),8/(3*a*pl));bl=r['conditional108_comparison']
ck(pair(bl['BL_over_old098_factor'])==rb,'conditional BL factor')
ratio=(rb[0]/ra[1],rb[1]/ra[0])
ck(pair(bl['BL_over_105_cap_factor'])==ratio and ratio[1]<1,'conditional strict dominance')
ck(bl['premise_accepted_or_proved_here'] is False and bl['full104_budget_recalculation_omitted'] is True,'conditional108 scope')
ck(s['actual_G_or_I_lower_bound'] is False and s['numeric_complete104_budget_recalculated'] is False,'budget scope')
ck(r['no_new_Cn'] is True and r['new_Cn_claim'] is False,'no Cn claim')
print(json.dumps({'status':'PASS_MINIMAL_SAVED_FINITE_OBJECTS','checks':n,
    'new_root_searches':0,'old104_100_or_kernel_executions':0,
    'scope':'Does not prove external Xie/classical input,104 acceptance,108 premise,globalg or fullI.'}))
