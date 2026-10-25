"""Seven registered scalar configurations; no frame/tuple/kernel enumeration."""
import hashlib,json,time,os
from fractions import Fraction as Q
from pathlib import Path
from flint import arb,fmpq,ctx
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def B(q):
 q=Q(q);return arb(fmpq(q.numerator,q.denominator))
def ends(v):
 l,e=v.lower().man_exp();u,f=v.upper().man_exp();e,f=int(e),int(f);g=min(e,f)
 return dict(lower_mantissa=str(int(l)*2**(e-g)),upper_mantissa=str(int(u)*2**(f-g)),common_binary_exponent=str(g))
def main():
 if (P/'results.json').exists():raise SystemExit('Refuse repeat finite diagnostic')
 start=time.monotonic();ctx.prec=512;reg=json.loads((P/'registration.json').read_text())
 for f,h in reg['source_sha256'].items():assert sha(Path(f))==h
 source=next(Path(f) for f in reg['source_sha256'] if f.endswith('exact_projection_peaks/results.json'))
 data=json.loads(source.read_text());R={v['id']:Q(v['R_exact']) for v in data['cases']}
 S=[10*R['r3_4110']+20*R['r3_3210'],2*R['r3_4110']+8*R['r3_3210']+6*R['r3_2220']+8*R['r3_3111']+6*R['r3_2211'],12*R['r3_3210']+18*R['r3_2211']]
 M=[];configs=[];pi=arb.pi()
 for k,s in enumerate(S,1):
  aa=64*R[f'r2_k{k}'];cc=4096*s;assert cc>960*aa
  A=B(aa)/pi**6;C=B(cc)/pi**6;z=(A/(C-960*A)).sqrt()
  m=C/B(960).sqrt()*(B(960).sqrt()*z).atan()+A/z
  M.append(m);configs.append(dict(configuration=k,kind='pair orbit coefficient',k=k,S_exact=str(s),A_times_pi6=str(aa),C_times_pi6=str(cc),C_minus960A_times_pi6=str(cc-960*aa),M=ends(m),third_only_over_min=(str((C*pi/(2*B(960).sqrt()))/m))))
 MT=192*M[0]+480*M[1]+320*M[2]
 for u in [Q(1,4),Q(1,2),Q(3,4)]:
  r2=3*(1+u*u);dv2=Q(3,2)*(1-u)**2;dw2=Q(3,2)*(1+u)**2
  lf=Q(4,9)/(1-u*u)**2;assert Q(1,r2)*(1/dv2+1/dw2)==lf
  L=B(lf)*M[2];collapsed=MT/(B(r2)*B(dv2))
  configs.append(dict(configuration=len(configs)+1,kind='representative middle',u=str(u),r_squared=str(r2),d_v_squared=str(dv2),d_w_squared=str(dw2),leading_coefficient_over_M3_exact=str(lf),pair_resolved_leading_coefficient=ends(L),collapsed_common_distance_leading_coefficient=ends(collapsed),collapsed_over_resolved=ends(collapsed/L)))
 delta=Q(1,64);a=B(Q(3,32)).sqrt()-B(3).sqrt()*B(delta);g=1/(4*B(2).sqrt())-B(3).sqrt()*B(delta)
 assert a.lower()>0 and g.lower()>0
 A=2*M[2]/a**2;BG=(MT-2*M[2])/(a**2*g);q=Q(6,5);v4=pi*pi/2
 # 034: face coarea sqrt(3*3), cone fraction1/3, middle length1/2.
 singular=A*(2*v4/B(4-q))**B(Q(5,6))*B(delta)**B(Q(7,3))
 background=BG*(v4*B(Q(1,2))*B(delta)**4)**B(Q(5,6))
 configs.append(dict(configuration=7,kind='one fixed uniform representative tube',delta=str(delta),line_gap_a=ends(a),other_plane_gap_g=ends(g),tube_A_over_rho=ends(A),background_B=ends(BG),singular_L6over5_upper=ends(singular),background_L6over5_upper=ends(background),total_L6over5_upper=ends(singular+background),relative_to_1_over11250=ends((singular+background)*11250),ordered_pair_count=992,target_ordered_pairs=2,other_ordered_pairs=990))
 assert len(configs)==reg['planned_configurations']==7<=reg['accurate_algebra_Arb_scalar_configuration_cap']
 result=dict(status='CONDITIONAL_ON_033_PENDING_F035',configurations=configs,precision_bits=512,first_level_Q5_coefficient=ends(32*B(R['r1'])/pi**6),global_M_ordered_sum=ends(MT),source_sha256=sha(Path(__file__)),pid=os.getpid(),elapsed_seconds=time.monotonic()-start,kernel_calls=0,old_enumerations=0,grids=0,T_searches=0,MATLAB_calls=0,new_C6_bound=False)
 (P/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],configurations=7,seconds=result['elapsed_seconds'],representative_norm_display=str(singular+background),coefficient_ratios_display=[str(MT/(M[2]*(1+B(Q(1-u,1+u))**2))) for u in [Q(1,4),Q(1,2),Q(3,4)]])))
if __name__=='__main__':main()
