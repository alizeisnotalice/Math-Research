from fractions import Fraction as F
from pathlib import Path
import json
p=Path(__file__).resolve().parent
r6=json.loads((p/'probe-R6-j3.json').read_text());r8=json.loads((p/'formal-R8-j7.json').read_text())
a=lambda row,key:list(map(F,row[key]))
I1=list(map(F,[0,'5/3','-13/4','31/12','-43/30','11/30','1/15',0,0]))
I2=list(map(F,[0,'5/3','-13/4','33/5','-227/30','287/40','-4079/630','1277/280','-2719/1120']))
checks={}
checks['R6_R8_same_first_layers']=all(a(r6['rows'][j],k)==a(r8['rows'][j],k)[:7] for j in range(4) for k in ('moment_I_coefficients','area_coefficients'))
checks['I1_matches_R067']=a(r8['rows'][1],'moment_I_coefficients')==I1
checks['I2_matches_R067']=a(r8['rows'][2],'moment_I_coefficients')==I2
checks['area_shift_all_layers']=all(a(r8['rows'][j+1],'area_coefficients')==[F(0),F(1)]+[-x for x in a(r8['rows'][j],'moment_I_coefficients')[1:8]] for j in range(7))
valuation=[]
for j in range(6):
 A=a(r8['rows'][j],'area_coefficients'); B=a(r8['rows'][j+2],'area_coefficients');
 val=next((i for i,(x,y) in enumerate(zip(A,B)) if x!=y),None);valuation.append(val)
checks['two_step_valuation_at_least_j_plus_2']=all(x is None or x>=j+2 for j,x in enumerate(valuation))
checks['c3_is_241_over_60']=a(r8['rows'][2],'moment_I_coefficients')[3]-a(r8['rows'][1],'moment_I_coefficients')[3]==F(241,60)
result={'checks':checks,'two_step_area_valuation':valuation,'canonical_area_candidate_coefficients':r8['rows'][7]['area_coefficients'],'all_pass':all(checks.values()),'qualification':'canonical fixed point through degree8 depends on the independent L1 contraction tail estimate'}
print(json.dumps(result,indent=2))
