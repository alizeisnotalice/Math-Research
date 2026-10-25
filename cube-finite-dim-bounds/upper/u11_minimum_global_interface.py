"""Locked 3-case whole-G rational interface. No g/I/cost computations on import.
G*pi^3=N17/D20; homogeneous pair a^3*N17,D20 is for max-coordinate face.
Boundary D=0 is a removable/true singular set; use exact allnode fallback,
never 0. Whole-cell coefficient ratios may be used after an affine substitution.
"""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,importlib.util
from flint import fmpq_mpoly_ctx,fmpq
DATA=Path(__file__).with_name('u11_minimum_global_results.json')
DATA_SHA='43be4de041c9da823d192ff6d8ef2440c3668f88b13fc52dbff320dd0011d96f'
assert hashlib.sha256(DATA.read_bytes()).hexdigest()==DATA_SHA
_raw=json.loads(DATA.read_text());CTX=fmpq_mpoly_ctx.get(('a','b','c','d'),ordering='lex');A=CTX.gens()[0]
def _poly(dic):return CTX.from_dict({tuple(map(int,e.split(','))):fmpq(v) for e,v in dic.items()})
PAIRS={r['case']:(_poly(r['numerator_dict']),_poly(r['denominator_dict'])) for r in _raw['physical_cases']}
CONE_RAYS={'I':((1,0,0,0),(1,1,0,0),(2,1,1,0),(3,1,1,1)),
           'II':((2,1,1,0),(1,1,0,0),(1,1,1,0),(1,1,1,1)),
           'III':((2,1,1,0),(1,1,0,0),(1,1,1,1),(3,1,1,1))}
def homogeneous_pair(case):
 n,d=PAIRS[case]
 return A**3*n,d

def coefficient_ratio_bounds(P,D):
 """Exact coefficientwise rational bound, provided every coefficient >=0.
 None upper means unmatched positive numerator support, NOT a zero integral.
 """
 pp=P.to_dict();dd=D.to_dict()
 if any(v<0 for v in pp.values()) or any(v<0 for v in dd.values()):
  return {'status':'negative polynomial coefficient, no ratio conclusion'}
 if not dd:return {'status':'zero denominator polynomial'}
 ratios=[Q(str(pp.get(e,0)))/Q(str(v)) for e,v in dd.items() if v>0]
 extras=[e for e,v in pp.items() if v>0 and not dd.get(e,0)]
 return {'status':'coefficientwise bound','lower':min(ratios),'upper':None if extras else max(ratios),'unmatched_numerator_terms':len(extras)}

def coefficient_exact_sorted(x):
 x=tuple(map(Q,x));a,b,c,d=x
 if not(a>=b>=c>=d>=0) or a==0:raise ValueError('need sorted nonnegative nonzero point')
 case='I' if a>=b+c+d else 'II' if b+c>=a+d else 'III'
 N,D=PAIRS[case];v=[fmpq(z.numerator,z.denominator) for z in x];den=D(*v)
 if den:return Q(str(N(*v)/den))
 src=Path(__file__).with_name('u11_green.py');sp=importlib.util.spec_from_file_location('minimum_allnode',src);g=importlib.util.module_from_spec(sp);sp.loader.exec_module(g)
 return g.green_allnodes_pi3_coefficient(x)
