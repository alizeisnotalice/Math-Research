"""Portable exact verification of the new algebraic C5 upper and comparisons."""
from fractions import Fraction as Q
from pathlib import Path
import json
P=Path(__file__).resolve().parent.parent
def read(n):return json.loads((P/n).read_bytes())
a=read('round016/primary_upper139_acceptance.json')
new=read('upper_increment139/actual/results.json')
ind=read('audit_increment140/exact_actual/results.json')
old=read('round015/primary_upper_acceptance.json')
g=Q(a['g']);I=Q(a['complete_I_upper'])
assert I==Q(new['complete_I_upper'])==Q(ind['complete_I_upper'])<Q(old['complete_I_upper'])
v=3125*I**4/g**5
assert v==Q(a['best_upper']['fourth_power'])<Q(old['best_upper']['fourth_power'])<32**4
l,u=map(Q,a['best_upper']['interval']);assert 0<l<u and l**4<=v<=u**4
assert u<=Q(a['simple_rational_upper'])<Q(old['best_upper']['interval'][0])
print('PASS exact new full upper and old/new comparison; no root generation or H')
