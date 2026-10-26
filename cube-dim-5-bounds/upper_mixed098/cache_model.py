"""Epsilon-aware cache keys. A same-point H64 entry is never an H128 hit."""
from fractions import Fraction as F
def epsilon_key(epsilon):return str(F(epsilon))
def point_key(point):return ','.join(map(str,map(F,point)))
def pair_key(epsilon,point):return (epsilon_key(epsilon),point_key(point))
def lookup(state,epsilon,point):
 e,p=pair_key(epsilon,point)
 return state['cache_by_epsilon'].get(e,{}).get(p)
def insert(state,epsilon,record):
 e,p=pair_key(epsilon,record['point'])
 assert str(record['epsilon'])==e,'explicit record epsilon must agree with bucket'
 bucket=state['cache_by_epsilon'].setdefault(e,{})
 assert p not in bucket,'completed records are immutable; precision resumes use own checkpoint'
 bucket[p]=record
