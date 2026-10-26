"""Typed new-record insertion; immutable source cache is never mutated."""
from fractions import Fraction as F
from binding import KERNEL_SHA,MODEL
def pair_key(epsilon,point):return (str(F(epsilon)),','.join(map(str,map(F,point))))
def validate_record(record,target):
 assert str(F(record['epsilon']))=='1/128'
 assert F(record['qualification_target'])==F(target)==F(11,100000)
 assert record['pair_cache_key']==list(pair_key(record['epsilon'],record['point']))
 assert record['kernel_sha256']==KERNEL_SHA and record['model']==MODEL
 assert len(record['values'])==5
def insert_new(new_cache,source_cache_by_epsilon,record,target):
 validate_record(record,target);e,p=pair_key(record['epsilon'],record['point'])
 assert p not in source_cache_by_epsilon[e] and p not in new_cache,'completed typed cache entries cannot be overwritten'
 new_cache[p]=record
