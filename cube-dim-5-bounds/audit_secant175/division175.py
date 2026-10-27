"""No producer import: cross-gcd cancellation and binary long division."""
from math import gcd

def gain(weight,delta,sn,sd,level):
    assert type(weight) is int and weight>0 and type(delta) is int and delta>=0
    assert type(sn) is int and type(sd) is int and sn>0 and sd>0 and 16*sn>=sd
    assert type(level) is int and level>=0
    numerators=[weight,delta,sn];denominators=[sd,1<<(4*level)]
    for i in range(3):
        for j in range(2):
            g=gcd(numerators[i],denominators[j]);numerators[i]//=g;denominators[j]//=g
    A=numerators[0]*numerators[1]*numerators[2];B=denominators[0]*denominators[1]
    assert A>=0 and B>0
    quotient=0;remainder=A
    for shift in range(max(0,A.bit_length()-B.bit_length()),-1,-1):
        threshold=B<<shift
        if remainder>=threshold:remainder-=threshold;quotient|=1<<shift
    assert 0<=remainder<B and quotient*B+remainder==A
    assert quotient*B<=A<(quotient+1)*B
    return quotient,A,B,remainder

def units160(text):
    fields=text.split('/');assert len(fields) in (1,2)
    n=int(fields[0]);d=int(fields[1]) if len(fields)==2 else 1;assert d>0
    n,r=divmod(n<<160,d);assert r==0;return n

def ratio160(n):
    d=1<<160;g=gcd(n,d);n//=g;d//=g
    return str(n) if d==1 else f'{n}/{d}'
