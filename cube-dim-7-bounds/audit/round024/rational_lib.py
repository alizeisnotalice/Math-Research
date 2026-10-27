"""Exact positive scalar intervals: rational log/atanh tail and integer square roots.
Atanh terms j<N have tail <2*x^(2N+1)/((2N+1)*(1-x*x)).
For sqrt(q), floor sqrt(floor(q*S*S)) / S is a rigorous lower bound;
add 1/S for an upper bound unless exact. No floating point or external library.
"""
from fractions import Fraction as Q
import math,sys
sys.set_int_max_str_digits(0)
N=128
S=10**160
def at(x):
 assert 0<=x<1
 v=2*sum((x**(2*j+1)/Q(2*j+1)for j in range(N)),Q())
 return v,v+2*x**(2*N+1)/(Q(2*N+1)*(1-x*x))
def log(q):
 q=Q(q);assert q>1;k=q.numerator.bit_length()-q.denominator.bit_length();r=q/Q(2)**k
 if r<1:k-=1;r*=2
 assert 1<=r<2 and k>=0
 l,u=at((r-1)/(r+1));a,b=at(Q(1,3));return k*a+l,k*b+u
def sqrt(q):
 q=Q(q);assert q>=0;m=math.isqrt(q.numerator*S*S//q.denominator)
 return Q(m,S),Q(m if m*m*q.denominator==q.numerator*S*S else m+1,S)
def sumroot(A,Bl,Bu):
 return A+Bl+2*sqrt(A*Bl)[0],A+Bu+2*sqrt(A*Bu)[1]
def F24(a):
 a=Q(a);l,u=log(a);A=(1+2/a)**7;b=(1+1/a)**7
 return sumroot(A,b*(1+7*l),b*(1+7*u))
def E24(a):
 a=Q(a);l,u=log(a);assert a-7*u>0;c=((a+1)/(a+2))**5
 return c*(a-7*u)**2-4*(1+7*u),c*(a-7*l)**2-4*(1+7*l)
def F25(t):
 t=Q(t);l,u=log(t);return sumroot((1+1/t)**7,1+7*l,1+7*u)
def G26(t):
 t=Q(t);l,u=log(t);c=t**7/(t+1)**5-1;return c-7*u,c-7*l
def grid(upper):return Q((upper.numerator*10**12+upper.denominator-1)//upper.denominator,10**12)
def encode(v):return [str(v[0]),str(v[1])]
