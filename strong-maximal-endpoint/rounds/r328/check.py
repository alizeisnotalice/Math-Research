from fractions import Fraction as F
for N,H in [(2,2),(5,3),(17,2),(32,8)]:
    r=F(1,H)
    v=[1]+[1-r]*(N-2)+[1]
    assert all(sum(r**abs(i-j)*v[j] for j in range(N))==1+r for i in range(N))
    union=N-(N-1)*r
    capacity=sum(v)/(1+r)
    assert 1<=union/capacity<1+r
    print(N,H,union,capacity,union/capacity)
