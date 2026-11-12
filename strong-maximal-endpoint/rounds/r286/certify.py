from fractions import Fraction as F
for m in range(1,41):
    total=F(0)
    soft=F(0)
    for j in range(1,m+1):
        d=2**j
        assert d*(F(1,4)+F(1,2))/2==F(3,8)*d
        soft+=F(1,4)+F(1,8*d)
        total+=F(3,8)*d
        assert F(5,8)*d-F(1,2)>=F(3,8)*d
        assert total<=F(3,4)*d
    assert soft/(3*m)>=F(1,12)
for j in range(21):
    near=0 if j==0 else 2**(j-1)
    far=2**j
    assert F(1+far,1+near)<=2
print('40 moving-peak configurations and 21 expansion factors passed')
