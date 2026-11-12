from itertools import product
q = 7
points = list(product(range(q), repeat=2))
pattern = {(0,1),(2,4),(3,0),(6,3)}
for x,y in points:
    hits=sum(((x-i)%q,(y-j)%q) in pattern for i,j in points)
    assert hits==len(pattern)
print(len(points), 'column cells; common occupation count',len(pattern))
