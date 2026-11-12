from itertools import combinations, product
# At each source point, enumerate exactly the subsets giving a dictionary hit.
for n in range(1,9):
    subsets=[set(s) for k in range(1,n+1) for s in combinations(range(n),k)]
    for bits in product([0,1],repeat=n):
        col={i for i,bit in enumerate(bits) if bit}
        for i in col:
            assert sum(i in s and s<=col for s in subsets)==2**(len(col)-1)
print('dictionary multiplicities exact through eight labels')
