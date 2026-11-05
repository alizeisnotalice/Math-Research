from fractions import Fraction as F
import json
cases=[{'c':F(1),'C':F(1),'cy':F(1),'Cy':F(1),'v':F(1,10),'v0':F(1,10),'kap':F(1,10),'z':F(1)},
       {'c':F(1),'C':F(2),'cy':F(1),'Cy':F(3),'v':F(1,20),'v0':F(1,10),'kap':F(1,10),'z':F(2)},
       {'c':F(2),'C':F(3),'cy':F(2),'Cy':F(5),'v':F(1,25),'v0':F(1,5),'kap':F(3,25),'z':F(5)}]
rows=[]
for d in cases:
 c,C,cy,Cy,v,v0,kap,z=(d[k] for k in ('c','C','cy','Cy','v','v0','kap','z'))
 lo=c*cy*v;x=kap*z;up=C*Cy*v0
 margin=x/(1+x)**2-lo/(1+up)**2
 rows.append({'input':{k:str(v) for k,v in d.items()},'assumptions':0<v<=v0 and c*v<=kap<=C*v and cy<=z<=Cy,'margin':str(margin),'pass':margin>=0})
print(json.dumps(rows,indent=2))
