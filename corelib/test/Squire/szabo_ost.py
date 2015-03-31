
from Sire.Squire import *
from Sire.Maths import *
from Sire.Units import *

zeta1 = 2.0925
zeta2 = 1.2400

coef = [0.444635, 0.535328, 0.154329 ]
expn = [0.109818, 0.405771, 2.227660 ] 

p0 = Vector(0,0,0)
p1 = Vector(1.4632,0,0)

A1 = [0,0,0]
A2 = [0,0,0]

for i in range(0,3):
    A1[i] = expn[i] * zeta1**2
    A2[i] = expn[i] * zeta2**2

#STO-3G
s0 = CS_GTO( A1, coef )
s1 = CS_GTO( A2, coef )

print(s0.alpha()) 
print(s0.scale())
print(s1.alpha())
print(s1.scale())
print(Vector.distance(p0,p1), Vector.distance2(p0,p1))

n0 = PointCharge( Vector(0,0,0), 2*mod_electron )
n1 = PointCharge( Vector(1.4632,0,0), 1*mod_electron )

print(overlap_integral( CSS_GTO(p0,s0, p0,s0) ))
print(overlap_integral( CSS_GTO(p0,s0, p1,s1) ))
print(overlap_integral( CSS_GTO(p1,s1, p1,s1) ))

hf = HF()

hf.add(p0, s0)
hf.add(p1, s1)
hf.add(n0.center(), n0.charge()*mod_electron)
hf.add(n1.center(), n1.charge()*mod_electron)

hf.solve()
