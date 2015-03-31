
from Sire.Squire import *
from Sire.Maths import *
from Sire.Units import *

S = S_GTO(0.5, 1)
P = P_GTO(0.5, 1)
C = PointCharge( Vector(0,0,0), 1*mod_electron )
D = PointCharge( Vector(0,0,1.4), 1*mod_electron )

x = Vector(0,0,0)
y = Vector(0.0,0.0,1.4)

ss = SS_GTOs( [S, S], [x, y] )

print("SS overlap")
print(ss.overlap_integral())

print("SS kinetic")
print(ss.kinetic_integral())

print("SS potential")
print(ss.potential_integral( [C, D] ))

print("SS potential 3")
print(ss.potential_integral( [C, D], 3 ))

ps = PS_GTOs( [P], [x], [S], [y] )

print("PS overlap")
print(ps.overlap_integral())

print("PS kinetic")
print(ps.kinetic_integral())

print("PS potential")
print(ps.potential_integral( [C, D] ))

print("PS potential 3")
print(ps.potential_integral( [C, D], 3 ))

pp = PP_GTOs( [P, P], [x, y] )

print("PP overlap")
print(pp.overlap_integral())

print("PP kinetic")
print(pp.kinetic_integral())

print("PP potential")
print(pp.potential_integral( [C, D] ))

print("PP potential 3")
print(pp.potential_integral( [C, D], 3 ))
