
from Sire.Test import *

r = RunTests()

r.run()

a = TestA()

assert(a.__class__ == TestA)

b = TestB()

b.printMe()

a.printMe()

print(b.__class__)
print(b.__class__.__base__)

assert(b.__class__ == TestB)

a.printMe()
b.printMe()

c = castToA(b)

print(c.__class__)
print(c.__class__.__base__)

assert(c.__class__ == TestB)

c.printMe()

c = createB()

print(c.__class__)
print(c.__class__.__base__)

assert(c.__class__ == TestB)

c.printMe()

c = createC()

print(c.__class__)
print(c.__class__.__base__)

assert(c.__class__ == TestC)

c.printMe()

d = createD()
print(d.__class__)
print(d.__class__.__base__)

assert(d.__class__== TestD)

d.printMe()

a.functionA()
b.functionA()
c.functionA()
d.functionA()

b.functionB()
d.functionB()

c.functionC()
d.functionD()

a = createType("A", 1)
b = createType("B", 2)
c = createType("C", 3)
d = createType("D", 4)
e = createType("E", 5)
d2 = createType("D", 6)

print(a.__class__)
print(b.__class__)
print(c.__class__)
print(d.__class__)
print(e.__class__)
print(d2.__class__)

assert(a.__class__ == TestA)
assert(b.__class__ == TestB)
assert(c.__class__ == TestC)
assert(d.__class__ == TestD)
assert(e is None)
assert(d2.__class__ == TestD)

a.functionA()
b.functionB()
c.functionC()
d.functionD()
print(e is None)
d2.functionD()

