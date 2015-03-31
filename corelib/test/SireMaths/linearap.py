
from Sire.Base import *
from Sire.Maths import *
from Sire.Qt import *

costs = NMatrix(8, 8)

rand = RanGenerator()

print("Generating the costs...")

t = QTime()
t.start()

for i in range(0,costs.nRows()):
    for j in range(0,costs.nColumns()):
        d = rand.rand(0,5)
        costs.set(i,j, d*d)

ms = t.elapsed()

print("...took %d ms" % ms)

print("\nFinding the combination with the lowest total cost...")

t.start()

rows_to_columns = solve_linear_assignment(costs)

ms = t.elapsed()

t.start()

solve_linear_assignment(costs, True)

ms2 = t.elapsed()

print("\nSolution:")
print(rows_to_columns)

print("\nSolution took %d ms (%d ms with checking)" % (ms, ms2))

print("Total cost = %f" % calculate_total_cost(costs, rows_to_columns))

t.start()
rows_to_columns2 = brute_force_linear_assignment(costs)
ms = t.elapsed()

print("\nBrute force solution:")
print(rows_to_columns2)

print("Total cost = %f" % calculate_total_cost(costs, rows_to_columns2))

print("\nSolution took %d ms" % ms)

assert( rows_to_columns == rows_to_columns2 )
