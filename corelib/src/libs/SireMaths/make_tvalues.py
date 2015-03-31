
lines = open("t_values", "r").readlines()

words = lines[0].split()

import string

def join(vals, typ):
    numbers = []

    if typ == "float":
        for val in vals:
            numbers.append( str( float(val.replace('%','')) ) )
    else:
        for val in vals:
            numbers.append( str( int(val.replace('%','')) ) )

    return string.join(numbers, ", ")

nlevels = len(words[1:])
print "    int nlevels = %d;" % nlevels
print "    double levels[%d] = { %s };" % (nlevels, join(words[1:], "float"))

counts = []

for line in lines[1:]:
    words = line.split()
    counts.append(words[0])

print "    int ncounts = %d;" % len(counts)
print "    int counts[%d] = { %s };" % (len(counts), join(counts, "int"))

vals = []

for line in lines[1:]:
    vals.append( line.split()[1:] )

print "    double values[%d][%d] = {" % (len(counts), nlevels)

a = []

for val in vals:
    a.append("        { %s }" % join(val, "float"))

print "%s\n    };" % string.join(a, ",\n")

