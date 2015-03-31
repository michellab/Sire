#!/bin/env python
# -*- coding: utf-8 -*-

import sys
import math

lams = None
all_grads_f = []
all_grads_b = []

outdir = sys.argv[1]

its = {}

for arg in sys.argv[2:]:
    smallest = None

    if arg.find("-") == -1:
        its[ int(arg) ] = 1

    else:
        words = arg.split("-")

        start = int(words[0])
        stop = int(words[1])

        if start > stop:
            tmp = start
            start = stop
            stop = tmp

        for i in range(start, stop+1):
            its[i] = 1

iterations = its.keys()
iterations.sort()

print >>sys.stderr,"Averaging over iterations %s" % iterations

def calculatePMF(gradients):
    """This function calculates and return the PMF given the passed series
       of lambda values and gradients"""

    pmf = {}

    lamvals = gradients.keys()
    lamvals.sort()

    if lamvals[0] != 0:
        #we need to start from 0
        gradients[0] = gradients[lamvals[0]]
        lamvals.insert(0, 0)

    if lamvals[-1] != 1:
        #we need to end with 1
        gradients[1] = gradients[lamvals[-1]]
        lamvals.append(1)

    #start at 0
    pmf[ lamvals[0] ] = 0.0

    for i in range(1,len(lamvals)):
        last_lam = lamvals[i-1]
        this_lam = lamvals[i]

        delta_lam = this_lam - last_lam

        pmf[this_lam] = pmf[last_lam] + (delta_lam * 0.5 * (gradients[this_lam] + \
                                                            gradients[last_lam]))

    return pmf

convergence = []

delta_lam = 0.001
T = 310.15  # 37 C in kelvin

for it in iterations:
    filename = "%s/results_%0004d.log" % (outdir,it)
    print >>sys.stderr,"Reading file %s" % filename

    try:
        lines = open(filename, "r").readlines()
    except:
        print >>sys.stderr,"No output for iteration %d" % it
        break

    i = 0
    got_grads = False

    grads_f = {}
    grads_b = {}

    lamvals = []

    while not got_grads:
        line = lines[i]

        i += 1

        if line.find("delta_lambda") != -1:
            words = line.split()
            delta_lambda = float(words[2])
        
        elif line.find("temperature") != -1:
            words = line.split()
            T = float(words[2])

        elif line.find("Lambda") != -1:
            while not got_grads:
                line = lines[i]
                words = line.split()

                if len(words) == 0:
                    got_grads = True
                else:
                    lamvals.append( float(words[0]) )
                    grads_f[ float(words[0]) ] = float(words[3])
                    grads_b[ float(words[0]) ] = float(words[4])
                    i += 1          

    lams = lamvals
    all_grads_f.append(grads_f)
    all_grads_b.append(grads_b)

    pmf_f = calculatePMF(grads_f)
    pmf_b = calculatePMF(grads_b)

    convergence.append( 0.5*(pmf_f[lamvals[-1]] + pmf_b[lamvals[-1]]) )

print "\nConvergence"
for i in range(0,len(iterations)):
    try:
        print "%d  %f" % (iterations[i], convergence[i])
    except:
        pass

# Histogram the convergence so that we can see if there is any skew
def histogram( values, mean = None, stdev = None ):
    minval = values[0]
    maxval = values[0]

    for value in values:
        if minval > value:
            minval = value
        if maxval < value:
            maxval = value

    binwidth = (maxval - minval) / 25

    hist = []

    for i in range(0,26):
        hist.append( [ minval + (i+0.5)*binwidth, 0, 0 ] )

    for value in values:
        bin = int( (value-minval) / binwidth )
        hist[bin][1] += 1

    if not mean is None:
        bin = int( (mean-minval) / binwidth )
        nmean = hist[bin][1]

        for bin in hist:
            bin[2] = nmean * math.exp( -( (bin[0]-mean)**2 / (2*stdev*stdev) ) )

    return hist

def getMeanVarError(values):
    n = 0                                         
    avg = 0
    avg2 = 0

    for value in values:
        n += 1
	avg += value
        avg2 += (value*value)

    avg /= n
    avg2 /= n

    stdev = math.sqrt(avg2 - (avg*avg))
    err = 2 * stdev / math.sqrt(n)

    return (avg, stdev, err)

(mean, stdev, error) = getMeanVarError(convergence)

hist = histogram(convergence, mean, stdev)

print "\nHistogram of free energies"
for bin in hist:
    print "%f  %d  %f" % (bin[0], bin[1], bin[2])


print "\nMean == %f  Standard Deviation == %f  Standard Error == %f" % (mean, stdev, error)

# We need to average all of these gradients together using the exponential
# average (free energy average)
k = 0.00198720650096  # in kcal mol-1 K-1

kT = k*T

avg_grad_f = {}
avg_grad_b = {}
avg_grad = {}

print "\nAverage gradients"
for lam in lams:
    avg_grad_f[lam] = 0
    avg_grad_b[lam] = 0

    count = 0

    for grads_f in all_grads_f:
        dG = grads_f[lam] * delta_lam
        avg_grad_f[lam] += math.exp( -dG / kT )
        count += 1 

    for grads_b in all_grads_b:
        dG = grads_b[lam] * delta_lam
        avg_grad_b[lam] += math.exp( -dG / kT )

    avg_grad_f[lam] = -kT * math.log( avg_grad_f[lam] / count ) / delta_lam
    avg_grad_b[lam] = -kT * math.log( avg_grad_b[lam] / count ) / delta_lam

    avg_grad[lam] = 0.5 * (avg_grad_f[lam] + avg_grad_b[lam])

    print "%f  %f  %f  %f" % (lam, avg_grad_f[lam], avg_grad_b[lam], avg_grad[lam])

avg_grad2 = {}
avg_grad3 = {}

for key in avg_grad.keys():
    avg_grad2[key] = avg_grad[key]
    avg_grad3[key] = avg_grad[key]

# Remove the first and last gradients from avg_grad2
lamvals = avg_grad2.keys()
lamvals.sort()

avg_grad2[lamvals[0]] = avg_grad2[lamvals[2]]
avg_grad2[lamvals[1]] = avg_grad2[lamvals[2]]
avg_grad2[lamvals[-1]] = avg_grad2[lamvals[-3]]
avg_grad2[lamvals[-2]] = avg_grad2[lamvals[-3]]

# Remove the first and second gradients from avg_grad3
avg_grad3[lamvals[0]] =	avg_grad3[lamvals[3]]
avg_grad3[lamvals[1]] =	avg_grad3[lamvals[3]]
avg_grad3[lamvals[2]] = avg_grad3[lamvals[3]]
avg_grad3[lamvals[-1]] = avg_grad3[lamvals[-4]]
avg_grad3[lamvals[-2]] = avg_grad3[lamvals[-4]]
avg_grad3[lamvals[-3]] = avg_grad3[lamvals[-4]]

pmf_f = calculatePMF(avg_grad_f)
pmf_b = calculatePMF(avg_grad_b)
pmf = calculatePMF(avg_grad)
pmf2 = calculatePMF(avg_grad2)
pmf3 = calculatePMF(avg_grad3)

lamvals = pmf_f.keys()
lamvals.sort()

print "\nAverage PMFs"
for lam in lamvals:
    print "%f  %f  %f  %f  %f  %f" % (lam, pmf[lam], pmf_f[lam], pmf_b[lam], pmf2[lam], pmf3[lam])

