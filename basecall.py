# copyright 2023 - Mel Davey

import model
import matplotlib.pyplot as plt
import numpy as np
import json
import sys

verbose = 0
wantPlots = False
gridsearch = False

bases = ['A', 'C', 'G', 'T']
base_colors = ['green', 'yellow', 'blue', 'red']

# ie - the percent of product that does not incorporate
# cf - the percent of product that can incorpporate subsequent positions during UV cleavage
# dr - ther percent of product lost at each UV cycle (modeled as a system-wide param)
ie = 0.08
cf = 0.06
dr = 0.02

# easy to create alternate models to represent the physical system
# state_model selects which one to use

# if set to true, we will perform a fit across the measured data, and correct out any signal loss
# note that the dr (droop) param is then set to 0.0 for the model predictions
correctLoss = False

argc = len(sys.argv)
argcc = 1
while argcc < argc:
    if sys.argv[argcc] == '-ie':
        argcc += 1
        ie = float(sys.argv[argcc])
    if sys.argv[argcc] == '-cf':
        argcc += 1
        cf = float(sys.argv[argcc])
    if sys.argv[argcc] == '-dr':
        argcc += 1
        dr = float(sys.argv[argcc])
    if sys.argv[argcc] == '-v':
        verbose += 1
    if sys.argv[argcc] == '-plots':
        wantPlots = True
    if sys.argv[argcc] == '-grid':
        gridsearch = True
    if sys.argv[argcc] == '-loss':
        correctLoss = True
    argcc += 1

#
# CallBases
#
# Uses a physical model to predict what the measured signal would be after a UV cycle is applied
# performs predictions with an initially blank DNA template, so on the first pass only
# incompletion and signal loss (droop) can be accounted for.  On the second and subsequent
# passes, the model is able to improve signal predictions because it has a rough idea of the DNA
# template being sequenced, carry-forward in particual can now be accounted for correctly
#

def CallBases(ie, cf, dr, numCycles, measuredSignal):
    dnaTemplate = ''
    m = model.Model()
    m.SetParams(ie=ie, cf=cf, dr=dr)
    cumulativeError = 0
    dyeIntensities = np.zeros((numCycles, 4))
    totalSignal = np.zeros(numCycles)
    errorPerCycle = np.zeros(numCycles)

    # perform base calling
    numIterations = 3
    for iteration in range(numIterations):
        m.Reset()
        for cycle in range(numCycles):
            best_error = 0
            best_base = -1
            best_signal = 0

            for base in bases:
                # insert a "what-if" base at the current position, and predict what the signal looks like
                testTemplate = dnaTemplate[:cycle] + base + dnaTemplate[cycle+1:]

                # the model will return the predicted signal across all 4 dyes
                # plus an unknown and extra component that the model is also tracking
                # for example when we are beyond the end of the template, signals get bucketed up differently at the end
                signal = m.GetSignal(testTemplate)
                signalSum = np.sum(signal[:4]) # total intensity of the 4 dyes

                # compare to measured at this cycle across all 4 dyes
                error = 0
                for i in range(4):
                    delta = (measuredSignal[cycle][i] - signal[i])/signalSum
                    error += delta*delta

                # keep track of the lowest error, this is the best predition
                if error < best_error or best_base == -1:
                    best_base = base
                    best_error = error
                    best_signal = signal

            # append/replace with best base at current position (cycle)
            dnaTemplate = dnaTemplate[:cycle] + best_base + dnaTemplate[cycle+1:]
            dyeIntensities[cycle] = best_signal[:4]
            totalSignal[cycle] = np.sum(best_signal[:4])
            errorPerCycle[cycle] = best_error

            # update the model - note that we do this after getting the measured signals, because this matches the physical
            # system where the first base is exposed to nucleotides prior to UV cleavage
            m.ApplyUV(numCycles)

        if verbose > 0:
            print('iteration %d basecalls: %s' % (iteration, dnaTemplate))

    print('basecalls: %s' % dnaTemplate)
    cumulativeError = np.sum(errorPerCycle)
    return {'err':cumulativeError, 'basecalls':dnaTemplate, 'intensites':dyeIntensities, 'signal':totalSignal, 'error':errorPerCycle}

def CorrectSignalLoss(measuredSignal):
    totalMeasuredSignal = np.sum(measuredSignal, axis=1)
    loss_dim = 2 # 1 for linear, 2 for quadratic, etc
    X = np.arange(len(totalMeasuredSignal))
    coef = np.polyfit(X, totalMeasuredSignal, loss_dim)
    print('measured loss: %s' % coef)
    fit_fn = np.poly1d(coef)
    lossCorrectedSignal = np.copy(measuredSignal)
    for cycle in range(len(totalMeasuredSignal)):
        lossCorrectedSignal[cycle] /= fit_fn(cycle)
    return lossCorrectedSignal


def GridSearch():
    cf1 = 0.05
    cf2 = 0.08
    cfnum = 11
    ie1 = 0.07
    ie2 = 0.11
    ienum = 11
    dr1 = 0.01
    dr2 = 0.025
    drnum = 4

    minerr = 99999
    bestie = 0
    bestcf = 0
    bestdr = 0
    for cf in np.linspace(cf1, cf2, cfnum):
        for ie in np.linspace(ie1, ie2, ienum):
            for dr in np.linspace(dr1, dr2, drnum):
                res = CallBases(ie, cf, dr, numCycles, data)
                if res['err'] < minerr:
                    minerr = res['err']
                    bestie = ie
                    bestcf = cf
                    bestdr = dr
    print('best err:%f ie:%f cf:%f dr:%f' % (minerr, bestie, bestcf, bestdr))
    return bestie, bestcf, bestdr

#
# main starts here
#

# load up the test dataset
with open('jmrdata1.json') as f:
    data = np.array(json.load(f)).astype(np.float64)
data /= 100.0 # hacky normalization
numCycles = data.shape[0]
print('cycles: %d' % numCycles)

if gridsearch:
    ie,cf,dr = GridSearch()

if correctLoss:
    measuredSignal = CorrectSignalLoss(data)
    dr = 0.0
else:
    measuredSignal = data

results = CallBases(ie, cf, dr, numCycles, measuredSignal)
print('cumulative error: %f' % results['err'])

# plots
if wantPlots:
    fig, ax = plt.subplots()
    fig.suptitle('predicted signals per cycle')
    for cycle in range(numCycles):
        for base in range(4):
            ax.bar(cycle + base*0.1, results['intensites'][cycle, base], color = base_colors[base], width = 0.1)
    plt.plot(range(numCycles), results['signal'], label='total intensity')
    plt.plot(range(numCycles),  results['error'], label='error')
    plt.legend(loc="upper right")


    totalMeasuredSignal = np.sum(measuredSignal, axis=1)
    fig, ax = plt.subplots()
    fig.suptitle('measured signals per cycle')
    for cycle in range(numCycles):
        for base in range(4):
            ax.bar(cycle + base*0.1, measuredSignal[cycle, base], color = base_colors[base], width = 0.1)
    plt.plot(range(numCycles), totalMeasuredSignal, label='total intensity')
    plt.show()

