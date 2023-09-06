# copyright 2023 - Mel Davey

import numpy as np

#
# The 'Model' class approximates the physical system of DNA strand incorporation and UV cleavage
# here, we account for the following:
# incomplete extention (ie) where not all available binding sites incorporate a nucleotide 
# carry-forward (cf) where further incorporation can occur after the tagged nucleotide is cleaved
# droop (dr) where some DNA strands can no longer support new incorporations
#

class Model:
    def __init__(self, strandLen:int = 100):
        self.strandLen = strandLen

        self.state = np.zeros(strandLen)
        self.prevState = np.zeros(strandLen)

        self.state[0] = 1.0
        self.dr = 0.01
        self.ie = 0.05
        self.cf = 0.07
        self.extraBucket = 0
        self.bases = ['A', 'C', 'G', 'T']

    def SetParams(self, ie:float = 0, cf:float = 0, dr:float = 0):
        self.ie = ie
        self.cf = cf
        self.dr = dr

    def ApplyUV(self, maxLen:int = 0):
        # save current state
        self.prevState = self.state[:]
        self.extraBucket = 0
        numExtensions = 3 # technically this goes on forever, but after 3 rounds there is not much left
        extendAmount = np.zeros(numExtensions)

        # apply incompletion, and carry-forward, and signal loss effects to our state
        # note - this runs in reverse so we don't double-count advancing product
        for i in range(self.strandLen-1, -1, -1):
            # amount of product available to extend
            amount = self.state[i] * (1.0 - self.ie)
            if maxLen > 0 and i >= (maxLen-1):
                amount = 0
            if amount == 0:
                continue

            # divide up the amount into various extensions
            extendAmount[0] = amount
            for s in range(1,numExtensions):
                cfAmount = extendAmount[s-1] * self.cf
                extendAmount[s] = cfAmount
                extendAmount[s-1] -= extendAmount[s]

            # update state
            for s in range(numExtensions, 0, -1):
                extendState = i + s
                if extendState >= (self.strandLen-1):
                    self.extraBucket += extendAmount[s-1]
                else:
                    self.state[extendState] += extendAmount[s-1]

            self.state[i] -= amount

        # droop is applied across all states
        self.state *= (1.0 - self.dr)


    def GetSignal(self, dnaTemplate):
        signal = np.zeros(6) # 4 bases plus unknown & extra
        templateLen = len(dnaTemplate)
        # each position within our state, sum the signal, binned by known DNA bases
        for i in range(self.strandLen):
            if i < templateLen:
                signal[self.bases.index(dnaTemplate[i])] += self.state[i]
            else:
                signal[4] += self.state[i]
        signal[5] = self.extraBucket
        return signal

    def Revert(self):
        # revert to previous state
        self.state = self.prevState[:]

    def GetState(self):
        return self.state

    def Reset(self):
        self.state.fill(0)
        self.state[0] = 1.0

