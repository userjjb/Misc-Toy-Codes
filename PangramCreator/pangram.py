import numpy as np                                                                                                                                      
from collections import Counter
import time
import os
from timeit import default_timer as timer


N = 4

AtoZ = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

NumbersText = ['ZERO', 'ONE', 'TWO', 'THREE', 'FOUR', 'FIVE', 'SIX', 'SEVEN', 'EIGHT', 'NINE', 'TEN', 'ELEVEN', 'TWELVE', 'THIRTEEN', 'FOURTEEN','FIFTEE    N', 'SIXTEEN', 'SEVENTEEN', 'EIGHTEEN', 'NINETEEN', 'TWENTY', 'TWENTYONE', 'TWENTYTWO', 'TWENTYTHREE', 'TWENTYFOUR', 'TWENTYFIVE']

NumbersCount = np.zeros((26,26))
for i,n in enumerate(NumbersText):
    cnt = Counter(AtoZ + n)
    NumbersCount[i,:] = np.array(list(cnt.values()))[np.argsort(list(cnt.keys()))] - 1
NumbersCount = NumbersCount.astype(int)

SeedCount = [34, 1, 33, 7, 194, 20, 9, 24, 63, 1, 1, 7, 2, 98, 87, 55, 1, 97, 50, 99, 11, 19, 15, 5, 2, 26]
SeedCount += np.sum(NumbersCount[np.random.randint(10, size=26), :], axis=0)

Percent = 100*SeedCount/np.sum(SeedCount)
Actual = Percent.copy()

nLetterF = [34, 1, 33, 7, 94, 3, 1, 6, 36, 1, 1, 7, 2, 61, 33, 55, 1, 56, 37, 65, 3, 1, 3, 1, 1, 3]

LoopL = 100
OuterEnd = np.int(np.round(9090/LoopL)*80)
nLetter = SeedCount.copy()
start=timer()
for Mutate in range(OuterEnd):
    for tries in range(LoopL):
        nLetterOld = nLetter.copy()
        nLetter = (nLetterF
                + np.sum(NumbersCount[np.floor(Percent).astype(int),:], axis=0)
                + np.sum(NumbersCount[np.floor(10*(Percent - np.floor(Percent))).astype(int),:], axis=0)
                + np.sum(NumbersCount[np.floor(10*(10*Percent - np.floor(10*Percent))).astype(int),:], axis=0)
                + np.sum(NumbersCount[np.floor(10*(100*Percent - np.floor(100*Percent))).astype(int),:], axis=0)
                + np.sum(NumbersCount[np.round(10*(1000*Percent - np.floor(1000*Percent))).astype(int),:], axis=0))
        Actual = 100*nLetter/np.sum(nLetter)

        if np.sum(np.abs(nLetter-nLetterOld))<3:
            print('Diff: ' + str(np.sum(np.abs(nLetter-nLetterOld))))
            print(nLetter)
            print(Percent)
            print(Actual)

        Percent = Actual.copy()

    replace = np.random.randint(10,size=1)
    delete = np.random.randint(10,size=1)
    while delete==replace:
        delete = np.random.randint(10,size=1)
    nLetter += np.squeeze(NumbersCount[replace, :] - NumbersCount[delete, :])
    Percent = 100*nLetter/np.sum(nLetter)
end = timer()
print(end-start)
#print('PROCID: ' + os.environ['SLURM_PROCID'] + ' complete!')
