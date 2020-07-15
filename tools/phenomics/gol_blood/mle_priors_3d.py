#!/usr/bin/env python


class MLEPriors3D:

    def __init__(self, cMean, cStd, rMean, rStd, eMean, eStd):
        self.xMean = cMean[0]
        self.yMean = cMean[1]
        self.zMean = cMean[2]
        self.xStd = cStd[0]
        self.yStd = cStd[1]
        self.zStd = cStd[2]
        self.rMean = rMean
        self.rStd = rStd
        self.xEMean = eMean[0]
        self.xEStd = eStd[0]
        # Need to confirm with Prof. John Liechty. It should be eMean[1], in my opinion
        self.yEMean = cMean[1]
        self.yEStd = eStd[1]
        # Need to confirm with Prof. John Liechty. It should be eMean[2], in my opinion
        self.zEMean = cMean[2]
        self.zEStd = eStd[2]
        self.rSigShape = 50
        self.rSigScale = 0.5
        # Outside of 3 std reference dist dominates 
        self.stdLevelRefference = 3
