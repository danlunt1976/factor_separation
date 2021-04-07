"""
Need to do: 
module add anaconda/2.7-5.0.1
"""

import netCDF4
from netCDF4 import Dataset

import numpy as np

e280 = 1
e400 = 2
ei280 = 3
ei400 = 4
eo280 = 5
eo400 = 6
eoi280 = 7
eoi400 = 8

AllClims = {}

AllClims[e280] = Dataset("E280.nc", "r")
AllClims[e400] = Dataset("E400.nc", "r")
AllClims[ei280] = Dataset("Ei280.nc", "r")
AllClims[ei400] = Dataset("Ei400.nc", "r")
AllClims[eo280] = Dataset("Eo280.nc", "r")
AllClims[eo400] = Dataset("Eo400.nc", "r")
AllClims[eoi280] = Dataset("Eoi280.nc", "r")
AllClims[eoi400] = Dataset("Eoi400.nc", "r")

var = "atm/TREFHT"

eoi400d=np.float64(AllClims[eoi400][var][:, :])
ei400d=np.float64(AllClims[ei400][var][:, :])
eo400d=np.float64(AllClims[eo400][var][:, :])
e400d=np.float64(AllClims[e400][var][:, :])
eoi280d=np.float64(AllClims[eoi280][var][:, :])
ei280d=np.float64(AllClims[ei280][var][:, :])
eo280d=np.float64(AllClims[eo280][var][:, :])
e280d=np.float64(AllClims[e280][var][:, :])

eoi400s=AllClims[eoi400][var][:, :]
ei400s=AllClims[ei400][var][:, :]
eo400s=AllClims[eo400][var][:, :]
e400s=AllClims[e400][var][:, :]
eoi280s=AllClims[eoi280][var][:, :]
ei280s=AllClims[ei280][var][:, :]
eo280s=AllClims[eo280][var][:, :]
e280s=AllClims[e280][var][:, :]


def factorize_Lunt2012(eoi400,ei400,eo400,e400,eoi280,ei280,eo280,e280):
    """
    Factorize an anomaly into contributions from each of the different changes to the boundary conditions
    using the methodology of Lunt et al. 2012 and Haywood et al. 2016
    
    Returns:
        1. The total change between pliocene and PI
        2. The factorized contribution from CO2
        3. The factorized contribution from topography
        4. The factorized contribution from ice
    """

    dT = eoi400 - e280

    dTCO2 = (e400 - e280 +
             eo400 - eo280 +
             ei400 - ei280 +
             eoi400 - eoi280) / 4

    dTtopo = (eo280 - e280 +
              eo400 - e400 +
              eoi280 - ei280 +
              eoi400 - ei400) / 4

    dTice = (ei280 - e280 +
             ei400 - e400 +
             eoi280 - eo280 +
             eoi400 - eo400) / 4

    return dT, dTCO2, dTtopo, dTice


def factorizeLinearSum(eoi400,ei400,eo400,e400,eoi280,ei280,eo280,e280):
    """
    Factorize an anomaly into contributions from each of the different changes to the boundary conditions
    using the all-linear and complete methodology of Lunt et al. 2020.
    
    Returns:
        1. The total change between pliocene and PI
        2. The factorized contribution from CO2
        3. The factorized contribution from topography
        4. The factorized contribution from ice
    """

    var = "atm/TREFHT"
    dT = eoi400 - e280

    dTCO2 = (2 * (e400 - e280) +
             (eo400 - eo280) +
             (ei400 - ei280) +
             2 * (eoi400 - eoi280)) / 6

    dTtopo = (2 * (eo280 - e280) +
              (eo400 - e400) +
              (eoi280 - ei280) +
              2 * (eoi400 - ei400)) / 6

    dTice = (2 * (ei280 - e280) +
             (ei400 - e400) +
             (eoi280 - eo280) +
             2 * (eoi400 - eo400)) / 6

    return dT, dTCO2, dTtopo, dTice


def factorizeScaledResidualRel(eoi400,ei400,eo400,e400,eoi280,ei280,eo280,e280):
    """
    Factorize an anomaly into contributions from each of the different changes to the boundary conditions
    using the all-linear and complete methodology of Lunt et al. 2020.
    
    Returns:
        1. The total change between pliocene and PI
        2. The factorized contribution from CO2
        3. The factorized contribution from topography
        4. The factorized contribution from ice
    """
    var = "atm/TREFHT"
    dT, dTCO2, dTtopo, dTice = factorize_Lunt2012(eoi400,ei400,eo400,e400,eoi280,ei280,eo280,e280)

    res = ((3 * eoi400 + ei400 +
            eoi280 + eo400) -
           (3 * e280 + ei280 +
            eo280 + e400)) / 4

    dTCO2new = dTCO2 * (dT / res)
    dTtoponew = dTtopo * (dT / res)
    dTicenew = dTice * (dT / res)

    return dT, dTCO2, dTtopo, dTice, dTCO2new, dTtoponew, dTicenew


def factorizeScaledResidualAbs(eoi400,ei400,eo400,e400,eoi280,ei280,eo280,e280):
    """
    Factorize an anomaly into contributions from each of the different changes to the boundary conditions
    using the all-linear and complete methodology of Lunt et al. 2020.
    
    Returns:
        1. The total change between pliocene and PI
        2. The factorized contribution from CO2
        3. The factorized contribution from topography
        4. The factorized contribution from ice
    """
    var = "atm/TREFHT"
    dT, dTCO2, dTtopo, dTice = factorize_Lunt2012(eoi400,ei400,eo400,e400,eoi280,ei280,eo280,e280)

    res = eoi400-e280 - (dTCO2+dTtopo+dTice)

    dTCO2new = dTCO2 + (res*abs(dTCO2) / (abs(dTCO2)+abs(dTtopo)+abs(dTice)))
    dTtoponew = dTtopo + (res*abs(dTtopo) / (abs(dTCO2)+abs(dTtopo)+abs(dTice)))
    dTicenew = dTice + (res*abs(dTice) / (abs(dTCO2)+abs(dTtopo)+abs(dTice)))


    return dT, dTCO2, dTtopo, dTice, dTCO2new, dTtoponew, dTicenew



dT, dTCO2new, dTtoponew, dTicenew = factorizeLinearSum(eoi400d,ei400d,eo400d,e400d,eoi280d,ei280d,eo280d,e280d)
print(np.amax(dTCO2new))
print(np.amax(dTicenew))
print(np.amax(dTtoponew))

dT, dTCO2, dTtopo, dTice, dTCO2new, dTtoponew, dTicenew = factorizeScaledResidualRel(eoi400d,ei400d,eo400d,e400d,eoi280d,ei280d,eo280d,e280d)
print(np.amax(dTCO2))
print(np.amax(dTice))
print(np.amax(dTtopo))
print(np.amax(dTCO2new))
print(np.amax(dTicenew))
print(np.amax(dTtoponew))

dT, dTCO2, dTtopo, dTice, dTCO2new, dTtoponew, dTicenew = factorizeScaledResidualAbs(eoi400d,ei400d,eo400d,e400d,eoi280d,ei280d,eo280d,e280d)
print(np.amax(dTCO2new))
print(np.amax(dTicenew))
print(np.amax(dTtoponew))

