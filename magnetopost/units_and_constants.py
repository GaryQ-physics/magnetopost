# units_and_constants
import numpy as np

# base units
hr = 1.
muA = 1.
R_e = 1.
nT = 1.

deg = (np.pi/180.)
amin = deg/60.
minn = hr/60.
s = minn/60.

A = 1e+6*muA
Tesla = 1e+9*nT
m = R_e/6.3781e+6  # R_e == 6.3781e+6*m  (note m < 10^-6 Re)        6371.0008?????
kg = Tesla*A*s**2
Gauss = 1e-4*Tesla
Coulomb = A/s
V = (kg*m**2/s**2)/Coulomb
mV = 1e-3*V
Ohm = V/A
Siemens = 1/Ohm

mu0_SI = 1.25663706212e-6 # almost 4*pi*1e-7
if False:
    mu0 = mu0_SI*kg*m/((s**2)*(A**2))
else:
    mu0 = 1.970237314e-10*nT*R_e/muA

phys = {
        'hr': hr,
        'muA': muA,
        'R_e': R_e,
        'nT': nT,

        'deg' : (np.pi/180.),
        'amin' : deg/60.,
        'minn' : hr/60.,
        's' : minn/60.,

        'A' : A,
        'T' : Tesla,
        'Gauss' : Gauss,
        'm' : m,
        'kg' : kg,
        'Coulomb' : Coulomb,
        'V'       : V,
        'mV'      : mV,
        'Ohm' : Ohm,
        'Siemens' : Siemens,


        'mu0_SI' : mu0_SI,
        'mu0' : mu0
    }

phys['c'] = 299792458.*phys['m']/phys['s']
