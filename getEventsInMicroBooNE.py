import math
import numpy as np
import Functions.EventRates as ER


POT = 6e20 # protons on target, three years livetime [uBooNE proposal]
mass = 350 # MeV
theta_mu = 1e-9

microBooneEvents = ER.EventRate_Nu2MuPi(mass,theta_mu,POT)[1]

print microBooneEvents
