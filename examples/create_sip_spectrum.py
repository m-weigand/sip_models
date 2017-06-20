#!/usr/bin/env python
import sip_models.res.cc
import numpy as np

frequencies = np.logspace(-3, 4, 22)
ccobj = sip_models.res.cc.cc(frequencies)

result = ccobj.response([100, 0.1, 0.04, 0.6])
result.plot('test.png')

# save to file
np.savetxt('frequencies.dat', frequencies)
np.savetxt('data.dat', result.to_one_line(result.rmag_rpha))

# create another one with two peaks
# note the order of parameters: rho0 m1 m2 tau1 tau2 c1 c2
result2 = ccobj.response([100, 0.1, 0.2, 0.04, 0.0001, 0.4, 0.8])
result2.plot('test2.png')

np.savetxt('data2.dat', result.to_one_line(result2.rmag_rpha))
