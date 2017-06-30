#!/usr/bin/env python

# This example shows how to create SIP-spectra with Cole-Cole conductivity model 'cond.cc' for different input formats.

import sip_models.cond.cc
import numpy as np

# Create a vector with needed frequency range
frequencies = np.logspace(-3, 4, 22)
ccobj = sip_models.cond.cc.cc(frequencies)
np.savetxt('cond_cc_frequencies.dat', frequencies)

## Example 1
# Create the SIP-spectrum for a given parameterset and the given frequency range
# Order of parameters: sigma0, m1, tau1, c1
result = ccobj.response([0.01, 0.1, 0.04, 0.6])

# plot the results and save the results to file
result.plot('cond_cc_01_plot.png')
np.savetxt('cond_cc_01_data.dat', result.to_one_line(result.rmag_rpha))

## Example 2: with two peaks
# Order of parameters: sigma0, m1, m2, tau1, tau2, c1, c2
result2 = ccobj.response([0.01, 0.1, 0.2, 0.04, 0.0001, 0.4, 0.8])

# plot the results and save the results to file
result2.plot('cond_cc_02_plot.png')
np.savetxt('cond_cc_02_data.dat', result.to_one_line(result2.rmag_rpha))
