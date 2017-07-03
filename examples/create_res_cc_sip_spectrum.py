#!/usr/bin/env python

# This example shows how to create SIP-spectra with Cole-Cole resistivity model 'res.cc' for different input formats.

import sip_models.res.cc
import numpy as np

# Create a vector with needed frequency range
frequencies = np.logspace(-3, 4, 22)
ccobj = sip_models.res.cc.cc(frequencies)
np.savetxt('res_cc_frequencies.dat', frequencies)

###############################################
## Example 1
# Create the SIP-spectrum for a given parameterset and the given frequency range
# Order of parameters: rho0, m1, tau1, c1
result = ccobj.response([100, 0.1, 0.04, 0.6])

# plot the results and save the results to file
result.plot('res_cc_01_plot.png')
np.savetxt('res_cc_01_data.dat', result.to_one_line(result.rmag_rpha))

###############################################
## Example 2: with two peaks
# Order of parameters: rho0 m1 m2 tau1 tau2 c1 c2
result2 = ccobj.response([100, 0.1, 0.2, 0.04, 0.0001, 0.4, 0.8])

# plot the results and save the results to file
result2.plot('res_cc_02_plot.png')
np.savetxt('res_cc_02_data.dat', result.to_one_line(result2.rmag_rpha))

###############################################
## Example 3: One peak, different parameter format
# Format of parameters: dict
cc_parameters_ex03 = {
	'rho0': 100,
	'm': 0.1,
	'tau': 0.04,
	'c': 0.6
}
result2 = ccobj.response(cc_parameters_ex03)

# plot the results and save the results to file
result2.plot('res_cc_03_plot.png')
np.savetxt('res_cc_03_data.dat', result.to_one_line(result2.rmag_rpha))

###############################################
## Example 4: Two peaks, different parameter format
# Format of parameters: dict
cc_parameters_ex04 = {
	'rho0': 100,
	'm': (0.1, 0.1),
	'tau': (0.4, 0.004),
	'c': (0.6, 0.6)
}
result2 = ccobj.response(cc_parameters_ex04)

# plot the results and save the results to file
result2.plot('res_cc_04_plot.png')
np.savetxt('res_cc_04_data.dat', result.to_one_line(result2.rmag_rpha))
