#!/usr/bin/env python
"""Simple example of using the  sip_response object for visualizing an SIP
signature
"""
import sip_models.res.cc
import numpy as np


def generate_sip_data1():
    # generate a SIP signature
    frequencies = np.logspace(-3, 4, 22)
    ccobj = sip_models.res.cc.cc(frequencies)
    result = ccobj.response([100, 0.1, 0.04, 0.6])
    return result


def generate_sip_data2():
    # generate a SIP signature
    frequencies = np.logspace(-3, 4, 22)
    ccobj = sip_models.res.cc.cc(frequencies)
    result = ccobj.response([100, 0.1, 0.4, 0.6])
    return result


if __name__ == '__main__':
    print('generate sip spectrum')
    sip_obj = generate_sip_data1()
    # plot spectrum to file
    sip_obj.plot(
        filename='plot_spectrum_simple_rho.png',
    )

    sip_obj.plot(
        filename='plot_spectrum_simple_R.png',
        dtype='R',
    )

    # generate another sip signature
    sip_obj2 = generate_sip_data2()
    sip_obj.plot(
        filename='plot_spectrum_reciprocal.png',
        reciprocal=sip_obj2
    )

    sip_obj.plot(
        filename='plot_spectrum_reciprocal_limits.png',
        reciprocal=sip_obj2,
        limits={
            'rmag_min': 80,
            'rpha_min': -40,
        },
    )
