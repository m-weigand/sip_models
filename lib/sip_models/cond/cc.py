# *-* coding: utf-8 *-*
"""
Cole-Cole model after Tarasov and Titov, 2013

[CITATION]
"""
import numpy as np
import sip_models.sip_response as sip_response


class cc(object):

    def __init__(self, frequencies):
        self.frequencies = frequencies

    def response(self, pars):
        """Return the forward response in base dimensions

        Parameters
        ----------
        pars:

        Returns
        -------
        response: Nx2 array, first axis denotes frequencies, seconds real and
                  imaginary parts
        """
        sigmai = pars[0]
        m = pars[1]
        tau = pars[2]
        c = pars[3]

        omega = 2 * np.pi * self.frequencies

        relterm = m / (1 + (1j * omega * tau) ** c)
        ccomplex = sigmai * (1 - relterm)

        response = sip_response.sip_response(
            self.frequencies, ccomplex=ccomplex
        )
        return response
