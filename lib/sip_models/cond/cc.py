# *-* coding: utf-8 *-*
"""
Cole-Cole model after Tarasov and Titov, 2013

Tarasov, A. and Titov, K. (2013). On the use of the Cole–Cole equations
in spectral induced polarization. Geophys J Int, 195(1):352-356.
doi: 10.1093/gji/ggt251
"""
import numpy as np
import sip_models.sip_response as sip_response

def _make_list(number_or_list):
    # return the object enclosed in a list if its not a tuple or list
    if isinstance(number_or_list, (tuple, list)):
        return number_or_list
    else:
        return [number_or_list, ]

class cc_base(object):
    """
    Base class for Cole-Cole objects (conductivity)
    """
    def __init__(self, frequencies):
        self.f = frequencies

    def _sort_parameters(self, parameters):
        # type 1
        if isinstance(parameters, (list, tuple, np.ndarray)):
            pars = np.atleast_1d(parameters)
            nr_pars = int((pars.shape[0] - 1) / 3)

            sigmai = pars[0]
            m = pars[1:nr_pars + 1]
            tau = pars[nr_pars + 1: 2 * nr_pars + 1]
            c = pars[2 * nr_pars + 1:]

        elif isinstance(parameters, dict):
            sigmai = parameters['sigmai']
            m = _make_list(parameters['m'])
            tau = _make_list(parameters['tau'])
            c = _make_list(parameters['c'])
        else:
            print(parameters)
            raise Exception('Input format not recognized')

        return sigmai, m, tau, c


    def _set_parameters(self, parameters):
        """Sort out the various possible parameter inputs and return a config
        object (dict)

        We have multiple input formats:

        1) a list, tuple, or numpy.ndarray, containing the linear parameters
        in the following order:
            for single term: sigmai, m1, tau1, c1
            for multiple terms: sigmai, m1, m2, ..., tau1, tau2, ..., c1, c2,
            ...

        2) a dictionary with the entries "sigmai", "m", "tau", "c"

        2b) if the dictionary entries for "m", "tau", and "c" are lists, the
        entries correspond to mulitple polarisazion terms

        """
        nr_f = self.f.size

        # sort out parameters
        sigmai, m, tau, c = self._sort_parameters(parameters)

        newsize = (nr_f, len(m))
        # sigmai_resized = np.resize(sigmai, newsize)
        m_resized = np.resize(m, newsize)
        tau_resized = np.resize(tau, newsize)
        c_resized = np.resize(c, newsize)

        omega = np.atleast_2d(2 * np.pi * self.f).T
        self.w = np.resize(omega, (len(m), nr_f)).T
        self.sigmai = sigmai
        self.m = m_resized
        self.tau = tau_resized
        self.c = c_resized

        # compute some common terms
        #self.otc = (self.w * self.tau) ** self.c
        #self.otc2 = (self.w * self.tau) ** (2 * self.c)
        #self.ang = self.c * np.pi / 2.0  # rad
        #self.denom = 1 + 2 * self.otc * np.cos(self.ang) + self.otc2


class cc(cc_base):

    def response(self, parameters):
        r"""Return the forward response in base dimensions
        :math:`\hat{\sigma }(\omega ) = \sigma _\infty \left(1 - \sum_i \frac
        {m_i}{1 + (j \omega \tau_i)^c_i}\right)`

        Parameters
        ----------
        pars:

        Returns
        -------
        response: Nx2 array, first axis denotes frequencies, seconds real and
                  imaginary parts
        """
        # get a config object
        self._set_parameters(parameters)
        terms = self.m / (1 + (1j * self.w * self.tau) ** self.c)
        # sum up terms
        specs = np.sum(terms, axis=1)
        ccomplex = self.sigmai * (1 - specs)

        response = sip_response.sip_response(self.f, ccomplex=ccomplex)
        
        return response