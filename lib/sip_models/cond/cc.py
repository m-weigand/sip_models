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
        self.sigma0 = (1 - self.m) * self.sigmai

        # compute some common terms
        self.otc = (self.w * self.tau) ** self.c
        self.otc1 = (self.w * self.tau) ** (self.c - 1)
        self.otc2 = (self.w * self.tau) ** (2 * self.c)
        self.ang = self.c * np.pi / 2.0  # rad
        # numerator and denominator
        self.num = 1 + self.otc * np.cos(self.ang)
        self.denom = 1 + 2 * self.otc * np.cos(self.ang) + self.otc2


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

    def dre_dsigmai(self, pars):
        r"""
        :math:Add formula
        """
        self._set_parameters(pars)
        terms = self.m * self.num / self.denom
        specs = np.sum(terms, axis=1)
        result = 1 - specs

        return result

    def dre_dlog10sigmai(self, pars):
        # first call the linear response to set the parameters
        linear_response = self.dre_dsigmai(pars)
        result = np.log(10) * self.sigmai * linear_response
        return result

    def dre_dm(self, pars):
        r"""
        :math:Add formula
        """
        self._set_parameters(pars)
        terms = self.num / self.denom
        result = - self.sigmai * terms

        return result

    def dre_dlog10m(self, pars):
        # first call the linear response to set the parameters
        lin_response = self.dre_dm(pars)
        result = np.log(10) * self.m * lin_response
        return result

    def dre_dtau(self, pars):
        r"""
        :math:Add formula
        """
        self._set_parameters(pars)
        # term 1
        num1 = self.c * self.w * self.otc1 * np.cos(self.ang)
        term1 = num1/self.denom

        # term 2
        num2a = self.otc * np.cos(self.ang)
        num2b = 1 + num2a
        denom2 = self.denom ** 2
        term2 = num2b / denom2

        # term 3
        term3 = 2 * self.c * self.w * self.otc1 * np.cos(self.ang) + self.otc2

        result = self.sigmai * self.m * (term1 + term2 * term3)

        return result

    def dre_dlog10tau(self, pars):
        # first call the linear response to set the parameters
        lin_response = self.dre_dtau(pars)
        result = np.log(10) * self.tau * lin_response
        return result

    def dre_dc(self, pars):
        r"""
        :math:Add formula
        """
        self._set_parameters(pars)
        # term 1
        num1a = np.log(self.w * self.tau) * self.otc * np.sin(self.ang)
        num1b = self.otc * np.cos(self.ang) * np.pi / 2.0
        term1 = (num1a + num1b) / self.denom

        # term 2
        num2 = self.otc * np.sin(self.c / np.pi) * 2
        denom2 = self.denom ** 2
        term2 = num2 / denom2

        # term 3
        num3a = 2 * np.log(self.w * self.tau) * self.otc * np.cos(self.ang)
        num3b = 2 * ((self.w * self.tau) ** 2) * np.pi / 2.0 * np.sin(self.ang)
        num3c = 2 * np.log(self.w * self.tau) * self.otc2
        term3 = num3a - num3b + num3c

        result = self.sigmai * self.m * (term1 + term2 * term3)

        return result

    def dim_dsigmai(self, pars):
        r"""
        :math:Add formula
        """
        self._set_parameters(pars)
        result = np.sum(- self.m * self.otc * np.sin(self.ang) / self.denom,
                        axis=1)

        return result

    def dim_dlog10sigmai(self, pars):
        # first call the linear response to set the parameters
        lin_response = self.dim_dsigmai(pars)
        result = np.log(10) * self.sigmai * lin_response
        return result

    def dim_dm(self, pars):
        r"""
        :math:Add formula
        """
        self._set_parameters(pars)
        num1 = self.otc * np.sin(self.ang)
        result = -self.sigmai * num1 / self.denom

        return result

    def dim_dlog10m(self, pars):
        # first call the linear response to set the parameters
        lin_response = self.dim_dm(pars)
        result = np.log(10) * self.m * lin_response
        return result

    def dim_dtau(self, pars):
        r"""
        :math:Add formula
        """
        self._set_parameters(pars)
        # term 1
        num1 = -self.m * (self.w ** self.c) * self.c\
            * (self.tau ** (self.c - 1)) * np.sin(self.ang)
        term1 = self.sigmai * num1 / self.denom

        # term 2
        num2a = -self.m * self.otc * np.sin(self.ang)
        num2b = 2 * (self.w ** 2.0) * self.c * (self.tau ** (self.c - 1)) *\
            np.cos(self.ang)
        num2c = 2 * self.c * (self.w ** (self.c * 2)) *\
            (self.tau ** (2 * self.c - 1))
        term2 = self.sigma0 * num2a * (num2b + num2c) / (self.denom ** 2)

        result = term1 + term2

        return result

    def dim_dlog10tau(self, pars):
        # first call the linear response to set the parameters
        lin_resp = self.dim_dtau(pars)
        result = np.log(10) * self.tau * lin_resp
        return result

    def dim_dc(self, pars):
        r"""
        :math:Add formula
        """
        self._set_parameters(pars)
        # term 1
        num1a = self.m * np.sin(self.ang) * np.log(self.w * self.tau)\
            * self.otc
        num1b = self.m * self.otc * np.pi / 2 * np.cos(np.pi / 2)
        term1 = self.sigma0 * (-num1a - num1b) / self.denom

        # term 2
        num2a = -self.m * self.otc * np.cos(self.ang)
        num2b = -2 * np.log(self.w * self.tau) * self.otc * np.cos(self.ang)
        num2c = 2 * self.otc * np.pi / 2 * np.cos(self.ang)
        num2d = 2 * np.log(self.w * self.tau) * self.otc2
        numerator = num2a * (num2b + num2c) + num2d
        term2 = self.sigma0 * numerator / (self.denom ** 2)

        result = term1 + term2

        return result

    def test_derivatives(self):
        parameters = {
            'sigmai': 0.01,
            'm': 0.1,
            'tau': 0.04,
            'c': 0.8
        }
        # parameters = {
        #     'sigmai': 0.01,
        #     'm': (0.15, 0.2),
        #     'tau': (0.4, 0.004),
        #     'c': (0.5, 0.8),
        # }
        print(self.dre_dsigmai(parameters))
        print(self.dre_dm(parameters))
        print(self.dre_dtau(parameters))
        print(self.dre_dc(parameters))

        print(self.dim_dsigmai(parameters))
        print(self.dim_dm(parameters))
        print(self.dim_dtau(parameters))
        print(self.dim_dc(parameters))
