# *-* coding: utf-8 *-*
""" Cole-Cole model (resistivity/resistance formulation) after Pelton et al.
1978

Pelton, W., Ward, S., Hallof, P., Sill, W., and Nelson, P. (1978). Mineral
discrimination and removal of inductive coupling with multifrequency ip.
Geophysics, 43(3):588–609.
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
    """ Base class for Cole-Cole objects (both resistivity and conductivity)
    """
    def __init__(self, frequencies):
        self.f = frequencies

    def _sort_parameters(self, parameters):
        # type 1
        if isinstance(parameters, (list, tuple, np.ndarray)):
            pars = np.atleast_1d(parameters)
            nr_pars = int((pars.shape[0] - 1) / 3)

            rho0 = pars[0]
            m = pars[1:nr_pars + 1]
            tau = pars[nr_pars + 1: 2 * nr_pars + 1]
            c = pars[2 * nr_pars + 1:]

        elif isinstance(parameters, dict):
            rho0 = parameters['rho0']
            m = _make_list(parameters['m'])
            tau = _make_list(parameters['tau'])
            c = _make_list(parameters['c'])
        else:
            print(parameters)
            raise Exception('Input format not recognized')

        return rho0, m, tau, c

    def _set_parameters(self, parameters):
        """Sort out the various possible parameter inputs and return a config
        object (dict)

        We have multiple input formats:

        1) a list, tuple, or numpy.ndarray, containing the linear parameters
        in the following order:
        * for single term: rho0, m1, tau1, c1
        * for multiple termss: rho0, m1, m2, ..., tau1, tau2, ..., c1, c2, ...

        2) a dictionary with the entries "rho0", "m", "tau", "c"

        2b) if the dictionary entries for "m", "tau", and "c" are lists, the
        entries correspond to mulitple polarisazion terms

        """
        nr_f = self.f.size

        # sort out parameters
        rho0, m, tau, c = self._sort_parameters(parameters)

        newsize = (nr_f, len(m))
        # rho0_resized = np.resize(rho0, newsize)
        m_resized = np.resize(m, newsize)
        tau_resized = np.resize(tau, newsize)
        c_resized = np.resize(c, newsize)

        omega = np.atleast_2d(2 * np.pi * self.f).T
        self.w = np.resize(omega, (len(m), nr_f)).T
        self.rho0 = rho0
        self.m = m_resized
        self.tau = tau_resized
        self.c = c_resized

        # compute some common terms
        self.otc = (self.w * self.tau) ** self.c
        self.otc2 = (self.w * self.tau) ** (2 * self.c)
        self.ang = self.c * np.pi / 2.0  # rad
        self.denom = 1 + 2 * self.otc * np.cos(self.ang) + self.otc2


class cc(cc_base):

    def response(self, parameters):
        r"""Complex response of the Cole-Cole model::
        :math:`\hat{\rho} = \rho_0 \left(1 - \sum_i m_i (1 - \frac{1}{1 + (j
        \omega \tau_i)^c_i})\right)`

        Parameters
        ----------
        parameters: list or tuple or numpy.ndarray
            Cole-Cole model parameters: rho0, m, tau, c (all linear)

        Returns
        -------
        response: :class:`sip_models.sip_response.sip_response`
            model response object
        """
        # get a config object
        self._set_parameters(parameters)
        terms = self.m * (1 - (1 / (1 + (1j * self.w * self.tau) ** self.c)))
        # sum up terms
        specs = np.sum(terms, axis=1)
        rcomplex = self.rho0 * (1 - specs)
        response = sip_response.sip_response(self.f, rcomplex=rcomplex)

        return response

    def dre_drho0(self, pars):
        r""" Compute partial derivative of real parts with respect to
        :math:`\rho_0`

        :math:`\frac{\partial \hat{\rho'}(\omega)}{\partial \rho_0} = 1 -
        \frac{m (\omega \tau)^c cos(\frac{c \pi}{2}) + (\omega \tau)^c}{1 + 2
        (\omega \tau)^c cos(\frac{c \pi}{2}) + (\omega \tau)^{2 c}}`

        Note that partial derivatives towards :math:`\rho_0` are 1D, in
        contrast to the other parameter derivatives, which usually return 2D
        arrays!

        Returns
        -------
        dre_drho0: :class:`numpy.ndarray`
            Size N (nr of frequencies) array with the derivatives

        """
        self._set_parameters(pars)
        numerator = self.m * self.otc * (np.cos(self.ang) + self.otc)
        term = numerator / self.denom
        specs = np.sum(term, axis=1)

        result = 1 - specs
        return result

    def dre_dlog10rho0(self, pars):
        """Compute partial derivative of real parts to log10(rho0)
        """

        # first call the linear response to set the parameters
        linear_response = self.dre_drho0(pars)
        result = np.log(10) * self.rho0 * linear_response
        return result

    def dre_dm(self, pars):
        r"""
        :math:`\frac{\partial \hat{\rho'}(\omega)}{\partial m} = - \rho_0 m
        (\omega \tau)^c \frac{(cos(\frac{c \pi}{2}) + (\omega \tau)^c)}{1 + 2
        (\omega \tau)^c cos(\frac{c \pi}{2}) + (\omega \tau)^{2 c}}`
        """
        self._set_parameters(pars)
        numerator = -self.otc * (np.cos(self.ang) + self.otc)
        result = numerator / self.denom
        result *= self.rho0
        return result

    def dre_dlog10m(self, pars):
        lin_response = self.dre_dm(pars)
        result = np.log(10) * self.m * lin_response
        return result

    def dre_dtau(self, pars):
        r"""
        :math:`\frac{\partial \hat{\rho'}(\omega)}{\partial \tau} = \rho_0
        \frac{-m \omega^c c \tau^{c-1} cos(\frac{c \pi}{2} - m \omega^{2 c} 2 c
        \tau^{2c - 1}}{1 + 2 (\omega \tau)^c cos(\frac{c \pi}{2}) + (\omega
        \tau)^{2 c}} +
        \rho_0 \frac{\left[m (\omega \tau)^c (cos(\frac{c \pi}{2}) + (\omega
        \tau)^c) \right] \cdot \left[ 2 \omega^c c \tau^{c-1} cos(\frac{c
        \pi}{2}) + 2 c \omega^{2 c} \tau^{2 c - 1}\right]}{\left[1 + 2 (\omega
        \tau)^c cos(\frac{c \pi}{2}) + (\omega \tau)^{2 c}\right]^2}`
        """
        self._set_parameters(pars)
        # term1
        nom1 = - self.m * self.c * self.w ** self.c * self.tau ** \
            (self.c - 1) *\
            np.cos(self.ang) - self.m * self.w ** (2 * self.c) *\
            2 * self.c * self.tau ** (2 * self.c - 1)
        term1 = nom1 / self.denom

        # term2
        nom2 = self.m * self.otc * (np.cos(self.ang) + self.otc) *\
            (2 * self.w ** self.c * self.c * self.tau ** (self.c - 1) *
                np.cos(self.ang) + 2 * self.c * self.w ** (2 * self.c) *
                self.tau ** (2 * self.c - 1))
        term2 = nom2 / self.denom ** 2

        result = term1 + term2
        result *= self.rho0
        return result

    def dre_dlog10tau(self, pars):
        lin_response = self.dre_dtau(pars)
        result = np.log(10) * self.tau * lin_response
        return result

    def dre_dc(self, pars):
        r"""
        :math:`\frac{\partial \hat{\rho'}(\omega)}{\partial c} = \rho_0
        \frac{-m ln(\omega \tau) (\omega \tau)^c cos(\frac{c \pi}{2}) + m
        (\omega\tau)^c \frac{\pi}{2} sin(\frac{c \pi}{2}) + ln(\omega
        \tau)(\omega \tau)^c}{1 + 2 (\omega \tau)^c cos(\frac{c \pi}{2}) +
        (\omega \tau)^{2 c}} +
        \rho_0 \frac{\left[-m (\omega \tau)^c (cos(\frac{c \pi}{2}) + (\omega
        \tau)^c) \right] \cdot \left[ -2 ln(\omega \tau) (\omega \tau)^c
        cos(\frac{c \pi}{2}) + 2 (\omega \tau)^c \frac{\pi}{2} cos(\frac{c
        \pi}{2} + 2 ln(\omega \tau) (\omega \tau)^{2 c}\right]}{\left[1 + 2
        (\omega \tau)^c cos(\frac{c \pi}{2}) + (\omega \tau)^{2 c}\right]^2}`
        """
        self._set_parameters(pars)
        # term1
        nom1 = - self.m * np.log(self.w * self.tau) * self.otc *\
            np.cos(self.ang) +\
            self.m * self.otc * (np.pi / 2.0) *\
            np.sin(self.ang) -\
            2 * self.m * np.log(self.w * self.tau) *\
            self.otc2
        term1 = nom1 / self.denom

        # term2
        nom2 = (self.m * self.otc * (np.cos(self.ang) + self.otc)) *\
            (2 * np.log(self.w * self.tau) * self.otc * np.cos(self.ang) -
             2 * self.otc * (np.pi / 2.0) * np.sin(self.ang) +
             2 * np.log(self.w * self.tau) * self.otc2)
        term2 = nom2 / self.denom ** 2

        result = term1 + term2
        result *= self.rho0
        return result

    def dim_drho0(self, pars):
        r"""
        :math:`\frac{\partial \hat{\rho}''(\omega)}{\partial \rho_0} = -
        \frac{m (\omega \tau)^c sin(\frac{c \pi}{2})}{1 + 2
        (\omega \tau)^c cos(\frac{c \pi}{2}) + (\omega \tau)^{2 c}}`
        """
        self._set_parameters(pars)

        result = np.sum(- self.m * self.otc * np.sin(self.ang) / self.denom,
                        axis=1)

        return result

    def dim_dlog10rho0(self, pars):
        lin_resp = self.dim_drho0(pars)
        result = np.log(10) * self.rho0 * lin_resp
        return result

    def dim_dm(self, pars):
        r"""
        :math:`\frac{\partial \hat{\rho''}(\omega)}{\partial m} = - \rho_0 m
        (\omega \tau)^c \frac{sin(\frac{c \pi}{2})}{1 + 2 (\omega \tau)^c
        cos(\frac{c \pi}{2}) + (\omega \tau)^{2 c}}`
        """
        self._set_parameters(pars)
        numerator = -self.otc * np.sin(self.ang)
        result = numerator / self.denom
        result *= self.rho0
        return result

    def dim_dlog10m(self, pars):
        lin_response = self.dim_dm(pars)
        result = np.log(10) * self.m * lin_response
        return result

    def dim_dtau(self, pars):
        r"""
        :math:`\frac{\partial \hat{\rho''}(\omega)}{\partial \tau} = \rho_0
        \frac{-m \omega^c c \tau^{c-1} sin(\frac{c \pi}{2} }{1 + 2 (\omega
        \tau)^c cos(\frac{c \pi}{2}) + (\omega \tau)^{2 c}} +
        \rho_0 \frac{\left[-m (\omega \tau)^c sin(\frac{c \pi}{2}
        \right] \cdot \left[ 2 \omega^c c \tau^{c-1} cos(\frac{c
        \pi}{2}) + 2 c \omega^{2 c} \tau^{2 c - 1}\right]}{\left[1 + 2 (\omega
        \tau)^c cos(\frac{c \pi}{2}) + (\omega \tau)^{2 c}\right]^2}`
        """
        self._set_parameters(pars)
        # term1
        nom1 = - self.m * np.sin(self.ang) * self.w ** self.c *\
            self.c * self.tau ** (self.c - 1)
        term1 = nom1 / self.denom

        # term2
        nom2 = (self.m * self.otc * np.sin(self.ang)) *\
            (2 * self.w ** self.c * self.c * self.tau ** (self.c - 1) *
             np.cos(self.ang) + 2 * self.c * self.w ** (2 * self.c) *
             self.tau ** (2 * self.c - 1))
        term2 = nom2 / self.denom ** 2

        result = term1 + term2
        result *= self.rho0
        return result

    def dim_dlog10tau(self, pars):
        lin_resp = self.dim_dtau(pars)
        result = np.log(10) * self.tau * lin_resp
        return result

    def dim_dc(self, pars):
        r"""
        :math:`\frac{\partial \hat{\rho''}(\omega)}{\partial c} = \rho_0
        \frac{-m sin(\frac{c \pi}{2}) ln(\omega \tau)(\omega \tau)^c - m
        (\omega \tau)^c \frac{\pi}{2} cos(\frac{\pi}{2}}{1 + 2 (\omega \tau)^c
        cos(\frac{c \pi}{2}) + (\omega \tau)^{2 c}} + \rho_0 \frac{\left[-m
        (\omega \tau)^c cos(\frac{c \pi}{2}) \right] \cdot \left[ -2 ln(\omega
        \tau) (\omega \tau)^c cos(\frac{c \pi}{2}) + 2 (\omega \tau)^c
        \frac{\pi}{2} cos(\frac{c \pi}{2}) \right] + \left[2 ln(\omega \tau)
        (\omega \tau)^{2 c}\right]}{\left[1 + 2 (\omega \tau)^c cos(\frac{c
        \pi}{2}) + (\omega \tau)^{2 c}\right]^2}`
        """
        self._set_parameters(pars)
        # term1
        nom1a = - self.m * np.log(self.w * self.tau) * self.otc *\
            np.sin(self.ang)
        nom1b = - self.m * self.otc * (np.pi / 2.0) * np.cos(self.ang)
        term1 = (nom1a + nom1b) / self.denom

        # term2
        nom2 = (self.m * self.otc * np.sin(self.ang)) *\
            (2 * np.log(self.w * self.tau) * self.otc * np.cos(self.ang) -
             2 * self.otc * (np.pi / 2.0) * np.sin(self.ang) +
             2 * np.log(self.w * self.tau) * self.otc2)
        term2 = nom2 / self.denom ** 2
        result = term1 + term2

        result *= self.rho0
        return result

    def Jacobian_re_im(self, pars):
        r"""
        :math:`J`

        >>> import sip_models.res.cc as cc
        >>> import numpy as np
        >>> f = np.logspace(-3, 3, 20)
        >>> pars = [100, 0.1, 0.04, 0.8]
        >>> obj = cc.cc(f)
        >>> J = obj.Jacobian_re_im(pars)
        """
        partials = []

        # partials.append(self.dre_dlog10rho0(pars)[:, np.newaxis, :])
        partials.append(self.dre_drho0(pars)[:, np.newaxis])
        partials.append(self.dre_dm(pars))
        # partials.append(self.dre_dlog10tau(pars))
        partials.append(self.dre_dtau(pars))
        partials.append(self.dre_dc(pars))
        # partials.append(self.dim_dlog10rho0(pars)[:, np.newaxis, :])
        partials.append(self.dim_drho0(pars)[:, np.newaxis])
        partials.append(self.dim_dm(pars))
        # partials.append(self.dim_dlog10tau(pars))
        partials.append(self.dim_dtau(pars))
        partials.append(self.dim_dc(pars))

        print('SHAPES')
        for x in partials:
            print(x.shape)

        J = np.concatenate(partials, axis=1)
        return J
