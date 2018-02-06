# test resistivity model
# *-* coding: utf-8 *-*
import pytest

import numpy as np
import numdifftools as nd

import sip_models.res.cc as cc


@pytest.fixture
def setup():
    s = {}
    s['f'] = np.logspace(-3, 3, 20)
    pars = [100, 0.1, 0.04, 0.8]
    # pars2 = [100, 0.1, 0.1, 0.5, 0.1, 0.4, 0.6]

    pars2 = [100, 0.1, 0.04, 0.8]
    s['p'] = [pars, pars2]
    s['obj'] = cc.cc(s['f'])
    return s


def test_derivates(setup):
    obj = setup['obj']
    for pars in setup['p']:
        def ffunc_re(pars):
            return obj.response(pars).rre
        def ffunc_im(pars):
            return obj.response(pars).rim

        Jfunc_re = nd.Jacobian(ffunc_re, order=4)
        Jfunc_im = nd.Jacobian(ffunc_im, order=4)

        # rho0
        assert np.allclose(Jfunc_re(pars)[:, 0], obj.dre_drho0(pars).squeeze())
        assert np.allclose(Jfunc_im(pars)[:, 0], obj.dim_drho0(pars).squeeze())

        # m
        assert np.allclose(Jfunc_re(pars)[:, 1], obj.dre_dm(pars).squeeze())
        assert np.allclose(Jfunc_im(pars)[:, 1], obj.dim_dm(pars).squeeze())

        # tau
        assert np.allclose(Jfunc_re(pars)[:, 2], obj.dre_dtau(pars).squeeze())
        assert np.allclose(Jfunc_im(pars)[:, 2], obj.dim_dtau(pars).squeeze())

        # c
        assert np.allclose(Jfunc_re(pars)[:, 3], obj.dre_dc(pars).squeeze())
        assert np.allclose(Jfunc_im(pars)[:, 3], obj.dim_dc(pars).squeeze())

        # log10tau
        def ffunc_re_log10tau(pars):
            pars2 = pars.copy()
            pars2[2] = 10 ** pars2[2]
            return obj.response(pars).rre

        Jfunc_re_log10tau = nd.Jacobian(ffunc_re_log10tau, order=4)

        pars_logtau = pars.copy()
        pars_logtau[2] = np.log10(pars_logtau[2])

        assert np.allclose(
            Jfunc_re_log10tau(pars_logtau)[:, 2],
            obj.dre_dlog10tau(pars).squeeze()
        )

