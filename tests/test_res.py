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
    pars2 = [1000, 0.1, 0.1, 0.2]
    s['p'] = [pars, pars2]
    s['obj'] = cc.cc(s['f'])

    s['responses'] = {
        0: {
            'rre': np.array([99.99591051, 99.99266375, 99.98681164,
                             99.97620375, 99.95678877, 99.9206855,
                             99.85188258, 99.71628544, 99.43929068,
                             98.86466816, 97.73315998, 95.87692785,
                             93.69912846, 91.98080919, 90.98416593,
                             90.48756482, 90.24828323, 90.1303687,
                             90.07009753, 90.03828985]),
            'rim': np.array([-0.01253267, -0.02240756, -0.04004258,
                             -0.07149009, -0.12741901, -0.22639435,
                             -0.3998905, -0.69832466, -1.19178546,
                             -1.93962096,
                             -2.86640995, -3.55906024, -3.46922655,
                             -2.67933078, -1.76839385, -1.07303918,
                             -0.62495881, -0.35681798, -0.20170429,
                             -0.11343277]),
        },
        1: {
            'rre': np.array([981.84641815, 979.54729305, 977.03507033,
                             974.31005743, 971.37783503, 968.2499524,
                             964.94435731, 961.48545676, 957.90372902,
                             954.23485241, 950.51837684, 946.79603085,
                             943.1098148, 939.50006701, 936.00369406,
                             932.65272892, 929.47332703, 926.48524501,
                             923.70178301, 921.13012163]),
            'rim': np.array([-4.75423444, -5.19849775, -5.64464455,
                             -6.0825728, -6.50076027, -6.88671331,
                             -7.22760444,
                             -7.51106524, -7.72606711, -7.86379349,
                             -7.91839164, -7.88749702, -7.77245071,
                             -7.5781757, -7.31273163, -6.98661624,
                             -6.61191526,
                             -6.20141333, -5.76776854, -5.32282679])
        },
    }
    return s


def test_make_list(setup):
    assert isinstance(cc._make_list(4), list)
    assert isinstance(cc._make_list([4, ]), list)


def test_sort_parameters(setup):
    """this is just a simple test to make sure the different inputs don't lead
    to crashes"""
    obj = setup['obj']
    obj._sort_parameters((100, 0.1, 0.04, 0.5))
    obj._sort_parameters((100, 0.1, 0.2, 0.04, 0.004, 0.5, 0.6))

    obj._sort_parameters(
        {'rho0': 100,
         'm': 0.1,
         'tau': 0.04,
         'c': 0.6
         }
    )

    obj._sort_parameters(
        {'rho0': 100,
         'm': (0.1, 0.2),
         'tau': (0.04, 0.004),
         'c': (0.6, 0.8),
         }
    )

    with pytest.raises(Exception):
        obj._sort_parameters(None)


def test_response(setup):
    obj = setup['obj']
    for nr, pars in enumerate(setup['p']):
        response = obj.response(pars)
        # check responses
        np.allclose(
            response.rre,
            setup['responses'][nr]['rre']
        )
        np.allclose(
            response.rim,
            setup['responses'][nr]['rim']
        )


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

        # log10rho0
        def ffunc_re_log10rho0(pars):
            """reparameterize the response for log10(rho0)"""
            pars2 = pars.copy()
            pars2[0] = 10 ** pars2[0]
            return obj.response(pars2).rre

        Jfunc_re_log10rho0 = nd.Jacobian(ffunc_re_log10rho0, order=4)
        pars_logrho = pars.copy()
        pars_logrho[0] = np.log10(pars_logrho[0])

        assert np.allclose(
            Jfunc_re_log10rho0(pars_logrho)[:, 0],
            obj.dre_dlog10rho0(pars).squeeze()
        )

        def ffunc_im_log10rho0(pars):
            """reparameterize the response for log10(rho0)"""
            pars2 = pars.copy()
            pars2[0] = 10 ** pars2[0]
            return obj.response(pars2).rim

        Jfunc_im_log10rho0 = nd.Jacobian(ffunc_im_log10rho0, order=4)
        pars_logrho = pars.copy()
        pars_logrho[0] = np.log10(pars_logrho[0])

        assert np.allclose(
            Jfunc_im_log10rho0(pars_logrho)[:, 0],
            obj.dim_dlog10rho0(pars).squeeze()
        )

        # log10m
        def ffunc_re_log10m(pars):
            """reparameterize the response for log10(m)"""
            pars2 = pars.copy()
            pars2[1] = 10 ** pars2[1]
            return obj.response(pars2).rre

        Jfunc_re_log10m = nd.Jacobian(ffunc_re_log10m, order=4)
        pars_logrho = pars.copy()
        pars_logrho[1] = np.log10(pars_logrho[1])

        assert np.allclose(
            Jfunc_re_log10m(pars_logrho)[:, 1],
            obj.dre_dlog10m(pars).squeeze()
        )

        def ffunc_im_log10m(pars):
            """reparameterize the response for log10(m)"""
            pars2 = pars.copy()
            pars2[1] = 10 ** pars2[1]
            return obj.response(pars2).rim

        Jfunc_im_log10m = nd.Jacobian(ffunc_im_log10m, order=4)
        pars_logrho = pars.copy()
        pars_logrho[1] = np.log10(pars_logrho[1])

        assert np.allclose(
            Jfunc_im_log10m(pars_logrho)[:, 1],
            obj.dim_dlog10m(pars).squeeze()
        )

        # log10tau
        def ffunc_re_log10tau(pars):
            """reparameterize the response for log10(tau)"""
            pars2 = pars.copy()
            pars2[2] = 10 ** pars2[2]
            return obj.response(pars2).rre

        Jfunc_re_log10tau = nd.Jacobian(ffunc_re_log10tau, order=4)

        def ffunc_im_log10tau(pars):
            """reparameterize the response for log10(tau)"""
            pars2 = pars.copy()
            pars2[2] = 10 ** pars2[2]
            return obj.response(pars2).rim

        Jfunc_im_log10tau = nd.Jacobian(ffunc_im_log10tau, order=4)

        pars_logtau = pars.copy()
        pars_logtau[2] = np.log10(pars_logtau[2])

        assert np.allclose(
            Jfunc_re_log10tau(pars_logtau)[:, 2],
            obj.dre_dlog10tau(pars).squeeze()
        )

        assert np.allclose(
            Jfunc_im_log10tau(pars_logtau)[:, 2],
            obj.dim_dlog10tau(pars).squeeze()
        )

        # test full Jacobian
        def ffunc_re(pars):
            return obj.response(pars).rre

        def ffunc_im(pars):
            return obj.response(pars).rim

        Jfunc_re = nd.Jacobian(ffunc_re, order=4)
        Jfunc_im = nd.Jacobian(ffunc_im, order=4)

        assert np.allclose(
            Jfunc_re(pars),
            obj.Jacobian_re_im(pars)[:, 0:4]
        )
        assert np.allclose(
            Jfunc_im(pars),
            obj.Jacobian_re_im(pars)[:, 4:8]
        )
