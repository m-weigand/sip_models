Introduction
============

Provide Python classes for the Cole-Cole model, both in resistivity and
conductivity formulations.

This package is licenced under the GPL-3 licence.

Status
======

* implemented: sip_models.res.cc

  This is the resistivity formulation of the Cole-Cole model, including
  derivatives and Jacobian matrix

* implemented: sip_models.sip_response

  Hold one spectrum and return it in various formats


Roadmap
=======

* implement conductivity models
* unit tests
* Planned is also the implementation of derived models such as the Debye
  decomposition or the Warburg model.

Usage
=====

::

    import numpy as np
    frequencies = np.logspace(-3, 4, 30)

    import sip_models.res.cc as cc_res
    cc_obj = cc_res.cc(frequencies)

    cc_paramers = {
        'rho0': 100,
        'm': 0.1,
        'tau': 0.04,
        'c': 0.6
    }

    # return an SIP spectrum object
    response = cc_obj.response(cc_parameters)

    # resistivity magnitude
    response.rmag

    # resistivity phase
    response.rpha

    # resistivity real part
    response.rre

    # resistivity imaginary part
    response.rim

    # resistivity magnitude and phase
    response.rmag_rpha

    # for other formats, please look at sip_models.sip_response.py


Planned Usage
=============

These snippets represent the planned usage of this library. When its
finnished...

::

    # resistivity formulation
    import sip_models.res.ccd
    import sip_models.res.cc
    import sip_models.res.dd
    import sip_models.res.warburg
    import sip_models.res.ccd_em
    import sip_models.res.dd_em
    import sip_models.res.dd_spline
    import sip_models.res.ccd_spline

    # conductivity formulation
    import sip_models.cond.cc
    import sip_models.cond.ccd
    import sip_models.cond.dd
    import sip_models.cond.warburg
    import sip_models.cond.ccd_em
    import sip_models.cond.dd_em
    import sip_models.cond.dd_spline
    import sip_models.cond.ccd_spline


    frequencies = np.logspace(-4, 4, 30)
    parameters = [
        100,  # rho0
        0.1,  # m
        0.04, # tau
        0.5,  # c
    ]

    # also possible:
    parameters = {
        'rho0': 100,
        'm': 0.1,
        'tau': 0.04,
        'c': 0.5,
    }

    # this also holds true for multiple pol. terms
    parameters = {
        'rho0': 100,
        'm': [0.1, 0.2, 0.3],
        'tau': [0.4, 0.04, 0.004],
        'c': [0.4, 0.6, 0.8],
    }

    response = sip_models.res.cc.forward(frequencies, parameters)

    # response is a numpy.ndarray...
    print response

    # with added functionality
    print response.rmag
    print response.rpha
    print response.cre
    print response.cim
    print response.frequencies

    Jacobian = sip_models.res.cc.Jacobian(frequencies, parameters)

    # or individual derivatives
    dcre_drho0 = sip_models.res.cc.dcre_drho0(frequencies, parameters)

    # derivatives can also be found in this dict
    sip_models.res.cc.derivatives

