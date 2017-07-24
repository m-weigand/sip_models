# *-* coding: utf-8 *-*
"""
Define a container for one SIP spectrum. Include converters and plot functions

"""
import sip_formats.convert as SC
import numpy as np
from crtomo.mpl_setup import *
# import matplotlib as mpl
# mpl.rcParams['font.size'] = 8.0
# import pylab as plt
# mpl.rcParams['font.size'] = 8.0


class sip_response():
    """
    Hold one spectrum and return it in various formats

    """
    def __init__(self, frequencies, rcomplex=None, ccomplex=None):
        if rcomplex is None and ccomplex is None:
            raise Exception('One initialisation array is required!')
        if rcomplex is not None and ccomplex is not None:
            raise Exception('Only one initialisation array is required!')

        self.frequencies = frequencies

        if rcomplex is not None:
            self.rcomplex = rcomplex
            self.ccomplex = SC.convert('rcomplex', 'ccomplex', rcomplex)
        elif ccomplex is not None:
            self.ccomplex = ccomplex
            self.rcomplex = SC.convert('ccomplex', 'rcomplex', ccomplex)

        self.rmag = np.abs(self.rcomplex)
        self.rpha = np.arctan2(
            np.imag(self.rcomplex),
            np.real(self.rcomplex)
        ) * 1000
        self.cmag = np.abs(self.ccomplex)
        self.cpha = np.arctan2(
            np.imag(self.ccomplex),
            np.real(self.ccomplex)
        ) * 1000

        self.rmag_rpha = np.vstack((self.rmag, self.rpha)).T
        self.cmag_cpha = np.vstack((self.cmag, self.cpha)).T

        self.rre = np.real(self.rcomplex)
        self.rim = np.imag(self.rcomplex)
        self.cre = np.real(self.ccomplex)
        self.cim = np.imag(self.ccomplex)

        self.rre_rim = np.vstack((self.rre, self.rim)).T
        self.cre_cim = np.vstack((self.cre, self.cim)).T

    def to_one_line(self, array):
        """flatten the array to one dimension using the 'F' (Fortran) style and
        return a 2D array
        """
        return np.atleast_2d(array.flatten(order='F'))

    def plot(self, filename, reciprocal=None):
        """Standard plot of spectrum
        """
        fig, axes = plt.subplots(
            2, 2, figsize=(10 / 2.54, 6 / 2.54), sharex=True
        )

        # resistivity magnitude
        ax = axes[0, 0]
        ax.semilogx(self.frequencies, self.rmag, '.-', color='k')
        ax.set_ylabel(r'$|\rho|~[\Omega m]$')

        # resistivity phase
        ax = axes[0, 1]
        ax.semilogx(self.frequencies, -self.rpha, '.-', color='k')
        ax.set_ylabel(r'$-\phi~[mrad]$')

        # conductivity real part
        ax = axes[1, 0]
        ax.loglog(self.frequencies, self.cre, '.-', color='k')
        ax.set_ylabel(r"$\sigma'~[S/m]$")

        # conductivity imaginary part
        ax = axes[1, 1]
        ax.loglog(self.frequencies, self.cim, '.-', color='k')
        ax.set_ylabel(r"$\sigma''~[S/m]$")

        if reciprocal is not None:
            axes[0, 0].semilogx(
                reciprocal.frequencies,
                reciprocal.rmag,
                '.-',
                color='k',
                linestyle='dashed',
            )
            axes[0, 1].semilogx(
                reciprocal.frequencies,
                -reciprocal.rpha,
                '.-',
                color='k',
                linestyle='dashed',
            )
            axes[1, 0].loglog(
                reciprocal.frequencies,
                reciprocal.cre,
                '.-',
                color='k',
                linestyle='dashed',
            )
            axes[1, 1].loglog(
                reciprocal.frequencies,
                reciprocal.cim,
                '.-',
                color='k',
                linestyle='dashed',
            )

        for ax in axes.flatten()[0:2]:
            ax.xaxis.set_major_locator(mpl.ticker.LogLocator(numticks=5))
            ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))

        for ax in axes.flatten()[2:]:
            ax.set_xlabel('frequency [Hz]')
            ax.xaxis.set_major_locator(mpl.ticker.LogLocator(numticks=5))
            ax.yaxis.set_major_locator(mpl.ticker.LogLocator(numticks=5))

        fig.tight_layout()
        fig.savefig(filename, dpi=300)
        plt.close(fig)
