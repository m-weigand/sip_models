# *-* coding: utf-8 *-*
"""
Define a container for one SIP spectrum. Include converters and plot functions

"""
import sip_formats.convert as SC
import numpy as np
import matplotlib as mpl
mpl.rcParams['font.size'] = 8.0
import pylab as plt


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

    def plot(self, filename):
        """Standard plot of spectrum
        """
        fig, axes = plt.subplots(2, 2, figsize=(10 / 2.54, 6 / 2.54))

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
        ax.semilogx(self.frequencies, self.cre, '.-', color='k')
        ax.set_ylabel(r"$\sigma'~[S/m]$")

        # conductivity imaginary part
        ax = axes[1, 1]
        ax.semilogx(self.frequencies, self.cim, '.-', color='k')
        ax.set_ylabel(r"$\sigma''~[S/m]$")

        for ax in axes.flatten():
            ax.set_xlabel('frequency [Hz]')
            ax.xaxis.set_major_locator(mpl.ticker.LogLocator(numticks=5))
            ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))

        fig.tight_layout()
        fig.savefig(filename, dpi=300)
        plt.close(fig)
