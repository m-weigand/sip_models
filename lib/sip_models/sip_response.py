# *-* coding: utf-8 *-*
"""
Define a container for one SIP spectrum. Include converters and plot functions

"""
import sip_formats.convert as SC
import numpy as np
import crtomo.mpl
plt, mpl = crtomo.mpl.setup()


class sip_response():
    """
    Hold one spectrum and return it in various formats

    """
    def __init__(self, frequencies, rcomplex=None, ccomplex=None):
        if rcomplex is None and ccomplex is None:
            raise Exception('One initialization array is allowed!')
        if rcomplex is not None and ccomplex is not None:
            raise Exception('Only one initialization array is allowed!')

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

    def _add_labels(self, axes):
        """Given a 2x2 array of axes, add x and y labels
        """
        ax = axes[0, 0]
        ax.set_ylabel(r'$|\rho|~[\Omega m]$')

        ax = axes[0, 1]
        ax.set_ylabel(r'$-\phi~[mrad]$')

        ax = axes[1, 0]
        ax.set_xlabel('frequency [Hz]')
        ax.set_ylabel(r"$\sigma'~[S/m]$")

        ax = axes[1, 1]
        ax.set_xlabel('frequency [Hz]')
        ax.set_ylabel(r"$\sigma''~[S/m]$")

    def _plot(self, title=None, reciprocal=None, limits=None):
        """Standard plot of spectrum

        Parameters
        ----------
        title: string|None, optional
            Title of plot
        reciprocal: sip_response object|None, optional
            If provided, plot this spectrum with another color
        limits: dict|None, optional
            used to set ylimits of the plots. Possible entries: rmag_min,
            rmag_max, rpha_min, rpha_max, cre_min, cre_max, cim_min, cim_max


        Returns
        -------
        fig: figure object
            the generated matplotlib figure
        axes: list
            matplotlib axes objects
        """
        if limits is None:
            limits = {}

        fig, axes = plt.subplots(
            2, 2, figsize=(10 / 2.54, 6 / 2.54), sharex=True
        )
        if title is not None:
            fig.suptitle(title)

        # resistivity magnitude
        ax = axes[0, 0]
        ax.semilogx(
            self.frequencies, self.rmag, '.-', color='k',
            label='normal',
        )
        ax.set_ylabel(r'$|\rho|~[\Omega m]$')
        ax.set_ylim(
            limits.get('rmag_min', None),
            limits.get('rmag_max', None)
        )

        # resistivity phase
        ax = axes[0, 1]
        ax.semilogx(self.frequencies, -self.rpha, '.-', color='k')
        ax.set_ylabel(r'$-\phi~[mrad]$')
        # note the switch of _min/_max because we change the sign while
        # plotting
        ymin = limits.get('rpha_max', None)
        if ymin is not None:
            ymin *= -1
        ymax = limits.get('rpha_min', None)
        if ymax is not None:
            ymax *= -1
        ax.set_ylim(
            ymin,
            ymax,
        )

        # conductivity real part
        ax = axes[1, 0]
        ax.loglog(self.frequencies, self.cre, '.-', color='k')
        ax.set_ylabel(r"$\sigma'~[S/m]$")
        ax.set_ylim(
            limits.get('cre_min', None),
            limits.get('cre_max', None)
        )

        # conductivity imaginary part
        ax = axes[1, 1]
        ax.loglog(self.frequencies, self.cim, '.-', color='k')
        ax.set_ylabel(r"$\sigma''~[S/m]$")
        ax.set_ylim(
            limits.get('cim_min', None),
            limits.get('cim_max', None)
        )

        fig.tight_layout()
        if reciprocal is not None:
            axes[0, 0].semilogx(
                reciprocal.frequencies,
                reciprocal.rmag,
                '.-',
                color='k',
                linestyle='dashed',
                label='reciprocal',
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

            fig.subplots_adjust(
                bottom=0.3,
            )

            axes[0, 0].legend(
                loc="lower center",
                ncol=4,
                bbox_to_anchor=(0, 0, 1, 1),
                bbox_transform=fig.transFigure
            )

        for ax in axes.flatten()[0:2]:
            ax.xaxis.set_major_locator(mpl.ticker.LogLocator(numticks=5))
            ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))

        for ax in axes.flatten()[2:]:
            ax.set_xlabel('frequency [Hz]')
            ax.xaxis.set_major_locator(mpl.ticker.LogLocator(numticks=5))
            ax.yaxis.set_major_locator(mpl.ticker.LogLocator(numticks=5))

        return fig, axes

    def plot(self, filename, title=None, reciprocal=None, limits=None):
        """Standard plot of spectrum
        """
        fig, axes = self._plot(
            reciprocal=reciprocal,
            limits=limits,
            title=title,
        )
        fig.savefig(filename, dpi=300)
        plt.close(fig)
