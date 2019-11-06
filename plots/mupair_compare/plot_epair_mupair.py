
import pyPROPOSAL as pp
import pyPROPOSAL.parametrization as parametrization

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
except ImportError:
    raise ImportError("Matplotlib not installed!")

try:
    import numpy as np
except ImportError:
    raise ImportError(
        "Numpy not installed! Needed to calculate the detector cylinder"
    )

import math


if __name__ == "__main__":

    mu = pp.particle.MuMinusDef.get()
    medium = pp.medium.StandardRock(1.0)  # With densitiy correction
    cuts = pp.EnergyCutSettings(-1, -1)  # ecut, vcut

    dEdx_photo = []
    energy = np.logspace(2, 9, 2000)

    interpolation_def = pp.InterpolationDef()

    # =========================================================
    # 	Constructor args for parametrizations
    #
    #   - particle
    #   - medium
    #   - cut
    #   - multiplier
    #   - lpm effect
    #   - interpolation definition
    # =========================================================

    param_defs = [
            mu,
            medium,
            cuts,
            1.0,
            True,
        ]

    param_defs_mupair = [
            mu,
            medium,
            cuts,
            1.0,
            False,
        ]

    params = [
        parametrization.pairproduction.KelnerKokoulinPetrukhin(
            *param_defs
        ),
        parametrization.mupairproduction.KelnerKokoulinPetrukhin(
            *param_defs_mupair
        )
    ]

    # =========================================================
    # 	Create x sections out of their parametrizations
    # =========================================================

    crosssections = []

    crosssections.append(pp.crosssection.EpairIntegral(
        params[0]
            ))

    crosssections.append(pp.crosssection.MupairIntegral(
        params[1]
    ))

    # =========================================================
    # 	Calculate DE/dx at the given energies
    # =========================================================

    for cross in crosssections:
        dEdx = []
        for E in energy:
            dEdx.append(cross.calculate_dEdx(E))

        dEdx_photo.append(dEdx)

    # =========================================================
    # 	Plot
    # =========================================================

    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[2, 1])

    ax = fig.add_subplot(gs[0])

    labels = [r'$e$-Paarproduktion', r'$\mu$-Paarproduktion']

    for dEdx, param, _label in zip(dEdx_photo, params, labels):
        ax.loglog(
            energy,
            dEdx,
            linestyle='-',
            label=_label
        )

    ax.set_ylabel(r'$\left(\frac{\mathrm{d}E}{\mathrm{d}X}\right) \,\left/\, \left( \rm{MeV} \cdot \rm{g}^{-1} \rm{cm}^2 \right) \right. $', fontsize=12)
    ax.xaxis.grid(which='major', ls="-")
    ax.yaxis.grid(which='major', ls="-")
    ax.legend(loc='best')

    # ====[ ratio ]============================================

    ax = fig.add_subplot(gs[1], sharex=ax)

    start = 0
    ax.loglog(
        energy[start:],
        np.array(dEdx_photo)[1][start:] / np.array(dEdx_photo[0][start:]),
        linestyle='-',
        label=""
    )

    ax.xaxis.grid(which='major', ls="-")
    ax.yaxis.grid(which='major', ls="-")

    ax.set_xlabel(r'$E$ / MeV', fontsize=12)
    ax.set_ylabel(r'$\left(\frac{\mathrm{d}E}{\mathrm{d}X}\right)_\mu \,\left/\, \left(\frac{\mathrm{d}E}{\mathrm{d}X}\right)_e \right.$', fontsize=12)
    plt.tight_layout()
    fig.savefig('mupair_compare.pdf')
    plt.show()
