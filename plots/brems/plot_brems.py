
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


def plot_brems(medium, plotname, energy):    

    eminus = pp.particle.EMinusDef.get()
    cuts = pp.EnergyCutSettings(-1, -1)  # ecut, vcut
    cuts_total = pp.EnergyCutSettings(1, -1)  # ecut, vcut


    dEdx_photo = []
    dNdx_photo = []

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
            eminus,
            medium,
            cuts,
            1.0,
            False,
        ]

    param_defs_total = [
            eminus,
            medium,
            cuts_total,
            1.0,
            False,
        ]

    # =========================================================
    #   Selection of cross sections
    # ========================================================= 

    params = []
    params_total = []

    def wrapper(func):
        params.append(func(*param_defs))
        params_total.append(func(*param_defs_total))

    wrapper(parametrization.bremsstrahlung.ElectronScreening)
    wrapper(parametrization.bremsstrahlung.CompleteScreening)
    #wrapper(parametrization.bremsstrahlung.KelnerKokoulinPetrukhin)
    #wrapper(parametrization.bremsstrahlung.PetrukhinShestakov)
    wrapper(parametrization.bremsstrahlung.AndreevBezrukovBugaev)
    wrapper(parametrization.bremsstrahlung.SandrockSoedingreksoRhode)

    # =========================================================
    # 	Create x sections out of their parametrizations
    # =========================================================

    crosssections = []
    crosssections_total = []

    for param in params:
        crosssections.append(pp.crosssection.BremsIntegral(
            param,
        ))

    for param in params_total:
     	crosssections_total.append(pp.crosssection.BremsIntegral(
            param,
        ))

    # print(crosssections[0])

    # =========================================================
    # 	Calculate DE/dx at the given energies
    # =========================================================

    for cross, cross_total in zip(crosssections, crosssections_total):
        dEdx = []
        dNdx = []
        for E in energy:
            dEdx.append(cross.calculate_dEdx(E)/E)
            dNdx.append(cross_total.calculate_dNdx(E))

        dEdx_photo.append(dEdx)
        dNdx_photo.append(dNdx)
        # print(dEdx)

    # =========================================================
    # 	Plot dEdx
    # =========================================================

    colors = np.array(['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'])

    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 1])

    ax = fig.add_subplot(gs[0])

    for dEdx, param in zip(dEdx_photo, params):
        ax.loglog(
            energy,
            dEdx,
            linestyle='-',
            label="".join([param.name[5:]])
        )

    ax.axvline(50, color='r', lw=0.5, ls='-.')
    #ax.axvline(0.5, color='b', lw=0.5, ls='-.')
    ax.set_ylabel(r'-1/E dEdx / $\rm{g}^{-1} \rm{cm}^2$')

    ax.legend(loc='best')

    # ====[ ratio ]============================================

    ax = fig.add_subplot(gs[1], sharex=ax)

    start = 0
    for i in np.linspace(1, len(crosssections)-1, len(crosssections)-1, dtype=int):
    	ax.semilogx(
    	    energy[start:],
    	    np.array(dEdx_photo[i][start:]) / np.array(dEdx_photo[0][start:]),
    	    linestyle='-',
    	    color=colors[i],
    	    label="{} / {}".format(
    	        "".join([c for c in params[i].name[1:] if c.isupper()]),
    	        "".join([c for c in params[0].name[1:] if c.isupper()])
    	    )
    	)

    ax.xaxis.grid(which='major', ls=":")
    ax.yaxis.grid(which='minor', ls=":")
    #ax.set_ylim(top=1.1, bottom=0.7)
    #ax.set_xlim(left=1e1, right=1e9)

    ax.set_xlabel(r'$E$ / MeV')
    ax.set_ylabel(r'ratio')

    ax.axhline(1.0, color='k', lw=0.5, ls='-.')

    ax.axvline(50, color='r', lw=0.5, ls='-.')

    ax.legend(loc='best')
    fig.tight_layout()
    fig.savefig('brems_dEdx_' + plotname + '.pdf')
    plt.show()

    plt.clf()

    # =========================================================
    # 	Plot dNdx
    # =========================================================

    colors = np.array(['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'])

    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 1])

    ax = fig.add_subplot(gs[0])

    for dNdx, param in zip(dNdx_photo, params_total):
        ax.loglog(
            energy,
            dNdx,
            linestyle='-',
            label="".join([param.name[5:]])
        )

    ax.axvline(50, color='r', lw=0.5, ls='-.')
    ax.set_ylabel(r'dNdx / $\rm{g}^{-1} \rm{cm}^2$')

    ax.legend(loc='best')

    # ====[ ratio ]============================================

    ax = fig.add_subplot(gs[1], sharex=ax)

    start = 0
    for i in np.linspace(1, len(crosssections_total)-1, len(crosssections_total)-1, dtype=int):
    	ax.semilogx(
    	    energy[start:],
    	    np.array(dNdx_photo[i][start:]) / np.array(dNdx_photo[0][start:]),
    	    linestyle='-',
    	    color=colors[i],
    	    label="{} / {}".format(
    	        "".join([c for c in params[i].name[1:] if c.isupper()]),
    	        "".join([c for c in params[0].name[1:] if c.isupper()])
    	    )
    	)

    ax.xaxis.grid(which='major', ls=":")
    ax.yaxis.grid(which='minor', ls=":")
    #ax.set_ylim(top=1.1, bottom=0.7)
    #ax.set_xlim(left=1e1, right=1e9)

    ax.set_xlabel(r'$E$ / MeV')
    ax.set_ylabel(r'ratio')

    ax.axhline(1.0, color='k', lw=0.5, ls='-.')

    ax.axvline(50, color='r', lw=0.5, ls='-.')

    ax.legend(loc='best')
    fig.tight_layout()
    fig.savefig('brems_dNdx_' + plotname + '.pdf')
    plt.show()

if __name__ == "__main__":
    energy = np.logspace(0, 8, 100)
    energy_high = np.logspace(3, 8, 100)
    medium = pp.medium.StandardRock(1.0)
    medium2 = pp.medium.Air(1.0)
    medium3 = pp.medium.Iron(1.0)
    medium4 = pp.medium.Hydrogen(1.0)
    plot_brems(medium, "sr", energy)
    plot_brems(medium2, "air", energy)
    plot_brems(medium3, "fe", energy)
    plot_brems(medium4, "h", energy)
    plot_brems(medium, "sr_high", energy_high)


