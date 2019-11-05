from __future__ import division

import sys
import os
import pyPROPOSAL as pp
import math
import time
import datetime

try:
    import matplotlib
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import AutoMinorLocator
    from mpl_toolkits.axes_grid1 import make_axes_locatable

except ImportError:
    print("Matplotlib not installed!")

try:
    import numpy as np
except ImportError:
    print("Numpy not installed!")

try:
    from sklearn.utils import check_random_state
except ImportError:
    print("SkLearn not installed!")


class ProgressBar(object):

    def __init__(self, loops, bar_lenght=50, start=0, **keywords):

        self._bar_lenght = bar_lenght
        self._bar = []
        self._loops = loops
        self._start = float(start)
        self._current_loop = start

        self._started_process = False
        self._start_time = None

        self._pacman = False

        self._status = ""
        self._text = "\rPercent: [{0}] {1}% Time: {2} Iteration: {3}/{4} {5}"

        self._bar_full = "="
        self._bar_empty = " "

        for key, value in keywords.iteritems():
            if key is "pacman":
                assert type(value) is bool
                self._pacman = value

        if self._pacman:
            self._bar_full = "-"
            self._bar_empty = "o"

            current = self._bar_empty
            for i in range(self._bar_lenght):
                if current is self._bar_empty:
                    current = " "
                    self._bar.append(current)
                else:
                    current = self._bar_empty
                    self._bar.append(current)
        else:
            for i in range(self._bar_lenght):
                self._bar.append(self._bar_empty)

        self._current_pac_state = "C"
        self._current_pac_block = 0

    def reset(self):
        self._current_loop = self._start
        self._status = ""
        self._started_process = False

    def start(self):
        self._started_process = True
        self._start_time = time.time()

    def update(self):
        if self._started_process is False:
            print("Pleas start ProgressBar before updating it!")
            return

        self._current_loop += 1.0
        progress = self._current_loop / self._loops

        if progress >= 1.0:
            self._status = "Done...\n"

        if self._pacman:
            block = int((self._bar_lenght - 1) * progress)

            if self._current_pac_block < block:
                self._current_pac_block = block
                if self._current_pac_state is "c":
                    self._current_pac_state = "C"
                else:
                    self._current_pac_state = "c"
            else:
                pass

            self._bar[block] = '\033[1m' + "\033[93m" + \
                               self._current_pac_state + '\033[0m'
            self._bar[:block] = block * [self._bar_full]
        else:
            block = int(self._bar_lenght * progress)
            self._bar[:block] = block * [self._bar_full]

        text = self._text.format(
            "".join(self._bar),
            progress*100,
            str(datetime.timedelta(seconds=(time.time() - self._start_time))),
            int(self._current_loop),
            self._loops,
            self._status
        )

        sys.stdout.write(text)
        sys.stdout.flush()
#

def power_law_sampler(gamma, xlow, xhig, n, random_state=None):
    r"""
    Sample n events from a power law with index gamma between xlow and xhig
    by using the analytic inversion method. The power law pdf is given by

    .. math::
       \mathrm{pdf}(\gamma) = x^{-\gamma} / \mathrm{norm}

    where norm ensures an area under curve of one. Positive spectral index
    gamma means a falling spectrum.

    Note: When :math:`\gamma=1` the integral is

    .. math::
       \int 1/x \mathrm{d}x = ln(x) + c

    This case is also handled.

    Sampling of power laws over multiple order of magnitude with the rejection
    method is VERY inefficient.

    Parameters
    ----------
    gamma : float
        Power law index.
    xlow, xhig : float
        Border of the pdf, needed for proper normalization.
    n : int
        Number of events to be sampled.
    random_state : seed, optional
        Turn seed into a np.random.RandomState instance. See
        `sklearn.utils.check_random_state`. (default: None)

    Returns
    -------
    sample : float array
        Array with the requested n numbers drawn distributed as a power law
        with the given parameters.
    """
    rndgen = check_random_state(random_state)
    # Get uniform random number which is put in the inverse function
    u = rndgen.uniform(size=int(n))

    if gamma == 1:
        return np.exp(u * np.log(xhig / xlow)) * xlow
    else:
        radicant = (u * (xhig**(1. - gamma) - xlow**(1. - gamma)) +
                    xlow**(1. - gamma))
        return radicant**(1. / (1. - gamma))
#


def propagate_muons():

    # start_time = time.time()

    mu_def = pp.particle.MuMinusDef.get()
    geometry = pp.geometry.Sphere(pp.Vector3D(), 1.e20, 0.0)
    ecut = 500
    vcut = 5e-2

    sector_def = pp.SectorDefinition()
    sector_def.cut_settings = pp.EnergyCutSettings(ecut, vcut)
    sector_def.medium = pp.medium.StandardRock(1.0)
    sector_def.geometry = geometry
    sector_def.scattering_model = pp.scattering.ScatteringModel.NoScattering
    sector_def.crosssection_defs.brems_def.lpm_effect = True
    sector_def.crosssection_defs.epair_def.lpm_effect = True

    sector_def.crosssection_defs.mupair_def.parametrization = pp.parametrization.mupairproduction.MupairParametrization.KelnerKokoulinPetrukhin
    sector_def.crosssection_defs.mupair_def.particle_output = False

    # sector_def.crosssection_defs.photo_def.parametrization = pp.PhotoParametrization.BezrukovBugaev
    # sector_def.do_stochastic_loss_weighting = True
    # sector_def.stochastic_loss_weighting = -0.1

    detector = geometry

    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = "~/.local/share/PROPOSAL/tables"

    prop = pp.Propagator(mu_def, [sector_def], detector, interpolation_def)

    statistics_log = 6
    statistics = int(10**statistics_log)
    propagation_length = 1e4 # cm
    E_min_log = 8.0
    E_max_log = 8.0
    spectral_index = 1
    pp.RandomGenerator.get().set_seed(1234)

    #muon_energies = np.logspace(E_min_log, E_max_log, statistics)
    muon_energies = power_law_sampler(spectral_index, 10**E_min_log, 10**E_max_log, statistics)
    #muon_energies = np.ones(statistics)*10**10
    epair_secondary_energy = []
    brems_secondary_energy = []
    ioniz_secondary_energy = []
    photo_secondary_energy = []
    mupair_secondary_energy = []

    plot_energy_histogram(np.log10(muon_energies), 'energyhist_{}_{}_Emin_{}_Emax_{}.pdf'.format(
                prop.particle.particle_def.name,
                sector_def.medium.name.lower(),
                E_min_log,
                E_max_log,
                ecut,
                vcut))

    progress = ProgressBar(statistics, pacman=True)
    progress.start()

    for mu_energy in muon_energies:
        progress.update()

        prop.particle.position = pp.Vector3D(0, 0, 0)
        prop.particle.direction = pp.Vector3D(0, 0, -1)
        prop.particle.propagated_distance = 0
        prop.particle.energy = mu_energy

        secondarys = prop.propagate(propagation_length)

        for sec in secondarys:
            log_sec_energy = math.log10(sec.energy)

            if sec.id == pp.particle.Data.Epair:
                epair_secondary_energy.append(log_sec_energy)
            elif sec.id == pp.particle.Data.Brems:
                brems_secondary_energy.append(log_sec_energy)
            elif sec.id == pp.particle.Data.DeltaE:
                ioniz_secondary_energy.append(log_sec_energy)
            elif sec.id == pp.particle.Data.NuclInt:
                photo_secondary_energy.append(log_sec_energy)
            elif sec.id == pp.particle.Data.MuPair:
                mupair_secondary_energy.append(log_sec_energy)
            else:
                print("Something unknown: ", sec.id)
    # =========================================================
    #   Write
    # =========================================================

    dir_prefix = ""
    np.savez(
        os.path.join(
            dir_prefix,
            'data_mupair_sec_dist_{}_{}_Emin_{}_Emax_{}'.format(
                prop.particle.particle_def.name,
                sector_def.medium.name.lower(),
                E_min_log,
                E_max_log,
                ecut,
                vcut)
        ),
        brems=brems_secondary_energy,
        epair=epair_secondary_energy,
        photo=photo_secondary_energy,
        ioniz=ioniz_secondary_energy,
        mupair=mupair_secondary_energy,
        statistics=[statistics],
        statistics_log=[statistics_log],
        E_min=[E_min_log],
        E_max=[E_max_log],
        spectral_index=[spectral_index],
        distance=[prop.particle.propagated_distance / 100],
        medium_name=[sector_def.medium.name.lower()],
        particle_name=[prop.particle.particle_def.name],
        ecut=[ecut],
        vcut=[vcut]
    )

    #statistics:
    num_all = len(brems_secondary_energy) + len(epair_secondary_energy) + len(photo_secondary_energy) + len(ioniz_secondary_energy) + len(mupair_secondary_energy)
    print(num_all)
    print("Brem: ", len(brems_secondary_energy), len(brems_secondary_energy)/num_all)
    print("Epair: ",  len(epair_secondary_energy), len(epair_secondary_energy)/num_all)
    print("photo: ", len(photo_secondary_energy), len(photo_secondary_energy)/num_all)
    print("Ioniz: ", len(ioniz_secondary_energy), len(ioniz_secondary_energy)/num_all)
    print("MuPair: ", len(mupair_secondary_energy), len(mupair_secondary_energy)/num_all)

    plot_secondary_spectrum('data_mupair_sec_dist_{}_{}_Emin_{}_Emax_{}.npz'.format(
                prop.particle.particle_def.name,
                sector_def.medium.name.lower(),
                E_min_log,
                E_max_log,
                ecut,
                vcut))
#


def plot_secondary_spectrum(data_name):

    # =========================================================
    #   Plot
    # =========================================================

    # tex_preamble = [
    #     r"\usepackage{amsmath}",
    #     r"\usepackage[utf8]{inputenc}",
    #     r"\usepackage[T1]{fontenc}",
    # ]

    # font_size = 10

    # params = {
    #     'backend': 'pdf',
    #     'font.family': 'serif',
    #     'font.size': 12,
    #     'text.usetex': True,
    #     'text.latex.preamble': tex_preamble,
    #     'axes.labelsize': font_size,
    #     'legend.numpoints': 1,
    #     'legend.shadow': False,
    #     'legend.fontsize': font_size,
    #     'xtick.labelsize': font_size,
    #     'ytick.labelsize': font_size,
    #     'axes.unicode_minus': True
    # }

    # plt.rcParams.update(params)

    inch_to_cm = 2.54
    golden_ratio = 1.61803
    width = 29.7  # cm

    # # =========================================================
    # #   All hists together $10^{{{:.0g}}}$
    # # =========================================================

    npzfile = np.load(data_name)

    ioniz_secondary_energy = npzfile['ioniz']
    brems_secondary_energy = npzfile['brems']
    photo_secondary_energy = npzfile['photo']
    epair_secondary_energy = npzfile['epair']
    mupair_secondary_energy = npzfile['mupair']

    statistics = npzfile['statistics'][0]
    statistics_log = npzfile['statistics_log'][0]
    E_min_log = npzfile['E_min'][0]
    E_max_log = npzfile['E_max'][0]
    spectral_index = npzfile['spectral_index'][0]
    distance = npzfile['distance'][0]
    medium_name = npzfile['medium_name'][0]
    particle_name = npzfile['particle_name'][0]
    ecut = npzfile['ecut'][0]
    vcut = npzfile['vcut'][0]


    fig_all = plt.figure(
        figsize=(width / inch_to_cm, width / inch_to_cm / golden_ratio)
    )
    #fig_all.suptitle(r"Propagating $10^{{{:.2g}}}$ {}, sampled from a $E^{{-{:.2g}}}$ power law with energies between $10^{{{:.2g}}}$ and $10^{{{:.2g}}}$ MeV." "\n" "Propagation length is {} m in {}".format(
    #    statistics_log,
    #    particle_name,
    #    spectral_index,
    #    E_min_log,
    #    E_max_log,
    #    distance,
    #    medium_name,
    # 
    #))

    ax_all = fig_all.add_subplot(111)
    ax_all.hist(
        [
            ioniz_secondary_energy,
            photo_secondary_energy,
            brems_secondary_energy,
            epair_secondary_energy,
            mupair_secondary_energy,
            np.concatenate((
                ioniz_secondary_energy,
                brems_secondary_energy,
                photo_secondary_energy,
                epair_secondary_energy,
                mupair_secondary_energy)
            )
        ],
        histtype='step',
        log=True,
        bins=100,
        label=['Ionisation', 'Photonuklear', 'Bremsstrahlung', r'$e$ Paarproduktion', r'$\mu$ Paarproduktion', 'Summe']
    )
    # ax_all.set_ylim(ymin=0)
    minor_locator = AutoMinorLocator()
    ax_all.xaxis.set_minor_locator(minor_locator)
    ax_all.legend(fontsize=14)
    ax_all.set_xlabel(r'$ \log\left( E \cdot v \,/\, \mathrm{MeV} \right)$', fontsize=16)
    ax_all.set_ylabel(r'$N$', fontsize=16)
    plt.xlim(left=2.5)

    fig_all.savefig("all_mupair_{}_stats_{}_Emin_{}_Emax_{}_index_{}.pdf".format(
        medium_name,
        statistics,
        E_min_log,
        E_max_log,
        spectral_index
    ))
#

def plot_energy_histogram(energies, dataname):
    plt.hist(energies)
    plt.xlabel(r'$log(E/\mathrm{MeV})$')
    plt.ylabel(r'$N$')
    plt.savefig(dataname)
    plt.clf()


if __name__ == "__main__":
    #propagate_muons()
    plot_secondary_spectrum("data_mupair_sec_dist_MuMinus_standardrock_Emin_8.0_Emax_8.0.npz")
    #plot_theory_curve()
#
