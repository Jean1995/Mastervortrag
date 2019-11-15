#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("Matplotlib not installed!")

try:
    import numpy as np
except ImportError:
    raise ImportError("Numpy not installed!")

try:
    import pyPROPOSAL as pp
except ImportError:
    raise ImportError("pyPROPOSAL not installed!")

from timeit import default_timer as timer

def muons(energy, statistics, vcut, do_continuous_randomization, dist):

    sec_def = pp.SectorDefinition()
    sec_def.medium = pp.medium.StandardRock(1.0)
    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), 1e20, 0)
    sec_def.particle_location = pp.ParticleLocation.inside_detector

    sec_def.scattering_model = pp.scattering.ScatteringModel.Moliere
    sec_def.do_continuous_randomization = do_continuous_randomization

    sec_def.cut_settings.ecut = 0
    sec_def.cut_settings.vcut = vcut

    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = ""

    prop = pp.Propagator(
            particle_def=pp.particle.MuMinusDef.get(),
            sector_defs=[sec_def],
            detector=pp.geometry.Sphere(pp.Vector3D(), 1e20, 0),
            interpolation_def=interpolation_def
    )

    mu = prop.particle

    mu_energies = []

    for i in range(statistics):

        mu.position = pp.Vector3D(0, 0, 0)
        mu.direction = pp.Vector3D(0, 0, -1)
        mu.energy = energy
        mu.propagated_distance = 0

        d = prop.propagate(dist * 100)

        mu_energies.append(mu.energy)

    return mu_energies


if __name__ == "__main__":

    # =========================================================
    # 	Save energies
    # =========================================================

    energy = 1e8
    statistics = int(5e4)
    dist = 300
    binning = 100


    start = timer()
    energies_1 = muons(energy, statistics, 0.05, False, dist)
    end1 = timer()
    energies_2 = muons(energy, statistics, 0.0001, False, dist)
    end2 = timer()
    energies_3 = muons(energy, statistics, 0.05, True, dist)
    end3 = timer()

    time1 = statistics / (end1-start)
    time2 = statistics / (end2 - end1)
    time3 = statistics / (end3 - end2)
    print("0.05, no rand: " + str(time1)) # all in per second
    print("0.0001, no rand: " + str(time2))
    print("0.05, rand: " + str(time3))

    tex_preamble = [
        r"\usepackage{amsmath}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
    ]

    font_size = 10

    params = {
        'backend': 'pdf',
        'font.family': 'serif',
        'font.size': 12,
        'text.usetex': True,
        'text.latex.preamble': tex_preamble,
        'axes.labelsize': font_size,
        'legend.numpoints': 1,
        'legend.shadow': False,
        'legend.fontsize': font_size,
        'xtick.labelsize': font_size,
        'ytick.labelsize': font_size,
        'axes.unicode_minus': True
    }

    plt.rcParams.update(params)

    # =========================================================
    # 	Plot energies
    # =========================================================

    binning = 100

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'Finale Myonenergie / MeV')
    ax.set_ylabel(r'Anzahl')

    fig.tight_layout()


    fig.savefig("cont_rand_0.pdf")


    ax.hist(
        energies_1,
        histtype="step",
        log=True,
        bins=binning,
        label=r"$v_\text{{cut}} = 0.05, \,  \text{{Teilchen pro Sekunde: }} {}$".format(int(round(time1)))
    )

    ax.legend(loc='upper left')
    fig.savefig("cont_rand_1.pdf")


    ax.hist(
        energies_2,
        histtype="step",
        log=True,
        bins=binning,
        label=r"$v_\text{{cut}} = 10^{{-4}}, \,  \text{{Teilchen pro Sekunde: }} {}$".format(int(round(time2)))
    )

    ax.legend(loc='upper left')
    fig.savefig("cont_rand_2.pdf")

    ax.hist(
        energies_3,
        histtype="step",
        log=True,
        bins=binning,
        label=r"$v_\text{{cut}} = 0.05, \text{{ mit kont.,}} \,  \text{{Teilchen pro Sekunde: }} {}$".format(int(round(time3)))
    )

    ax.legend(loc='upper left')
    fig.savefig("cont_rand_3.pdf")



