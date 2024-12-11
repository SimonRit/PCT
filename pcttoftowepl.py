#!/usr/bin/env python
import argparse
import json
import sys

import matplotlib.pyplot as plt
import opengate as gate
import uproot
import numpy as np


def tof_to_wepl_mc(
    output='pcttoftowepl',
    number_of_particles=1e4,
    visu=False,
    verbose=False
):

    phantom_length_cm = 20
    phantom_width_cm = 40
    source_energy_mev = 200
    detector_interval_cm = 2

    # Units
    nm = gate.g4_units.nm
    mm = gate.g4_units.mm
    cm = gate.g4_units.cm
    m = gate.g4_units.m
    sec = gate.g4_units.second
    MeV = gate.g4_units.MeV

    # Simulation
    sim = gate.Simulation()

    sim.random_engine = 'MersenneTwister'
    sim.random_seed = 'auto'
    sim.run_timing_intervals = [[0 * sec, 1 * sec]]
    sim.check_volumes_overlap = False
    sim.visu = visu
    sim.visu_type = 'vrml'
    sim.g4_verbose = False
    sim.progress_bar = verbose
    sim.number_of_threads = 1

    # Misc
    yellow = [1, 1, 0, 1]
    blue = [0, 0, 1, 1]

    # Geometry
    sim.volume_manager.add_material_database(gate.utility.get_contrib_path() / 'GateMaterials.db')
    sim.world.material = 'Vacuum'
    sim.world.size = [4 * m, 4 * m, 4 * m]
    sim.world.color = [0, 0, 0, 0]

    # Phantom
    phantom = sim.add_volume('Box', name='Phantom')
    phantom.size = [phantom_length_cm * cm, phantom_width_cm * cm, phantom_width_cm * cm]
    phantom.material = 'Water'
    phantom.color = blue

    # Beam
    source = sim.add_source('GenericSource', 'mybeam')
    source.particle = 'proton'
    source.energy.mono = source_energy_mev * MeV
    source.energy.type = 'mono'
    source.position.type = 'box'
    source.position.size = [1 * nm, 1 * nm, 1 * nm]
    source.position.translation = [((phantom_length_cm / 2) + 10) * cm, 0 * mm, 0 * mm]
    source.direction.type = 'momentum'
    source.direction.momentum = [-1, 0, 0]
    source.n = number_of_particles

    # Physics list
    sim.physics_manager.physics_list_name = 'QGSP_BIC_EMZ'

    # Phase spaces

    def add_detector(name, translation):
        plane = sim.add_volume('Box', 'PlanePhaseSpace' + name)
        plane.mother = phantom.name
        plane.size = [1 * nm, phantom_width_cm * cm, phantom_width_cm * cm]
        plane.translation = translation
        plane.material = 'Vacuum'
        plane.color = yellow

        phase_space = sim.add_actor('PhaseSpaceActor', 'PhaseSpace' + name)
        phase_space.attached_to = plane.name
        phase_space.output_filename = f'{output}/ps{name}.root'
        phase_space.attributes = [
            'EventID',
            'TrackID',
            'Position',
            'PreGlobalTime'
        ]
        particle_filter = sim.add_filter('ParticleFilter', 'Filter' + name)
        particle_filter.particle = 'proton'

        phase_space.filters.append(particle_filter)

    for x in np.arange(phantom_length_cm // 2, -phantom_length_cm // 2 - 1, -detector_interval_cm, dtype=int):
        add_detector(str(x), [x * cm, 0 * mm, 0 * mm])

    # Particle stats
    stat = sim.add_actor('SimulationStatisticsActor', 'stat')
    stat.output_filename = f'{output}/stats.txt'

    sim.run()


def tof_to_wepl_fit(
    output='pcttoftowepl',
    poly_deg=2,
    path_type='simple',
    display=False,
    savefig=False,
    verbose=False
):
    def print_verbose(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    data = uproot.concatenate(f'{output}/*.root', library='np')
    print_verbose("Loaded", len(data['EventID']), "events")

    # Sort data if needed
    ws = data['Position_X']
    if not np.all(ws[:-1] < ws[1:]):
        print_verbose("Sorting input data…")
        index_sorted = np.argsort(-ws)
        for key in data.keys():
            data[key] = data[key][index_sorted]

    # Remove nuclear interactions
    no_nuclears = data['TrackID'] == 1
    for key in data.keys():
        data[key] = data[key][no_nuclears]
    print_verbose(len(data['EventID']), 'events remain after nuclear interaction filtering')

    tofs = []
    wepls = []

    for n in np.unique(data['EventID']):
        event_mask = data['EventID'] == n
        if np.sum(event_mask) == 0:
            continue

        us = data['Position_Y'][event_mask]
        vs = data['Position_Z'][event_mask]
        ws = data['Position_X'][event_mask]
        times = data['PreGlobalTime'][event_mask]

        tofs_event = [times[k] - times[0] for k in range(len(times))]

        if path_type == 'simple':
            # Straight line between interaction position in first plane and in plane k
            wepls_event = [
                np.sqrt((us[k] - us[0])**2 + (vs[k] - vs[0])**2 + (ws[k] - ws[0])**2)
                for k in range(len(ws))
            ]
        elif path_type == 'realistic':
            # Path length through all detectors
            wepls_event = [
                np.sum([
                    np.sqrt((us[l] - us[l - 1])**2 + (vs[l] - vs[l - 1])**2 + (ws[l] - ws[l - 1])**2)
                    for l in range(1, k)
                ])
                for k in range(len(ws))
            ]
        else:
            sys.exit(f"Invalid path time {path_type}!")

        tofs.extend(tofs_event)
        wepls.extend(wepls_event)

    p = np.polyfit(tofs, wepls, deg=poly_deg)
    print_verbose("Fitted coefficients:", p)

    if display or savefig:
        tof_fit = np.linspace(0, np.max(tofs), 100)
        wepl_fit = np.polyval(p, tof_fit)

        plt.figure()
        plt.plot(tofs, wepls, '+', label="Detected events")
        plt.plot(tof_fit, wepl_fit, label=f"Polynomial fit (degree {poly_deg})")
        plt.xlabel("TOF")
        plt.ylabel("WEPL [mm]")
        plt.legend()

        if savefig:
            plt.savefig(f'{output}/tof_to_wepl_fit.pdf')
        if display:
            plt.show()

    with open(f'{output}/tof_to_wepl_fit.json', 'w', encoding='utf-8') as f:
        json.dump(list(p), f)

    return p


def pcttoftowepl(
    output='pcttoftowepl',
    number_of_particles=1e4,
    poly_deg=2,
    path_type='simple',
    visu=False,
    display=False,
    savefig=False,
    verbose=False
):
    tof_to_wepl_mc(output, number_of_particles, visu, verbose)
    p = tof_to_wepl_fit(output, poly_deg, path_type, display, savefig, verbose)
    return p


def main():

    parser = argparse.ArgumentParser(description='Convert TOF to WEPL using a fit on Monte Carlo data')
    parser.add_argument('--output', help="Path of outputs", default='pcttoftowepl')
    parser.add_argument('-n', '--number-of-particles', help="Number of generated particles", default=1e4, type=int)
    parser.add_argument('--poly-deg', help="Degrees of polynomial fit", default=2, type=int)
    parser.add_argument('--path-type', help="How to compute proton path", choices=['simple', 'realistic'], default='simple')
    parser.add_argument('--visu', help="Visualize Monte Carlo simulation", default=False, action='store_true')
    parser.add_argument('--display', help="Display polynomial fit plot", default=False, action='store_true')
    parser.add_argument('--savefig', help="Write polynomial fit plot to disk", default=False, action='store_true')
    parser.add_argument('--verbose', '-v', help="Verbose execution", default=False, action='store_true')
    args_info = parser.parse_args()

    pcttoftowepl(**vars(args_info))

if __name__ == '__main__':
    main()
