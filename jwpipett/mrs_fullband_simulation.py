#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# @author: Patrick Kavanagh (DIAS)
#
"""

Run a fullband MRS sim

python mrs_fullband.py <flux in uJy> <noise True/False>

"""
import os, glob, sys, shutil
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from mirisim import MiriSimulation
from mirisim.skysim import *
from mirisim.config_parser import *

from mirisim.skysim.sed import PLSed

class MRSFullband:
    """
    Run full band MRS
    """
    def __init__(self, exps, ints, groups, noise, dither, output_dir):
        self.exps = exps
        self.ints = ints
        self.groups = groups
        self.noise = noise
        self.dither = dither
        self.output_dir = output_dir

    def make_sim_config(self, band='SHORT'):
        mrs_sim_config = SimConfig.makeSim(name="MRS Simulation",
                                           rel_obsdate=0.0,
                                           scene="scene.ini",
                                           POP='MRS',
                                           ConfigPath='MRS_1SHORT',
                                           Dither=self.dither,
                                           StartInd=1,
                                           NDither=4,
                                           DitherPat="mrs_recommended_dither.dat",
                                           filter="F1130W",
                                           readDetect= 'FULL',
                                           ima_mode= 'FAST',
                                           ima_exposures=1,
                                           ima_integrations=1,
                                           ima_frames=50,
                                           disperser= band,
                                           detector='BOTH',
                                           mrs_mode='FAST',
                                           mrs_exposures=self.exps,
                                           mrs_integrations=self.ints,
                                           mrs_frames=self.groups)

        return mrs_sim_config

    def make_simulator_config(self):
        mrs_simulator_config = SimulatorConfig.from_default()
        #mrs_simulator_config['SCASim']['include_refpix'] = 'F'
        #mrs_simulator_config['SCASim']['include_badpix'] = 'F'
        #mrs_simulator_config['SCASim']['include_dark'] = 'F'
        #mrs_simulator_config['SCASim']['include_flat'] = 'F'
        #mrs_simulator_config['SCASim']['include_gain'] = 'F'
        #mrs_simulator_config['SCASim']['include_nonlinearity'] = 'F'
        #mrs_simulator_config['SCASim']['include_drifts'] = 'F'
        #mrs_simulator_config['SCASim']['include_latency'] = 'F'
        #mrs_simulator_config['SCASim']['cosmic_ray_mode'] = 'NONE'

        if not self.noise:
            mrs_simulator_config['SCASim']['include_poisson'] = 'F'
            mrs_simulator_config['SCASim']['include_readnoise'] = 'F'

        return mrs_simulator_config

    def make_scene_config(self):
        # set background
        background = Background(level= 'low',gradient=5.,pa=15.0,centreFOV=(0., 0.))

        # set PL SED
        SED1 = PLSed(alpha=1.0, wref=10., flux=1.e5)
        Point1 = Point(Cen=(0., 0.), vel=0.0)
        Point1.set_SED(SED1)

        Point2 = Point(Cen=(0., 0.), vel=0.0)
        SED2 = LinesSed(wavels=[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
                        fluxes=[5e+04, 4e+04, 3e+04, 2e+04, 2e+04, 1e+04, 1e+04, 9e+04, 8e+03, 7e+03, 6e+03, 5e+03,
                                4e+03, 3e+03, 3e+03, 3e+03, 3e+03, 2e+03, 2e+03, 2e+03, 2e+03],
                        fwhms=[0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, ])
        Point2.set_SED(SED2)

        # load point to targets list
        targets = [Point1,Point2]

        # create the object
        mrs_scene_config = SceneConfig.makeScene(loglevel=0, targets=targets, background=background)

        return mrs_scene_config

    def run_sim(self, simcfg, scenecfg, simulatorcfg):
        """
        Take the scene, simulation, and simulator config objects and run
        simulation
        """
        mysim = MiriSimulation(simcfg, scenecfg, simulatorcfg)
        mysim.run()

    def run(self):

        simtor_cfg = self.make_simulator_config()
        scene_cfg = self.make_scene_config()

        # setup the simulation configs for each of the dispersers
        sim_cfgs = []
        for band in ('SHORT', 'MEDIUM', 'LONG'):
            sim_config = self.make_sim_config(band=band)
            sim_cfgs.append(sim_config)

        # run the fullband simualtions
        for sim_cfg in sim_cfgs:
            self.run_sim(sim_cfg, scene_cfg, simtor_cfg)

        # move to output_dir
        out_dirs = glob.glob('20*')
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        for d in out_dirs:
            shutil.move(d, os.path.join(self.output_dir, d))



my_sims = MRSFullband(exps=1, ints=1, groups=50, noise=True,
                      dither=True, output_dir='MRS_fullband')
my_sims.run()
