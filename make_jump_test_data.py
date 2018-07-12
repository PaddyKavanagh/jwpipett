# -*- coding: utf-8 -*-
"""
Make test data to investigate the behaviour of pipeline
dealing with jumps.

Can be run from the command line as:

> python make_jump_test_data.py

"""
import os
import shutil
import time
import glob
import optparse
import numpy as np
import matplotlib.pyplot as plt

from mirisim import MiriSimulation
from mirisim.skysim import BBSed, Galaxy, Background
from mirisim.config_parser import SimConfig, SceneConfig, SimulatorConfig

from jwst import datamodels
from jwst.pipeline import Detector1Pipeline


class JumpTestData:
    """
    Class to create IMA data with jumps set for specfic groups to inspect pipeline
    behaviour in dealing with ramp_fitting of data with jumps.

    The simulated data is of a large galaxy so all pixels are well illuminated. We will
    run the simulation with MIRISim, run Detector1Pipeline up to before the
    jump step, and then manually edit the groupdq array to set jumps manually.

    Only one exposure with one integration of 50 groups is simulated.

    Parameters
    ----------
    noise :  boolean
        when noise=True, readnoise, Poisson noise, and cosmic rays are included

    output_dir :
        the location to save the simulated and pipelined output

    """
    def __init__(self, noise=True, output_dir='jump_test_data', mode=1):

        self.output_dir = output_dir
        self.noise = noise
        self.mode = mode

        self.sim_file = None
        self.jump_file = None

        # if the output directory exists, delete it and create a new one
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)

    def simulate_test_data(self):
        """
        Run the MIRISim simulations
        """

        # simulation config
        ima_sim_config = SimConfig.makeSim(name="IMA Simulation", rel_obsdate=0.0, scene="scene.ini", POP='IMA',
                                                ConfigPath='IMA_FULL', Dither=False,StartInd=1, NDither=4,
                                                DitherPat="ima_recommended_dither.dat", filter="F1130W",
                                                readDetect= 'FULL', ima_mode= 'FAST', ima_exposures=1,
                                                ima_integrations=1, ima_frames=50, disperser= 'SHORT', detector= 'SW',
                                                mrs_mode= 'SLOW', mrs_exposures=5,mrs_integrations=4,mrs_frames=10)

        # simulator config
        if self.noise:
            simulator_config = SimulatorConfig.makeSimulator(max_fsm=0.050, max_dither=20.0, mrs_ref_channel=1,
                                                                  mrs_ref_band="SHORT", tau_telescope=0.88,tau_eol=0.8,
                                                                  telescope_area=25.032, telescope_pupil_diam=6.6052,
                                                                  take_webbPsf=False, include_refpix=True,
                                                                  include_poisson=True, include_readnoise=True,
                                                                  include_badpix=True, include_dark=True,
                                                                  include_flat=True, include_gain=True,
                                                                  include_nonlinearity=True, include_drifts=True,
                                                                  include_latency=False, cosmic_ray_mode='SOLAR_MIN')
        else:
            simulator_config = SimulatorConfig.makeSimulator(max_fsm=0.050, max_dither=20.0, mrs_ref_channel=1,
                                                                  mrs_ref_band="SHORT", tau_telescope=0.88, tau_eol=0.8,
                                                                  telescope_area=25.032, telescope_pupil_diam=6.6052,
                                                                  take_webbPsf=False, include_refpix=True,
                                                                  include_poisson=False, include_readnoise=False,
                                                                  include_badpix=True, include_dark=True,
                                                                  include_flat=True, include_gain=True,
                                                                  include_nonlinearity=True, include_drifts=True,
                                                                  include_latency=False, cosmic_ray_mode='NONE')

        # scene config
        background = Background(level='low',gradient=5.,pa=15.0,centreFOV=(0., 0.))

        SED1 = BBSed(Temp=300., wref=10., flux=1.e8)
        Gal1 = Galaxy(Cen=(0.,0.), n=1., re=200, q=0.99, pa=0.1)
        Gal1.set_SED(SED1)
        targets = [Gal1]

        ima_scene_config = SceneConfig.makeScene(loglevel=0, background=background, targets=targets)

        # run the simulation
        simulation = MiriSimulation(sim_config=ima_sim_config, scene_config=ima_scene_config,
                                    simulator_config=simulator_config, loglevel='DEBUG')
        simulation.run()

        # we only need the sim file so move it to output_dir and remove everthing else
        det_image_file = glob.glob(os.path.join(simulation.path_out, 'det_images', '*.fits'))[0]
        self.sim_file = os.path.join(self.output_dir, os.path.basename(det_image_file))
        shutil.move(det_image_file, self.sim_file)
        shutil.rmtree(simulation.path_out)

    def pipeline_test_data(self):
        """
        pipeline the test data to before the jump step, i.e., up to refpix
        """
        Detector1Pipeline.call(self.sim_file,
                               steps={'ipc': {'skip': True},
                                      'rscd': {'skip': True},
                                      'refpix': {'save_results': True,
                                                 'output_dir': self.output_dir},
                                      'jump': {'skip': True},
                                      'ramp_fit': {'skip': True}})

        self.jump_file = os.path.join(self.output_dir, 'step_refpix.fits')
        better_name = os.path.join(self.output_dir, 'jump_file.fits')
        shutil.move(self.jump_file, better_name)
        self.jump_file = better_name

    def set_jump_flags(self):
        """
        manually set the jump flags. We want to set:

        MODE 1=======================================
        top left quadrant:      no jumps
        top right quadrant:     jump at group 2  (after first frame)
        bottom left quadrant:   jump at group 24 (midpoint)
        bottom right quadrant:   jump at group 48 (before last frame)

        MODE 2=======================================
        top left quadrant:      jump at group 3  (after second frame - segment with only 1 group)
        top right quadrant:     jump at group 4  (after third frame - segment with 2 groups)
        bottom left quadrant:   jump at group 46 (before second last frame - segment with 2 groups)
        bottom right quadrant:   jump at group 47 (before second last frame - segment with only 1 group)

        MODE 3=======================================
        all:      one half has jump at group 2
        """

        with datamodels.MIRIRampModel(self.jump_file) as dm:

            nints, ngroups, nrows, ncols = dm.groupdq.shape
            mid_cols = int(ncols / 2)
            mid_rows = int(nrows / 2)

            if self.mode == 1:
                for n in range(ngroups):
                    if n == 1:
                        dm.groupdq[0, n, 0:mid_rows, mid_cols:] = 4
                    if n == 24:
                        dm.groupdq[0, n, mid_rows:, 0:mid_cols] = 4
                    if n == 48:
                        dm.groupdq[0, n, mid_rows:, mid_cols:] = 4

            elif self.mode == 2:
                for n in range(ngroups):
                    if n == 2:
                        dm.groupdq[0, n, 0:mid_rows, 0:mid_cols] = 4
                    if n == 3:
                        dm.groupdq[0, n, 0:mid_rows, mid_cols:] = 4
                    if n == 46:
                        dm.groupdq[0, n, mid_rows:, 0:mid_cols] = 4
                    if n == 47:
                        dm.groupdq[0, n, mid_rows:, mid_cols:] = 4

            elif self.mode == 2:
                for n in range(ngroups):
                    if n == 2:
                        dm.groupdq[0, n, mid_rows:, :] = 4

            dm.save(self.jump_file)

    def plot_jump_flags(self):
        """
        quick plot of the collapse groupdq to check this script
        """
        dm = datamodels.MIRIRampModel(self.jump_file)
        fig, axs = plt.subplots(2, 2, figsize=(12, 8))
        axs = axs.ravel()

        axs[0].imshow(dm.groupdq[0, 2], cmap='seismic_r', interpolation='nearest', origin='lower', vmin=0, vmax=4)
        axs[1].imshow(dm.groupdq[0, 1], cmap='seismic_r', interpolation='nearest', origin='lower', vmin=0, vmax=4)
        axs[2].imshow(dm.groupdq[0, 24], cmap='seismic_r', interpolation='nearest', origin='lower', vmin=0, vmax=4)
        axs[3].imshow(dm.groupdq[0, 49], cmap='seismic_r', interpolation='nearest', origin='lower', vmin=0, vmax=4)

        plt.tight_layout()
        plt.show()

    def run(self):
        """
        create the test data
        """
        self.simulate_test_data()
        self.pipeline_test_data()
        self.set_jump_flags()
        #self.plot_jump_flags()


if __name__ == "__main__":
    # Parse arguments
    help_text = ""
    usage = "\n\n%prog \n"
    usage += "\nMake test data simulating jump output for ramp_fit testing"

    parser = optparse.OptionParser(usage)
    (options,args) = parser.parse_args()

    my_test_data = JumpTestData(noise=False, mode=3)
    my_test_data.run()
