# -*- coding: utf-8 -*-
"""
Make some test data with ONLY dark effect and see how the
pipeline handles the dark subrtaction and jump detection

Can be run from the command line as:

> python simulation_and_pipeline_dark_handling.py

"""
import os
import shutil
import time
import glob
import optparse
import numpy as np
import matplotlib.pyplot as plt
import logging

from mirisim import MiriSimulation
from mirisim.skysim import BBSed, Galaxy, Background
from mirisim.config_parser import SimConfig, SceneConfig, SimulatorConfig
from miri.datamodels import MiriRampModel

from jwst import datamodels
from jwst.pipeline import Detector1Pipeline
from jwst.stpipe import crds_client

# disable the logger output
logging.disable(logging.INFO)


class DarkTests:
    """
    Class to create MIRI ramp data with dark current only. Can optionally add noise.

    The simulated data is of a large galaxy so all pixels are well illuminated. We will
    run the simulation with MIRISim and run Detector1Pipeline. We will have a look at
    the group_dq after the jump steps to see how many jumps were flagged. If no noise
    or cosmic rays added this should be minimal.

    Only one exposure with one integration of 50 groups is simulated.

    Parameters
    ----------
    instrument : string
        for when simulations are required. Options are IMA or MRS

    input_dir : string
        the location of simulated and pipelined output when only plots are required. Must contain
        a MIRISim output file, and various pipeline outputs

    noise :  boolean
        when noise=True, readnoise, Poisson noise, and cosmic rays are included

    output_dir : string
        the location to save the simulated and pipelined output

    """
    def __init__(self, instrument=None, input_dir=None, noise=False, linearity=True, output_dir='dark_handling_tests'):

        self.instrument = instrument
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.linearity = linearity
        self.noise = noise

        self.ramp_file = None
        self.jump_file = None
        self.rate_file = None
        self.pre_dark_file = None
        self.post_dark_file = None

        # if the output directory exists, delete it and create a new one
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)

    def simulate_test_data(self):
        """
        Run the MIRISim simulations
        """
        # simulation config
        if self.instrument == 'IMA':
            sim_config = SimConfig.makeSim(name="IMA Simulation", rel_obsdate=0.0, scene="scene.ini", POP='IMA',
                                                    ConfigPath='IMA_FULL', Dither=False,StartInd=1, NDither=4,
                                                    DitherPat="ima_recommended_dither.dat", filter="F1130W",
                                                    readDetect= 'FULL', ima_mode= 'FAST', ima_exposures=1,
                                                    ima_integrations=1, ima_frames=50, disperser= 'SHORT', detector= 'SW',
                                                    mrs_mode= 'SLOW', mrs_exposures=5,mrs_integrations=4,mrs_frames=10)

            # scene config
            background = Background(level='low', gradient=5., pa=15.0, centreFOV=(0., 0.))

            SED1 = BBSed(Temp=300., wref=10., flux=1.e8)
            Gal1 = Galaxy(Cen=(0., 0.), n=1., re=200, q=0.99, pa=0.1)
            Gal1.set_SED(SED1)
            targets = [Gal1]

            scene_config = SceneConfig.makeScene(loglevel=0, background=background, targets=targets)

        elif self.instrument == 'MRS':
            sim_config = SimConfig.makeSim(name="MRS Simulation", rel_obsdate=0.0, scene="scene.ini",
                                                    POP='MRS', ConfigPath='MRS_1SHORT', Dither=False, StartInd=1,
                                                    NDither=4, DitherPat="mrs_recommended_dither.dat", filter="F1130W",
                                                    readDetect='FULL', ima_mode='FAST', ima_exposures=1,
                                                    ima_integrations=1, ima_frames=20, disperser='SHORT', detector='SW',
                                                    mrs_mode='FAST', mrs_exposures=1,  mrs_integrations=1, mrs_frames=50)

            # scene config
            background = Background(level='low',gradient=5.,pa=15.0,centreFOV=(0., 0.))

            SED1 = BBSed(Temp=300., wref=10., flux=5.e6)
            Gal1 = Galaxy(Cen=(0.,0.), n=1., re=2, q=0.99, pa=0.1)
            Gal1.set_SED(SED1)
            targets = [Gal1]

            scene_config = SceneConfig.makeScene(loglevel=0, background=background, targets=targets)

        # simulator config
        if self.noise:
            simulator_config = SimulatorConfig.makeSimulator(max_fsm=0.050, max_dither=20.0, mrs_ref_channel=1,
                                                                  mrs_ref_band="SHORT", tau_telescope=0.88,tau_eol=0.8,
                                                                  telescope_area=25.032, telescope_pupil_diam=6.6052,
                                                                  take_webbPsf=False, include_refpix=True,
                                                                  include_poisson=True, include_readnoise=True,
                                                                  include_badpix=True, include_dark=True,
                                                                  include_flat=False, include_gain=True,
                                                                  include_nonlinearity=self.linearity, include_drifts=True,
                                                                  include_latency=False, cosmic_ray_mode='NONE')
        else:
            simulator_config = SimulatorConfig.makeSimulator(max_fsm=0.050, max_dither=20.0, mrs_ref_channel=1,
                                                                  mrs_ref_band="SHORT", tau_telescope=0.88, tau_eol=0.8,
                                                                  telescope_area=25.032, telescope_pupil_diam=6.6052,
                                                                  take_webbPsf=False, include_refpix=True,
                                                                  include_poisson=False, include_readnoise=False,
                                                                  include_badpix=True, include_dark=True,
                                                                  include_flat=False, include_gain=True,
                                                                  include_nonlinearity=self.linearity, include_drifts=True,
                                                                  include_latency=False, cosmic_ray_mode='NONE')


        # run the simulation
        simulation = MiriSimulation(sim_config=sim_config, scene_config=scene_config,
                                    simulator_config=simulator_config, loglevel='DEBUG')
        simulation.run()

        # we only need the sim file so move it to output_dir and remove everthing else
        det_image_file = glob.glob(os.path.join(simulation.path_out, 'det_images', '*.fits'))[0]
        self.ramp_file = os.path.join(self.output_dir, os.path.basename(det_image_file))
        shutil.move(det_image_file, self.ramp_file)
        shutil.rmtree(simulation.path_out)

    def pipeline_test_data(self):
        """
        pipeline the test data, save the jump output
        """
        if self.linearity:
            Detector1Pipeline.call(self.ramp_file, save_results=True, output_dir=self.output_dir, output_use_model=True,
                                   steps={'ipc': {'skip': True},
                                          'rscd': {'skip': True},
                                          'lastframe': {'save_results': True,
                                                   'output_dir': self.output_dir},
                                          'dark_current': {'save_results': True,
                                                        'output_dir': self.output_dir},
                                          #'linearity': {'skip': True},
                                          'jump': {'save_results': True,
                                                   'output_dir': self.output_dir}})
        else:
            Detector1Pipeline.call(self.ramp_file, save_results=True, output_dir=self.output_dir, output_use_model=True,
                                   steps={'ipc': {'skip': True},
                                          'rscd': {'skip': True},
                                          'lastframe': {'save_results': True,
                                                        'output_dir': self.output_dir},
                                          'dark_current': {'save_results': True,
                                                           'output_dir': self.output_dir},
                                          'linearity': {'skip': True},
                                          'jump': {'save_results': True,
                                                   'output_dir': self.output_dir}})

        self.pre_dark_file = os.path.join(self.output_dir, 'step_lastframe.fits')
        self.post_dark_file = os.path.join(self.output_dir, 'step_dark_current.fits')
        self.jump_file = os.path.join(self.output_dir, 'step_jump.fits')
        self.rate_file = os.path.join(self.output_dir, 'step_rate.fits')

    def plot_jump_flags_image(self):
        """
        We want to plot the location of any jumps, i.e., anywhere collapsed group_dq > 2
        """
        dm = datamodels.MIRIRampModel(self.jump_file)
        fig, axs = plt.subplots(1, 1, figsize=(8, 8))

        group_dq = np.squeeze(dm.groupdq, axis=0)
        group_dq = np.sum(group_dq, axis=0)

        axs.imshow(group_dq, cmap='gray', interpolation='nearest', origin='lower', vmin=2, vmax=4)

        plt.tight_layout()

        plot_name = os.path.join(self.output_dir, 'jump_flags.pdf')
        try:
            os.remove(plot_name)
        except:
            pass

        fig.savefig(plot_name, dpi=100)

    def _get_dark_cdp(self):
        """
        Convenience function to determine DARK CDP used in the simulations from the HISTORY in MIRISim file metadata
        """
        dm = MiriRampModel(self.ramp_file)

        for row in dm.history['entries']:
            my_entry = (row['description'])
            if 'DARK' in my_entry:
                strings = my_entry.split()
                for s in strings:
                    if 'DARK' in s:
                        cdp_file = s
                        cdp_file = cdp_file.replace("'", "")
                        #print("Found: {}".format(cdp_file))

        cdp = datamodels.DarkMIRIModel(os.path.join(os.environ['CDP_DIR'], cdp_file))

        return cdp

    def plot_ramps_pre_post_correction(self, pixel=[500, 500]):
        """
        We want to plot a test ramp pre and post dark correction and the dark itself.
        Also print the properties of the dark file to check that they match.

        Parameters
        ----------
        pixel : list
            the pixel to plot the ramp for
        """
        in_dm = datamodels.MIRIRampModel(self.pre_dark_file)
        out_dm = datamodels.MIRIRampModel(self.post_dark_file)

        nints, ngroups, nrows, ncols = in_dm.data.shape
        group = range(1, ngroups + 1)

        fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        axs[0].plot(group, in_dm.data[0, :, pixel[1], pixel[0]], c='r', marker='^', markersize=0, linestyle='-',
                    linewidth=1, label='pre-dark')
        axs[0].plot(group, out_dm.data[0, :, pixel[1], pixel[0]], c='b', marker='o', markersize=0, linestyle='-',
                    linewidth=1, label='post-dark')
        axs[0].set_title('dark current correction', fontsize=15)
        axs[0].set_ylabel('DN', fontsize=15)
        axs[0].set_xlim(-1, max(group) + 1)
        axs[0].legend(prop={'size': 12}, loc=0)

        # get CDP
        cdp = self._get_dark_cdp()

        # get CRDS
        crds_filename = crds_client.get_reference_file(in_dm, 'dark')
        crds = datamodels.DarkMIRIModel(crds_filename)

        axs[1].plot(group, in_dm.data[0, :, pixel[1], pixel[0]] - out_dm.data[0, :, pixel[1], pixel[0]], c='g',
                    marker='s', markersize=0, linestyle='-', linewidth=1, label='difference')
        axs[1].plot(group, cdp.data[0, 0:ngroups, pixel[1], pixel[0]], c='orange',
                    marker='o', markersize=5, linestyle='-', linewidth=0, label='CDP')
        axs[1].plot(group, crds.data[0, 0:ngroups, pixel[1], pixel[0]], c='red',
                    marker='x', markersize=9, linestyle='-', linewidth=0, label='CRDS')
        axs[1].set_ylabel('DN', fontsize=15)
        axs[1].set_xlabel('group', fontsize=15)
        axs[1].legend(prop={'size': 12}, loc=0)

        plt.tight_layout(h_pad=0)

        plot_name = os.path.join(self.output_dir, 'dark_correction.pdf')
        try:
            os.remove(plot_name)
        except:
            pass

        fig.savefig(plot_name, dpi=100)

    def plot_groupdq_flags(self, pixel=[500, 500]):
        """
        We want to plot the ramp and the location of the jumps via the groupdq.

        Parameters
        ----------
        pixel : list
            the pixel to plot the ramp for
        """
        dm = datamodels.MIRIRampModel(self.jump_file)
        nints, ngroups, nrows, ncols = dm.data.shape
        group = range(1, ngroups + 1)

        # plot--------------------------------------
        fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # first integration for input/output ramps
        axs[0].plot(group, dm.data[0, :, pixel[1], pixel[0]], c='b', marker='o', linestyle='-', linewidth=1,
                    markersize=0, label='ramp')
        axs[0].set_title('jumps flagged', fontsize=15)
        axs[0].set_ylabel('DN', fontsize=15)
        axs[0].set_xlim(-1, max(group) + 1)
        axs[0].set_ylim(min(dm.data[0, :, pixel[1], pixel[0]]) - 200, max(dm.data[0, :, pixel[1], pixel[0]]) + 200)
        axs[0].legend(prop={'size': 12}, loc=2)

        # input and output flag values (setting a slight offset for the output flags)
        axs[1].plot(group, dm.groupdq[0, :, pixel[1], pixel[0]], c='g', marker='o', linestyle='-', linewidth=1,
                    markersize=0, label='groupdq')
        axs[1].plot([-10, 100], [4, 4], linestyle='--', linewidth=2, c='r', alpha=0.7, label='jump threshold')
        axs[1].set_ylabel('flag values', fontsize=15)
        axs[1].set_ylabel('group dq', fontsize=15)
        axs[1].set_xlim(-1, max(group) + 1)
        axs[1].set_ylim(-0.5, 6.5)
        axs[1].legend(prop={'size': 12}, loc=2)

        # draw lines to show the groups which have been flagged as jumps
        for n, val in enumerate(group):
            if (dm.groupdq[0, n, pixel[1], pixel[0]] >= 4):
                axs[0].plot([n + 1, n + 1], [min(dm.data[0, :, pixel[1], pixel[0]]) - 200,
                                             max(dm.data[0, :, pixel[1], pixel[0]]) + 200], linestyle='--',
                            linewidth=0.3, c='k')
                axs[1].plot([n + 1, n + 1], [-1, 6], linestyle='--', linewidth=0.3, c='k')

        plt.tight_layout()

        plot_name = os.path.join(self.output_dir, 'jump_flags_ramp.pdf')
        try:
            os.remove(plot_name)
        except:
            pass

        fig.savefig(plot_name, dpi=100)

    def run(self):
        """
        run simulations, pipeline and all test plots
        """
        self.simulate_test_data()
        self.pipeline_test_data()
        self.plot_jump_flags_image()
        self.plot_groupdq_flags(pixel=[884, 550])
        self.plot_ramps_pre_post_correction(pixel=[884, 550])

    def run_plots(self):
        """
        run test plots only
        """
        # load the files
        self.pre_dark_file = os.path.join(self.input_dir, 'step_lastframe.fits')
        self.post_dark_file = os.path.join(self.input_dir, 'step_dark_current.fits')
        self.jump_file = os.path.join(self.input_dir, 'step_jump.fits')
        self.rate_file = os.path.join(self.input_dir, 'step_rate.fits')
        self.ramp_file = glob.glob(os.path.join(self.input_dir, '*.fits'))[0]

        # plots
        self.plot_jump_flags_image()
        self.plot_groupdq_flags(pixel=[884, 550])
        self.plot_ramps_pre_post_correction(pixel=[884, 550])


if __name__ == "__main__":
    # Parse arguments
    help_text = ""
    usage = "\n\n%prog <instrument>\n"
    usage += "\nTest the dark addition (MIRISim) and correction (cal pipeline). Instrument must be IMA or MRS"

    parser = optparse.OptionParser(usage)
    (options,args) = parser.parse_args()

    assert len(args) > 0, "Must supply either instrument (IMA or MRS) or directory containing MIRISim/pipeline output"

    # check what the input is and run appropriate class method
    if (args[0] == 'IMA') or (args[0] == 'MRS'):
        output_dir = args[0] + '_dark_handling_tests'
        my_dark_tests = DarkTests(instrument=args[0], noise=False, linearity=False, output_dir=output_dir)
        my_dark_tests.run()

    else:
        assert os.path.exists(args[0]), "Can't find directory {}".format(args[0])

        output_dir = args[0] + '_dark_handling_test_plots'
        my_dark_tests = DarkTests(input_dir=args[0], noise=False, output_dir=output_dir)
        my_dark_tests.run_plots()