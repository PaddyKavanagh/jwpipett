# -*- coding: utf-8 -*-
"""
make_demo_data.py

produce a set of simulations for use with pipeline notebooks

pkavanagh
"""
import os, shutil, time, glob
import optparse
import numpy as np
from mirisim import MiriSimulation
from mirisim.skysim import *
from mirisim.config_parser import *
from scene_from_cat import make_simple_cat_scene_obj
from jwst import datamodels


class IMADemoData():
    """
    Class to create IMA demo data
    """
    def __init__(self):
        """
        setup the MIRISim config objects
        """
        self.ima_sim_config = SimConfig.makeSim(name="IMA Simulation",rel_obsdate=0.0,scene="scene.ini",
                                       POP='IMA',ConfigPath='IMA_FULL',Dither=True,StartInd=1,
                                       NDither=4,DitherPat="ima_recommended_dither.dat",filter="F1130W",
                                       readDetect= 'FULL',ima_mode= 'FAST',ima_exposures=1,ima_integrations=2,
                                       ima_frames=50,disperser= 'SHORT',detector= 'SW',mrs_mode= 'SLOW',
                                       mrs_exposures=5,mrs_integrations=4,mrs_frames=10)

        self.simulator_config = SimulatorConfig.makeSimulator(max_fsm=0.050,max_dither=20.0,mrs_ref_channel=1,
                                                         mrs_ref_band="SHORT",tau_telescope=0.88,tau_eol=0.8,
                                                         telescope_area=25.032,telescope_pupil_diam=6.6052,
                                                         take_webbPsf=False,include_refpix=True,include_poisson=True,
                                                         include_readnoise=True,include_badpix=True,include_dark=True,
                                                         include_flat=True,include_gain=True,include_nonlinearity=True,
                                                         include_drifts=True,include_latency=False,
                                                              cosmic_ray_mode='SOLAR_MIN')

        self.ima_scene_config = make_simple_cat_scene_obj(target_coords=[1.0,1.0],random=True, source_num=500,
                                                    centre_coords=[1.0,1.0])


    def run(self):
        """
        run the simulation and rename output dir
        """
        self.ima_sim = MiriSimulation(sim_config=self.ima_sim_config,scene_config=self.ima_scene_config,
                                simulator_config=self.simulator_config,loglevel= 'DEBUG')

        self.ima_sim.run()

        new_dir = "IMA"
        if os.path.isdir(new_dir):
            shutil.rmtree(new_dir)

        os.rename(self.ima_sim.path_out , new_dir)



class LRSSlitDemoData():
    """
    Class to create LRS Slit demo data
    """
    def __init__(self):
        """
        setup the MIRISim config objects
        """
        self.lrs_slit_sim_config = SimConfig.makeSim(name="LRS Simulation",rel_obsdate=0.0,scene="scene.ini",
                                       POP='IMA',ConfigPath='LRS_SLIT',Dither=True, StartInd=1,
                                       NDither=2,DitherPat="lrs_recommended_dither.dat",filter="P750L",
                                       readDetect= 'FULL',ima_mode= 'FAST',ima_exposures=1,ima_integrations=2,
                                       ima_frames=50,disperser= 'SHORT',detector= 'SW',mrs_mode= 'SLOW',
                                       mrs_exposures=5,mrs_integrations=4,mrs_frames=10)

        self.lrs_slit_bkg_sim_config = SimConfig.makeSim(name="LRS Simulation", rel_obsdate=0.0, scene="scene.ini",
                                                     POP='IMA', ConfigPath='LRS_SLIT', Dither=False, StartInd=2,
                                                     NDither=2, DitherPat="lrs_recommended_dither.dat", filter="P750L",
                                                     readDetect='FULL', ima_mode='FAST', ima_exposures=1,
                                                     ima_integrations=2,
                                                     ima_frames=50, disperser='SHORT', detector='SW', mrs_mode='SLOW',
                                                     mrs_exposures=5, mrs_integrations=4, mrs_frames=10)

        self.simulator_config = SimulatorConfig.makeSimulator(max_fsm=0.050,max_dither=20.0,mrs_ref_channel=1,
                                                         mrs_ref_band="SHORT",tau_telescope=0.88,tau_eol=0.8,
                                                         telescope_area=25.032,telescope_pupil_diam=6.6052,
                                                         take_webbPsf=False,include_refpix=True,include_poisson=True,
                                                         include_readnoise=True,include_badpix=True,include_dark=True,
                                                         include_flat=True,include_gain=True,include_nonlinearity=True,
                                                         include_drifts=True,include_latency=False,
                                                              cosmic_ray_mode='SOLAR_MIN')

        # set background
        background = Background(level='low',gradient=5.,pa=15.0,centreFOV=(0., 0.))

        SED1 = PLSed(alpha=0.0, wref=10., flux=1.e4)
        Point1 = Point(Cen=(0., 0.))
        Point1.set_SED(SED1)

        Point2 = Point(Cen=(0., 0.))

        SED2 = LinesSed(wavels=[6, 7, 8, 9, 10, 11, 12],
                        fluxes=[5e+03, 4e+03, 3e+03, 3e+03, 3e+03, 3e+03, 3e+03],
                        fwhms=[0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
        Point2.set_SED(SED2)

        # set a background SED and object
        SED3 = BBSed(Temp=300., wref=10., flux=1.e-1)
        Gal1 = Galaxy(Cen=(0.,0.), n=2., re=50, q=0.9, pa=0.1)
        Gal1.set_SED(SED3)

        # load point to targets list
        targets = [Point1, Point2, Gal1]
        bkg_targets = [Gal1]


        # create the object
        self.lrs_scene_config = SceneConfig.makeScene(loglevel=0, background=background, targets=targets)

        self.lrs_bkg_scene_config = SceneConfig.makeScene(loglevel=0, background=background, targets=bkg_targets)


    def run(self):
        """
        run the simulation and rename output dir
        """
        self.lrs_sim = MiriSimulation(sim_config=self.lrs_slit_sim_config,scene_config=self.lrs_scene_config,
                                simulator_config=self.simulator_config,loglevel= 'DEBUG')

        self.lrs_sim.run()

        new_dir = "LRS-SLIT"
        if os.path.isdir(new_dir):
            shutil.rmtree(new_dir)

        os.rename(self.lrs_sim.path_out , new_dir)

        # background
        self.lrs_bkg_sim = MiriSimulation(sim_config=self.lrs_slit_bkg_sim_config,scene_config=self.lrs_bkg_scene_config,
                                simulator_config=self.simulator_config,loglevel= 'DEBUG')

        self.lrs_bkg_sim.run()

        new_dir = "LRS-SLIT_BKG"
        if os.path.isdir(new_dir):
            shutil.rmtree(new_dir)

        os.rename(self.lrs_bkg_sim.path_out, new_dir)


class LRSSlitlessDemoData():
    """
    Class to create LRS Slitless demo data. Need to simulate n integrations where there is variation
    in an absorption line and combine these into a combined 4d array
    """
    def __init__(self):
        """
        setup the MIRISim config objects
        """
        self.lrs_slitless_sim_config = SimConfig.makeSim(name="LRS Simulation",rel_obsdate=0.0,scene="scene.ini",
                                       POP='IMA',ConfigPath='LRS_SLITLESS',Dither=False,StartInd=0,NDither=2,
                                       DitherPat="lrs_recommended_dither.dat",filter="P750L",readDetect= 'FULL',
                                       ima_mode= 'FAST',ima_exposures=3,ima_integrations=20,ima_frames=20,
                                       disperser= 'SHORT',detector= 'SW',mrs_mode= 'SLOW',mrs_exposures=5,
                                       mrs_integrations=4,mrs_frames=10)

        self.simulator_config = SimulatorConfig.makeSimulator(max_fsm=0.050,max_dither=20.0,mrs_ref_channel=1,
                                                         mrs_ref_band="SHORT",tau_telescope=0.88,tau_eol=0.8,
                                                         telescope_area=25.032,telescope_pupil_diam=6.6052,
                                                         take_webbPsf=False,include_refpix=True,include_poisson=True,
                                                         include_readnoise=True,include_badpix=True,include_dark=True,
                                                         include_flat=True,include_gain=True,include_nonlinearity=True,
                                                         include_drifts=True,include_latency=False,
                                                              cosmic_ray_mode='SOLAR_MIN')


    def make_int_scene(self, factor):
        """
        return a scene config where the strength of an absorption line varies according to factor.
        """

        # set background
        background = Background(level= 'low',gradient=5.,pa=15.0,centreFOV=(0., 0.))

        # set SED
        Point1 = Point(Cen=(0., 0.))
        Point2 = Point(Cen=(0., 0.))
        Point3 = Point(Cen=(0., 0.))

        bb = BBSed(Temp=300., wref=10., flux=1.e4)

        lines = LinesSed(wavels=[6, 7, 8, 9, 10, 11, 12],
                        fluxes=[5e+03, 4e+03, 3e+03, 3e+03, 3e+03, 3e+03, 3e+03],
                        fwhms=[0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])

        abs_line = LinesSed(wavels=[9.0], fluxes=[-1e4 * factor], fwhms=[1.0])

        Point1.set_SED(bb)
        Point2.set_SED(lines)
        Point3.set_SED(abs_line)

        # load point to targets list
        targets = [Point1, Point2, Point3]

        # return the scene object
        return SceneConfig.makeScene(loglevel=0, background=background, targets=targets)

    def sim_ints(self, factors):
        """
        simulate the ints
        """
        for n, f in enumerate(factors):

            int_scene = self.make_int_scene(f)

            int_sim = SimConfig.makeSim(name="LRS Simulation", rel_obsdate=0.0, scene="scene.ini",
                              POP='IMA', ConfigPath='LRS_SLITLESS', Dither=False, StartInd=0, NDither=2,
                              DitherPat="lrs_recommended_dither.dat", filter="P750L", readDetect='FULL',
                              ima_mode='FAST', ima_exposures=1, ima_integrations=1, ima_frames=20,
                              disperser='SHORT', detector='SW', mrs_mode='SLOW', mrs_exposures=5,
                              mrs_integrations=4, mrs_frames=10)

            lrs_sim = MiriSimulation(sim_config=int_sim, scene_config=int_scene,
                                          simulator_config=self.simulator_config, loglevel='DEBUG')

            lrs_sim.run()

            new_dir = "int%d" % n
            if os.path.isdir(new_dir):
                shutil.rmtree(new_dir)

            os.rename(lrs_sim.path_out, new_dir)

    def replace_ints(self, dir, int_dirs):
        """
        replace the integrations in the full array with the int arrays
        """
        my_output_files = glob.glob(os.path.join(dir, 'det_images', '*.fits'))

        for f in my_output_files:
            out_dm = datamodels.open(f)

            for n, i in enumerate(int_dirs):
                my_int_file = glob.glob(os.path.join(i, 'det_images', '*.fits'))[0]
                int_dm = datamodels.open(my_int_file)

                out_dm.data[n] = int_dm.data

            out_dm.save(f)


    def run(self):
        """
        run the simulation and rename output dir
        """
        # simulate the ints
        factors1 = np.arange(10)
        factors2 = np.flipud(factors1)
        factors = np.concatenate((factors1, factors2))


        self.sim_ints(factors)

        # simulate with len(factors) ints
        self.lrs_scene_config = SceneConfig.from_default()
        self.lrs_sim = MiriSimulation(sim_config=self.lrs_slitless_sim_config,scene_config=self.lrs_scene_config,
                                simulator_config=self.simulator_config,loglevel= 'DEBUG')

        self.lrs_sim.run()

        new_dir = "LRS-SLITLESS"
        if os.path.isdir(new_dir):
            shutil.rmtree(new_dir)

        os.rename(self.lrs_sim.path_out , new_dir)

        # get the names of the int directors
        int_dirs = sorted(glob.glob('int*'))

        # replace the
        self.replace_ints(new_dir, int_dirs)

        # clean up
        # self.add_var_abs(det_images)


class MRSDemoData():
    """
    Class to create MRS demo data
    """
    def __init__(self):
        """
        setup the MIRISim config objects
        """
        self.mrs_sim_config = SimConfig.makeSim(name="MRS Simulation",rel_obsdate=0.0,scene="scene.ini",
                                   POP='MRS',ConfigPath='MRS_1SHORT',Dither=True,StartInd=1,NDither=4,
                                   DitherPat="mrs_recommended_dither.dat",filter="F1130W",readDetect= 'FULL',
                                   ima_mode= 'FAST',ima_exposures=1,ima_integrations=1,ima_frames=20,
                                   disperser= 'SHORT',detector= 'SW',mrs_mode= 'FAST',mrs_exposures=1,
                                   mrs_integrations=2,mrs_frames=50)

        self.simulator_config = SimulatorConfig.makeSimulator(max_fsm=0.050,max_dither=20.0,mrs_ref_channel=1,
                                                         mrs_ref_band="SHORT",tau_telescope=0.88,tau_eol=0.8,
                                                         telescope_area=25.032,telescope_pupil_diam=6.6052,
                                                         take_webbPsf=False,include_refpix=True,include_poisson=False,
                                                         include_readnoise=False,include_badpix=True,include_dark=False,
                                                         include_flat=True,include_gain=True,include_nonlinearity=True,
                                                         include_drifts=True,include_latency=False,
                                                              cosmic_ray_mode='SOLAR_MIN')

        # set background
        background = Background(level= 'low',gradient=5.,pa=15.0,centreFOV=(0., 0.))

        # set PL SED
        SED1 = PLSed(alpha=0.0, wref=10., flux=1.e5)
        Point1 = Point(Cen=(0., 0.))
        Point1.set_SED(SED1)

        Point2 = Point(Cen=(0., 0.))
        SED2 = LinesSed(wavels=[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
                        fluxes=[5e+04, 4e+04, 3e+04, 2e+04, 2e+04, 1e+04, 1e+04, 9e+04, 8e+03, 7e+03, 6e+03, 5e+03,
                                4e+03, 3e+03, 3e+03, 3e+03, 3e+03, 2e+03, 2e+03, 2e+03, 2e+03],
                        fwhms=[0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, ])
        Point2.set_SED(SED2)

        # load point to targets list
        targets = [Point1,Point2]

        # create the object
        self.mrs_scene_config = SceneConfig.makeScene(loglevel=0,background=background,targets=targets)


    def run(self):
        """
        run the simulation and rename output dir
        """
        self.mrs_sim = MiriSimulation(sim_config=self.mrs_sim_config,scene_config=self.mrs_scene_config,
                                simulator_config=self.simulator_config,loglevel= 'DEBUG')

        self.mrs_sim.run()

        new_dir = "MRS"
        if os.path.isdir(new_dir):
            shutil.rmtree(new_dir)

        os.rename(self.mrs_sim.path_out , new_dir)


if __name__ == "__main__":
    # Parse arguments
    help_text = ""
    usage = "\n\n%prog <instruments>\n"
    usage += "\nGenerate demo data for given instruments (IMA,MRS,LRS-SLIT,LRS-SLITLESS). "
    usage += "If no instrument provided, simulates all."

    parser = optparse.OptionParser(usage)
    (options,args) = parser.parse_args()

    instruments = args

    # check for correct instruments
    if len(args) >= 1:
        try:
            for instrument in instruments: assert instrument in ['IMA','MRS','LRS-SLIT', 'LRS-SLITLESS']

        except AssertionError:
            print(help_text)
            time.sleep(1) # Ensure help text appears before error messages.
            parser.error("Instrument not recognised, must be 'IMA', 'MRS', 'LRS-SLIT' or 'LRS-SLITLESS")
            sys.exit(1)

    else:
        instruments = ['IMA', 'MRS', 'LRS-SLIT', 'LRS-SLITLESS']

    for instrument in instruments:
        if instrument == 'IMA':
            print("Starting IMA demo simulations...")
            ima_demo = IMADemoData()
            ima_demo.run()
        elif instrument == 'LRS-SLIT':
            print("Starting LRS-SLIT demo simulations...")
            lrs_slit_demo = LRSSlitDemoData()
            lrs_slit_demo.run()
        elif instrument == 'LRS-SLITLESS':
            print("Starting LRS-SLITLESS demo simulations...")
            lrs_slitless_demo = LRSSlitlessDemoData()
            lrs_slitless_demo.run()
        elif instrument == 'MRS':
            print("Starting MRS demo simulations...")
            mrs_demo = MRSDemoData()
            mrs_demo.run()
