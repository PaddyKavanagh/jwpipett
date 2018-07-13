# -*- coding: utf-8 -*-
"""
test handling of jumps in ramp fitting

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
from jwst.ramp_fitting import ramp_fit_step

from plot_ramp_fits import RampFits


def test_jump_handling(jump_file, output_dir='jump_handling_tests'):
    """
    make plots of the ramp fits to assess jump handling
    """
    # if the output directory exists, delete it and create a new one
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    # load the jump file to datamodel
    jdm = datamodels.MIRIRampModel(jump_file)

    # run ramp fitting
    ramp_fit_step.RampFitStep.call(jdm, save_results=True, save_opt=True, output_file='rate.fits',
                                   opt_name='fitopt.fits', output_dir=output_dir)

    # get the output filenames
    fitopt_file = os.path.join(output_dir, 'fitopt_fitopt.fits')
    rate_file = os.path.join(output_dir, 'rate_0_rampfitstep.fits')

    # plot the fits using the RampFits class
    my_plots = RampFits(fitopt_file, jump_file, rate_file, xlim=[500, 530], ylim=[500, 530])
    my_plots.plot_all()

if __name__ == "__main__":
    # Parse arguments
    help_text = ""
    usage = "\n\n%prog <jump_file>\n"
    usage += "\nCheck the handling of jumps in ramp fitting"

    parser = optparse.OptionParser(usage)
    (options,args) = parser.parse_args()

    # check for file
    if len(args) == 1:
        try:
            assert os.path.isfile(args[0])

        except AssertionError:
            print(help_text)
            time.sleep(1)  # Ensure help text appears before error messages.
            parser.error("jump file does not exist...")
            sys.exit(1)

    test_jump_handling(jump_file=args[0])

