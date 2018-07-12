# -*- coding: utf-8 -*-
"""
Tools to visualise the ramp fit results in the optional
fitopt output of the ramp_fit step. Can be run from the command
line as:

> python plot_ramp_fits.py <fitopt file>

"""
import os
import glob
import shutil
import optparse
import time

import numpy as np
from astropy.io import fits
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

from jwst import datamodels


class RampPlots():
    """
    Make plots from the a ramp file or ramp data model.

    Parameters
    ----------
    ramp_filename :
        a ramp file from any pipeline step preceding ramp_fit

    xlim : list
        the limits of the xrange to plot

    ylim : list
        the limits of the yrange to plot

    """
    def __init__(self, ramp, xlim=[550, 650], ylim=[550, 650]):

        # check if ramp is a file or datamodel
        if type(ramp) is str:
            self.ramp_file = os.path.basename(ramp)
            self.ramp_dm = datamodels.open(ramp)

        else:
            self.ramp_file = ramp.meta.filename
            self.ramp_dm = ramp

        self.xlim = xlim
        self.ylim = ylim

    def plot_ramps(self):
        """
        plot the ramps
        """
        my_plot_ramps = self.ramp_dm.data[0, :, self.ylim[0]:self.ylim[1], self.xlim[0]:self.xlim[1]]
        (z, y, x) = my_plot_ramps.shape

        fig, axs = plt.subplots(1, 1, figsize=(12, 8))
        axs.set_xlabel('time (s)', fontsize=15)
        axs.set_ylabel('DN', fontsize=15)
        axs.set_title("ramps in {}".format(self.ramp_file), fontsize=12)

        for n in range(x):
            for m in range(y):
                axs.plot(my_plot_ramps[:, m, n], marker='.', markersize=0, linestyle='-', linewidth=0.5)
                axs.set_xlabel('groups')
                axs.set_ylabel('DN')

        plt.tight_layout()
        plot_name = 'ramps.pdf'
        try:
            os.remove(plot_name)
        except:
            pass
        fig.savefig(plot_name, dpi=100)

    def plot_groupdq(self):
        """
        plot the groupdq vector for each pixel and the location of pixels flagged with jump/do_not_use
        """
        my_plot_gdqs = self.ramp_dm.groupdq[0, :, self.ylim[0]:self.ylim[1], self.xlim[0]:self.xlim[1]]
        (z, y, x) = my_plot_gdqs.shape

        # plot the vectors
        fig, axs = plt.subplots(1, 1, figsize=(12, 8))
        axs.set_xlabel('time (s)', fontsize=15)
        axs.set_ylabel('GROUP_DQ', fontsize=15)
        axs.set_title("group_dq in {}".format(self.ramp_file), fontsize=12)

        jump_area = np.zeros((y, x))
        for n in range(x):
            for m in range(y):
                if np.any(my_plot_gdqs[:, m, n] > 1):
                    axs.plot(my_plot_gdqs[:, m, n], marker='.', markersize=0, linestyle='-', linewidth=0.5, c='r')
                    if self.ramp_dm.pixeldq[m, n] == 1:
                        jump_area[m, n] = 2
                    else:
                        jump_area[m, n] = 0
                else:
                    axs.plot(my_plot_gdqs[:, m, n], marker='.', markersize=0, linestyle='-', linewidth=0.5, c='b')
                    if self.ramp_dm.pixeldq[m, n] == 1:
                        jump_area[m, n] = 2
                    else:
                        jump_area[m, n] = 1

        plt.tight_layout()
        dq_vector_name = 'group_dq_vectors.pdf'
        try:
            os.remove(dq_vector_name)
        except:
            pass
        fig.savefig(dq_vector_name, dpi=100)

        # location of jumps/do not use
        fig, axs = plt.subplots(1, 1, figsize=(12, 8))

        axs.imshow(jump_area, cmap='seismic_r', interpolation='nearest', origin='lower', vmin=0, vmax=2)
        axs.set_title('group_dq including bad pixels', fontsize=15)

        plt.tight_layout()
        dq_image_name = 'group_dq_image.pdf'
        try:
            os.remove(dq_image_name)
        except:
            pass
        fig.savefig(dq_image_name, dpi=100)

    def plot_all(self):
        """
        run all plot methods
        """
        self.plot_ramps()
        self.plot_groupdq()


if __name__ == "__main__":
    # Parse arguments
    help_text = ""
    usage = "\n\n%prog <ramp file>\n"
    usage += "\nMakes a series of plots from a ramp file"

    parser = optparse.OptionParser(usage)
    (options,args) = parser.parse_args()

    # check for file
    if len(args) == 1:

        try:
            assert os.path.isfile(args[0])

        except AssertionError:
            print(help_text)
            time.sleep(1)  # Ensure help text appears before error messages.
            parser.error("ramp file does not exist...")
            sys.exit(1)

    else:
        print(help_text)
        time.sleep(1)  # Ensure help text appears before error messages.
        parser.error("must supply a rate file...")
        sys.exit(1)

    my_ramp_fits = RampPlots(args[0], xlim=[540, 560], ylim=[640, 660])
    my_ramp_fits.plot_all()
