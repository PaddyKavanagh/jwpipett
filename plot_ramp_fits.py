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
from plot_ramps import RampPlots


class RampFits:
    """
    Make plots from the optional ramp_fit output file.

    Parameters
    ----------
    fitopt_filename :
        the fitopt file produced from save_opt=True in ramp_fit

    rate_filename :
        the output file of the ramp fitting

    ramp_filename :
        the output file of the jump step

    xlim : list
        the limits of the xrange to plot

    ylim : list
        the limits of the yrange to plot

    """
    def __init__(self, fitopt_filename, ramp_filename, rate_filename, xlim=[550, 650], ylim=[550, 650]):

        self.fitopt_file = os.path.basename(fitopt_filename)
        self.ramp_file = os.path.basename(ramp_filename)
        self.rate_file = os.path.basename(rate_filename)

        self.fitopt_dm = datamodels.RampFitOutputModel(fitopt_filename)
        self.ramp_dm = datamodels.open(ramp_filename)
        self.rate_dm = datamodels.open(rate_filename)

        self.xlim = xlim
        self.ylim = ylim

    def plot_fits(self):
        """
        plot the best fits in the fitopt file with the ramps
        """
        group = np.arange(1, self.ramp_dm.data[0, :, 512, 512].shape[0] + 1, 1)

        my_plot_area = self.ramp_dm.data[0, :, self.ylim[0]:self.ylim[1], self.xlim[0]:self.xlim[1]]
        (z, y, x) = my_plot_area.shape

        fig, axs = plt.subplots(1, 1, figsize=(12, 8))
        axs.set_xlabel('time (s)', fontsize=15)
        axs.set_ylabel('DN', fontsize=15)
        axs.set_title("best fits in {}, plotted on {} ramps".format(self.fitopt_file, self.ramp_file), fontsize=12)

        my_slopes = self.fitopt_dm.slope[0, 0, self.ylim[0]:self.ylim[1], self.xlim[0]:self.xlim[1]]
        my_ints = self.fitopt_dm.yint[0, 0, self.ylim[0]:self.ylim[1], self.xlim[0]:self.xlim[1]]

        for i in range(x):
            for j in range(y):

                slope = my_slopes[j, i]
                yint = my_ints[j, i]

                # convert group to time
                frame_time = (group - 1) * self.ramp_dm.meta.exposure.frame_time
                line = slope * frame_time + yint

                axs.plot(frame_time, my_plot_area[:, j, i], c='b', marker='o', markersize=2.0, linestyle='-',
                         linewidth=0)
                axs.plot(frame_time, line, c='r', marker='h', markersize=0, linestyle='-', linewidth=0.5)

        plt.tight_layout()
        plot_name = 'ramp_fits.pdf'
        try:
            os.remove(plot_name)
        except:
            pass
        fig.savefig(plot_name, dpi=100)

    def plot_hist(self):
        """
        plot a histogram of rate.fits and fitopt results
        """
        fig, axs = plt.subplots(1, 1, figsize=(8, 5))
        axs.set_xlabel('rate (DN/s)', fontsize=15)
        axs.set_ylabel('N', fontsize=15)
        axs.set_title("pixel histograms for {} and {}".format(self.ramp_file, self.fitopt_file), fontsize=12)

        my_slopes_fitopt = self.fitopt_dm.slope[:, :, self.ylim[0]:self.ylim[1], self.xlim[0]:self.xlim[1]]
        my_slopes_rate = self.rate_dm.data[self.ylim[0]:self.ylim[1], self.xlim[0]:self.xlim[1]]

        n, bins, patches = axs.hist(my_slopes_fitopt.ravel(), bins='auto', alpha=0.5, histtype='bar', ec='black',
                                    color='g', label='fitopt.fits')
        axs.hist(my_slopes_rate.ravel(), bins=bins, alpha=0.5, histtype='bar', ec='black',
                 color='b', label='rate.fits')

        axs.legend(prop={'size': 10}, loc=0)
        axs.grid(True, linewidth=0.5)

        plot_name = 'rate_fitopt_histogram.pdf'
        plt.tight_layout()
        try:
            os.remove(plot_name)
        except:
            pass
        fig.savefig(plot_name, dpi=200)

    def plot_pixeldq(self):
        """
        plot the pixel mask histogram and if more than 1 unique value plot the mask itself
        """
        fig, axs = plt.subplots(1, 1, figsize=(8, 4))

        my_ramps_pixdq = self.ramp_dm.pixeldq[self.ylim[0]:self.ylim[1], self.xlim[0]:self.xlim[1]]

        # histogram
        bin_edges = np.arange(-0.5, 2.5, 1.)
        axs.hist(my_ramps_pixdq.ravel(), bins=bin_edges, alpha=0.5, histtype='bar', ec='black', color='b',
                 label='PIXEL_DQ')
        axs.set_xlabel('PIXEL_DQ')
        axs.set_ylabel('N')
        axs.set_title("pixel_dq histogram for {}".format(self.ramp_file), fontsize=12)
        axs.annotate('Unique mask values:', xy=(0.7, 0.95), xycoords='axes fraction', fontsize=12,
                        color='k')
        axs.annotate(np.unique(my_ramps_pixdq), xy=(0.7, 0.89), xycoords='axes fraction', fontsize=12,
                     color='k')

        plt.tight_layout()
        plot_name = 'pixel_mask_histogram.pdf'
        try:
            os.remove(plot_name)
        except:
            pass
        fig.savefig(plot_name, dpi=200)

    def plot_ramps(self):
        """
        use the RampPlots class to plot the ramps in the jump output
        """
        my_ramps = RampPlots(self.ramp_dm, xlim=self.xlim, ylim=self.ylim)
        my_ramps.plot_ramps()

    def plot_groupdq(self):
        """
        use the RampPlots class to plot the groupdq from the jump output
        """
        my_ramps = RampPlots(self.ramp_dm, xlim=self.xlim, ylim=self.ylim)
        my_ramps.plot_groupdq()

    def plot_jump_pixels(self):
        """
        plot the pixels with a jump flagged
        """
        my_slopes_rate = self.rate_dm.data[self.ylim[0]:self.ylim[1], self.xlim[0]:self.xlim[1]]
        my_gdqs = self.ramp_dm.groupdq[0, :, self.ylim[0]:self.ylim[1], self.xlim[0]:self.xlim[1]]
        my_gdqs_image = np.sum(my_gdqs, axis=0)

        fig, axs = plt.subplots(1, 2, figsize=(12, 8))

        axs[0].imshow(my_slopes_rate, cmap='jet', interpolation='nearest', origin='lower')
        axs[1].imshow(my_gdqs_image, cmap='gray', interpolation='nearest', origin='lower', vmin=2, vmax=5)

        axs[0].annotate('rate image', xy=(0.03, 0.97), xycoords='axes fraction', fontsize=12, fontweight='bold',
                        color='w')
        axs[1].annotate('group_dq summed', xy=(0.03, 0.97), xycoords='axes fraction', fontsize=12,
                        fontweight='bold', color='w')

        plt.tight_layout()
        plot_name = 'rate_and_jump_images.pdf'
        try:
            os.remove(plot_name)
        except:
            pass
        fig.savefig(plot_name, dpi=100)

    def plot_all(self):
        """
        run all plot methods
        """
        self.plot_fits()
        self.plot_hist()
        self.plot_pixeldq()
        self.plot_ramps()
        self.plot_groupdq()
        self.plot_jump_pixels()


if __name__ == "__main__":
    # Parse arguments
    help_text = ""
    usage = "\n\n%prog <fitopt file> <ramp file> <rate file>\n"
    usage += "\nMakes a series of plots from the optional fitopt output of ramp_fit"

    parser = optparse.OptionParser(usage)
    (options,args) = parser.parse_args()

    # check for file
    if len(args) == 3:
        try:
            assert os.path.isfile(args[0])

        except AssertionError:
            print(help_text)
            time.sleep(1) # Ensure help text appears before error messages.
            parser.error("fitopt file does not exist...")
            sys.exit(1)

        try:
            assert os.path.isfile(args[1])

        except AssertionError:
            print(help_text)
            time.sleep(1) # Ensure help text appears before error messages.
            parser.error("ramp file does not exist...")
            sys.exit(1)

        try:
            assert os.path.isfile(args[2])

        except AssertionError:
            print(help_text)
            time.sleep(1)  # Ensure help text appears before error messages.
            parser.error("rate file does not exist...")
            sys.exit(1)

    else:
        print(help_text)
        time.sleep(1)  # Ensure help text appears before error messages.
        parser.error("must supply a fitopt, ramp and rate file...")
        sys.exit(1)

    my_ramp_fits = RampFits(args[0], args[1], args[2], xlim=[540, 560], ylim=[640, 660])
    my_ramp_fits.plot_all()
