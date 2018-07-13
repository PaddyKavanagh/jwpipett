# plot_ramps.py

Tools to visualise ramps in a specified region of a MIRIRampModel file. 


# Description
Takes a 4D ramp file and limits in the cols and rows dimensions and plots the ramps
and groupdq vectors for each pixel.


# Command Line

`python plot_ramps.py ramp_file.fits`


# Python Class

## Usage
```python
my_ramps = RampPlots(ramp_file.fits, xlim=[500, 600], ylim=[500, 600])
my_ramps.plot_all()
```

## Class
### class RampPlots():
    """
    Make plots from the a ramp file.

    Parameters
    ----------
    ramp_filename :
        a ramp file from any pipeline step preceding ramp_fit

    xlim : list
        the limits of the xrange to plot

    ylim : list
        the limits of the yrange to plot

    """

## Methods:
### plot_ramps()
	plot the ramp vectors for each pixel
	
### plot_groupdq()   
    plot the groupdq vector for each pixel and the location of pixels flagged with jump/do_not_use

### plot_all()
	run all plotting methods
	