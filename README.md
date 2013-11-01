Xmaps
=====

Python scripts to make maps of spectral parameters with Chandra data

The make_maps.py script will extract and fit spectra from regions defined 
in an input binmap. The binmap is simply a FITS image with the values that
lie in each bin equal to the bin's number. Binmaps can be generated in a 
number of ways, such as with the Contbin algorithm of Sanders (2006) or 
the WVT algorithm of Diehl & Statler (2007).

CIAO must be initialized before starting the script.
