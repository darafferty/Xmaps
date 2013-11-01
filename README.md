Xmaps
=====

Python scripts to make maps of spectral parameters with Chandra data

The make_maps.py script will extract and fit spectra from regions defined 
in an input binmap. The binmap is simply a FITS image with the values that
lie in each bin equal to the bin's number. Binmaps can be generated in a 
number of ways, such as the Contbin algorithm of Sanders (2006) or the WVT 
algorithm of Diehl & Statler (2007).


The make_maps.py script is called as follows:

Usage: make_maps.py [options] <binmap> <evt2_file> <bg_file> <pbk_file> <asol_file> <msk_file> <bpix_file> <redshift> <nH_Gal> <root>

Arguments:
  <binmap>      map of bins, with values equal to bin number
  <evt2_file>   events file or file of list of events files (e.g., @evt2.list)
  <bg_file>     background file or file of list of background files (e.g., @bg.list)
  <pbk_file>    pbk0 file or file of list of pbk0 files (e.g., @pbk0.list)
  <asol_file>   asol1 file or file of list of asol1 files (e.g., @asol1.list).
                If there are more than one asol1 files for an observation,
                list them on one line, separated with a comma and ordered by time
  <msk_file>    msk1 file or file of list of msk1 files (e.g., @msk1.list)
  <bpix_file>   bad pixel file or file of list of msk1 files (e.g., @bpix.list)
  <redshift>    redshift of source
  <nH_Gal>      Galactic N_H (10^22 cm^-2)
  <root>        root name for output map(s)

Options:
  -h, --help            show this help message and exit
  --vars_to_map=VAL     Variable(s) to map (kT, Z, nH, norm, fkT, fZ, fnH, chi2);
                        default = "kT, Z"
  --bin=VAL             Binning for spectra (counts); default = 25
  --kT=VAL              Initial guess for kT (keV); default = 3
  --Ab=VAL              Initial guess for abundance (solar); default = 0.3
  --plindx=VAL          Initial guess for power-law index; default = 1.8
  --lo_energy=VAL       Low energy bound for spectra (keV); default = 0.5
  --hi_energy=VAL       Upper energy bound for spectra (keV); default = 7
  --plasma_model=STR    plasma model to use in fit (mekal or apec); default =
                        mekal
  --fix_nh=VAL          Freeze nH (yes/no); default = yes
  --fix_abund=VAL       Freeze abundance (yes/no); default = no
  --find_err=VAL        Calculate errors (yes/no); default = no
  --binmap_bin=VAL      Binning for binmap (pixels); default = None
  --binmap_minx=VAL     Minimum sky x for binmap (sky coords); default = None
  --binmap_miny=VAL     Minimum sky y for binmap (sky coords); default = None
  -e                    Skip extraction
  -f                    Skip fitting
  -b                    Do binning during extraction instead of during fitting
  -v                    Enable verbose mode
  -p                    Enable plotting
  -c                    Clobber any existing files

A subdirectory in the current directory is made to store the
(many) region files, spectra, responses, and fits. The sub-
directory is named 'root_spectra' where root is specified
in the call to this script. The output maps will be placed
in the current directory and named as follows:

   output maps = root_kT_map.fits and
                 root_Z_map.fits
   fit results = root_spectra/root_wabs_mekal.dat or
                 root_spectra/root_wabs_apec.dat

The input evt2_file, bg_file, and pbk_file can be text files
containing a list of files. If so, they should be prepended by
"@". Note that the WCS of the first observation must match that
of the binmap.

CIAO must be initialized before starting the script.
