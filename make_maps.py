"""
Script to make maps of spectral properties.

Usage: make_maps.py [options] <binmap> <evt2_file> <bg_file> <pbk_file> <asol_file> <msk_file> <redshift> <nH_Gal> <root>

Arguments:
  <binmap>      map of bins, with values equal to bin number
  <evt2_file>   events file or file of list of events files (e.g., @evt2.list)
  <bg_file>     background file or file of list of background files (e.g., @bg.list)
  <pbk_file>    pbk0 file or file of list of pbk0 files (e.g., @pbk0.list)
  <asol_file>   asol1 file or file of list of asol1 files (e.g., @asol1.list).
                If there are more than one asol1 files for an observation,
                list them on one line, separated with a comma and ordered by time
  <msk_file>    msk1 file or file of list of msk1 files (e.g., @msk1.list)
  <redshift>    redshift of source
  <nH_Gal>      Galactic N_H (10^22 cm^-2)
  <root>        root name for output map(s)

Options:
  -h, --help            show this help message and exit
  --vars_to_map=VAL     Variable(s) to map (kT, Z, nH, norm, plindx, mdot, fkT, fZ, fnH, fplindx, fmdot, chi2);
                        default = "kT, Z"
  --add_comp=VAL        Add a component to the default single-temperature model
                        (pow, mekal/apec, mkcflow); default = None'
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
  --binmap_bin=VAL      Binning for binmap (pixels); default = None
  --binmap_minx=VAL     Minimum sky x for binmap (sky coords); default = None
  --binmap_miny=VAL     Minimum sky y for binmap (sky coords); default = None
  -e                    Skip extraction
  -f                    Skip fitting
  -b                    Do binning during extraction instead of during fitting
  -v                    Enable verbose mode
  -p                    Enable plotting
  -c                    Clobber any existing files
  -a                    Add spectra together for fitting (instead of fitting the
                        spectra simultaneously)

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

Version 0.6: 3/5/2011 - Updated for Ciao 4.3
Version 0.7: 1/11/2013 - Updated for Ciao 4.5
Version 0.8: 26/11/2013 - Removed find_err option as redundant; updated to
    work with Python 3

"""

import os
import sys
import pyfits
# Check that CIAO was initialized
if os.environ.get("CALDB") == None:
    sys.exit('Please initalize CIAO before running this script.')
from extract_spectra import *
from fit_spectra import *
import numpy
from sherpa.astro.ui import *
from misc_functions import *


def paint_map(binmap_file, fit_file, vars_to_map, root=None, fit2_file=None, second_comp=None, best_fit=None, Fprob=None, clobber=False):
    """
    Paints the binmap with the values from spectral fitting.

    Inputs:  binmap_file - fits file of map of bins (pixel values = bin numbers)
             fit_file - output of fit_spectra.py with fit results
             vars_to_map - variables to map. A map is created for each variable
             root - root name of ouput map(s); defaults to "output"
             fit2_file - output of fit_spectra.py with fit results for 2-comp model
             second_comp - name of component added to single-temperature model to
                           make the 2-comp model whose results are given in fit2_file
             best_fit - index (1 or 2) of "best" fit; used only when fit2_file and
                        second_comp are specified
             Fprob - F-test probability; used only when fit2_file and
                     second_comp are specified
             clobber - if True, overwrite any existing files

    Outputs: Writes maps using the output of "fit_spectra.py". The maps are called:
             root_kT_map.fits, root_Z_map.fits, etc.

    """
    # Check if min bin is negative or starts or ends on the image boundary.
    # If so, assume it is not wanted (e.g., for wvt bin maps).
    import pyfits
    binmap = pyfits.open(binmap_file, mode="readonly")
    binimage = binmap[0].data
    minbin = int(binimage.min())
    maxbin = int(binimage.max())
    if minbin < 0:
        minbin = 0
    inbin = numpy.where(binimage == minbin)
    if 0 in inbin[0] or numpy.size(binimage,0)-1 in inbin[0]:
        minbin += 1
    nbins = maxbin - minbin + 1

    # Read in the fit results file and
    # calculate errors and check for upper limits
    data1 = read_fit_results(fit_file)
    if 'kT' in data1.dtype.names:
        kT_err = numpy.sqrt(data1['kT_lo']**2 + data1['kT_hi']**2)
        upper_limits = numpy.where(kT_err/data1['kT'] >= 1.0)
        kT_err[upper_limits] = data1['kT'][upper_limits]
    if 'Z' in data1.dtype.names:
        Z_err = numpy.sqrt(data1['Z_lo']**2 + data1['Z_hi']**2)
        upper_limits = numpy.where(Z_err/data1['Z'] >= 1.0)
        Z_err[upper_limits] = data1['Z'][upper_limits]
    if 'nH' in data1.dtype.names:
        nH_err = numpy.sqrt(data1['nH_lo']**2 + data1['nH_hi']**2)
        upper_limits = numpy.where(nH_err/data1['nH'] >= 1.0)
        nH_err[upper_limits] = data1['nH'][upper_limits]
    if 'norm' in data1.dtype.names:
        norm_err = numpy.sqrt(data1['norm_lo']**2 + data1['norm_hi']**2)
        upper_limits = numpy.where(norm_err/data1['norm'] >= 1.0)
        norm_err[upper_limits] = data1['norm'][upper_limits]

    if fit2_file != None:
        data2 = read_fit_results(fit2_file, second_comp=second_comp)
        if 'kT' in data2.dtype.names:
            kT2_err = numpy.sqrt(data2['kT_lo']**2 + data2['kT_hi']**2)
            upper_limits = numpy.where(kT2_err/data2['kT'] >= 1.0)
            kT2_err[upper_limits] = data2['kT'][upper_limits]
        if 'kT1' in data2.dtype.names:
            kT1_err = numpy.sqrt(data2['kT1_lo']**2 + data2['kT1_hi']**2)
            upper_limits = numpy.where(kT1_err/data2['kT1'] >= 1.0)
            kT1_err[upper_limits] = data2['kT1'][upper_limits]
        if 'kT2' in data2.dtype.names:
            kT2_err = numpy.sqrt(data2['kT2_lo']**2 + data2['kT2_hi']**2)
            upper_limits = numpy.where(kT2_err/data2['kT2'] >= 1.0)
            kT2_err[upper_limits] = data2['kT2'][upper_limits]
        if 'Z' in data2.dtype.names:
            Z2_err = numpy.sqrt(data2['Z_lo']**2 + data2['Z_hi']**2)
            upper_limits = numpy.where(Z2_err/data2['Z'] >= 1.0)
            Z2_err[upper_limits] = data2['Z'][upper_limits]
        if 'Z1' in data2.dtype.names:
            Z1_err = numpy.sqrt(data2['Z1_lo']**2 + data2['Z1_hi']**2)
            upper_limits = numpy.where(Z1_err/data2['Z1'] >= 1.0)
            Z1_err[upper_limits] = data2['Z1'][upper_limits]
        if 'Z2' in data2.dtype.names:
            Z2_err = numpy.sqrt(data2['Z2_lo']**2 + data2['Z2_hi']**2)
            upper_limits = numpy.where(Z2_err/data2['Z2'] >= 1.0)
            Z2_err[upper_limits] = data2['Z2'][upper_limits]
        if 'nH' in data2.dtype.names:
            nH2_err = numpy.sqrt(data2['nH_lo']**2 + data2['nH_hi']**2)
            upper_limits = numpy.where(nH2_err/data2['nH'] >= 1.0)
            nH2_err[upper_limits] = data2['nH'][upper_limits]
        if 'norm' in data2.dtype.names:
            norm2_err = numpy.sqrt(data2['norm_lo']**2 + data2['norm_hi']**2)
            upper_limits = numpy.where(norm2_err/data2['norm'] >= 1.0)
            norm2_err[upper_limits] = data2['norm'][upper_limits]
        if 'norm1' in data2.dtype.names:
            norm1_err = numpy.sqrt(data2['norm1_lo']**2 + data2['norm1_hi']**2)
            upper_limits = numpy.where(norm1_err/data2['norm1'] >= 1.0)
            norm1_err[upper_limits] = data2['norm1'][upper_limits]
        if 'norm2' in data2.dtype.names:
            norm2_err = numpy.sqrt(data2['norm2_lo']**2 + data2['norm2_hi']**2)
            upper_limits = numpy.where(norm2_err/data2['norm2'] >= 1.0)
            norm2_err[upper_limits] = data2['norm2'][upper_limits]
        if 'plindx' in data2.dtype.names:
            plindx_err = numpy.sqrt(data2['plindx_lo']**2 + data2['plindx_hi']**2)
        if 'mdot' in data2.dtype.names:
            mkcnorm_err = numpy.sqrt(data2['mkcnorm_lo']**2 + data2['mkcnorm_hi']**2)
            upper_limits = numpy.where(mkcnorm_err/data2['mkcnorm'] >= 1.0)
            mkcnorm_err[upper_limits] = data2['mkcnorm'][upper_limits]

    nreg1 = len(data1)
    if fit2_file != None:
        nreg2 = len(data2)
    else:
        nreg2 = nreg1

    # Make sure both have the same length, and they match the number of bins in the binmap
    if nreg1 != nreg2:
        sys.exit('ERROR: The two fits have a different number of regions. Please check the fitting results')
    if nreg1 != nbins or nreg2 != nbins:
        sys.exit('ERROR: Number of regions does not match the number of bins. Please check the fitting results')

    if best_fit == None:
        best_fit = numpy.ones(nreg1, dtype=int)

    # import pdb; pdb.set_trace()
    # make copies of the binmap as needed
    if 'kT' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binimage_kT1 = numpy.zeros(binimage.shape, dtype=float)
            binimage_kT2 = numpy.zeros(binimage.shape, dtype=float)
        else:
            binimage_kT = numpy.zeros(binimage.shape, dtype=float)
    if 'Z' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binimage_Z1 = numpy.zeros(binimage.shape, dtype=float)
            binimage_Z2 = numpy.zeros(binimage.shape, dtype=float)
        else:
            binimage_Z = numpy.zeros(binimage.shape, dtype=float)
    if 'plindx' in vars_to_map: binimage_plindx = numpy.zeros(binimage.shape, dtype=float)
    if 'mdot' in vars_to_map: binimage_mkcnorm = numpy.zeros(binimage.shape, dtype=float)
    if 'nH' in vars_to_map: binimage_nH = numpy.zeros(binimage.shape, dtype=float)
    if 'norm' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binimage_norm1 = numpy.zeros(binimage.shape, dtype=float)
            binimage_norm2 = numpy.zeros(binimage.shape, dtype=float)
        else:
            binimage_norm = numpy.zeros(binimage.shape, dtype=float)
    if 'fkT' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binimage_fkT1 = numpy.zeros(binimage.shape ,dtype=float)
            binimage_fkT2 = numpy.zeros(binimage.shape, dtype=float)
        else:
            binimage_fkT = numpy.zeros(binimage.shape, dtype=float)
    if 'fZ' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binimage_fZ1 = numpy.zeros(binimage.shape, dtype=float)
            binimage_fZ2 = numpy.zeros(binimage.shape, dtype=float)
        else:
            binimage_fZ = numpy.zeros(binimage.shape, dtype=float)
    if 'fplindx' in vars_to_map: binimage_fplindx = numpy.zeros(binimage.shape, dtype=float)
    if 'fmdot' in vars_to_map: binimage_fmkcnorm = numpy.zeros(binimage.shape, dtype=float)
    if 'fnH' in vars_to_map: binimage_fnH = numpy.zeros(binimage.shape, dtype=float)
    if 'fnorm' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binimage_fnorm1 = numpy.zeros(binimage.shape, dtype=float)
            binimage_fnorm2 = numpy.zeros(binimage.shape, dtype=float)
        else:
            binimage_fnorm = numpy.zeros(binimage.shape, dtype=float)
    if 'chi2' in vars_to_map: binimage_chi2 = numpy.zeros(binimage.shape, dtype=float)
    if Fprob != None: binimage_Fprob = numpy.zeros(binimage.shape, dtype=float)

    for k in range(nreg1):
        # First, make sure the loop index matches the reg_id of the region of interest (i.e. to catch entries that are out order or missing)
        if k + minbin == data1['reg_id'][k]:
            i = k
        else:
            i = data1['reg_id'][k]

        inbin = numpy.where(binimage == i + minbin) # find all pixels in region of interest
        if 'kT' in vars_to_map:
            if second_comp == 'mekal' or second_comp == 'apec':
                if data2['kT1'][i] <= data2['kT2'][i]:
                    binimage_kT1[inbin] = data2['kT1'][i]
                    binimage_kT2[inbin] = data2['kT2'][i]
                else:
                    binimage_kT1[inbin] = data2['kT2'][i]
                    binimage_kT2[inbin] = data2['kT1'][i]
                if best_fit[i] == 1: # single-temp model preferred
                    binimage_kT1[inbin] = data1['kT'][i]
                    binimage_kT2[inbin] = data1['kT'][i]
            else:
                if best_fit[i] == 1:
                    binimage_kT[inbin] = data1['kT'][i]
                else:
                    binimage_kT[inbin] = data2['kT'][i]
        if 'Z' in vars_to_map:
            if second_comp == 'mekal' or second_comp == 'apec':
                if data2['kT1'][i] <= data2['kT2'][i]:
                    binimage_Z1[inbin] = data2['Z1'][i]
                    binimage_Z2[inbin] = data2['Z2'][i]
                else:
                    binimage_Z1[inbin] = data2['Z2'][i]
                    binimage_Z2[inbin] = data2['Z1'][i]
                if best_fit[i] == 1:
                    binimage_Z1[inbin] = data1['Z'][i]
                    binimage_Z2[inbin] = data1['Z'][i]
            else:
                if best_fit[i] == 1:
                    binimage_Z[inbin] = data1['Z'][i]
                else:
                    binimage_Z[inbin] = data2['Z'][i]
        if 'plindx' in vars_to_map:
            if best_fit[i] == 1:
                binimage_plindx[inbin] = 0.0
            else:
                binimage_plindx[inbin] = data2['plindx'][i]
        if 'mdot' in vars_to_map:
            if best_fit[i] == 1:
                binimage_mkcnorm[inbin] = 0.0
            else:
                binimage_mkcnorm[inbin] = data2['mkcnorm'][i]
        if 'nH' in vars_to_map:
            if best_fit[i] == 1:
                binimage_nH[inbin] = data1['nH'][i]
            else:
                binimage_nH[inbin] = data2['nH'][i]
        if 'norm' in vars_to_map:
            if second_comp == 'mekal' or second_comp == 'apec':
                if data2['kT1'][i] <= data2['kT2'][i]:
                    binimage_norm1[inbin] = data2['norm1'][i]
                    binimage_norm2[inbin] = data2['norm2'][i]
                else:
                    binimage_norm1[inbin] = data2['norm2'][i]
                    binimage_norm2[inbin] = data2['norm1'][i]
                if best_fit[i] == 1: # single-temp model preferred
                    binimage_norm1[inbin] = data1['norm'][i]
                    binimage_norm2[inbin] = data1['norm'][i]
            else:
                if best_fit[i] == 1:
                    binimage_norm[inbin] = data1['norm'][i]
                else:
                    binimage_norm[inbin] = data2['norm'][i]
        if 'fkT' in vars_to_map:
            if second_comp == 'mekal' or second_comp == 'apec':
                if data2['kT1'][i] <= data2['kT2'][i]:
                    binimage_fkT1[inbin] = kT1_err[i]
                    binimage_fkT2[inbin] = kT2_err[i]
                else:
                    binimage_fkT1[inbin] = kT2_err[i]
                    binimage_fkT2[inbin] = kT1_err[i]
                if best_fit[i] == 1: # single-temp model preferred
                    binimage_fkT1[inbin] = kT_err[i]
                    binimage_fkT2[inbin] = kT_err[i]
            else:
                if best_fit[i] == 1:
                    binimage_fkT[inbin] = kT_err[i]
                else:
                    binimage_fkT[inbin] = kT2_err[i]
        if 'fZ' in vars_to_map:
            if second_comp == 'mekal' or second_comp == 'apec':
                if data2['kT1'][i] <= data2['kT2'][i]:
                    binimage_fZ1[inbin] = Z1_err[i]
                    binimage_fZ2[inbin] = Z2_err[i]
                else:
                    binimage_fZ1[inbin] = Z2_err[i]
                    binimage_fZ2[inbin] = Z1_err[i]
                if best_fit[i] == 1: # single-temp model preferred
                    binimage_fZ1[inbin] = Z_err[i]
                    binimage_fZ2[inbin] = Z_err[i]
            else:
                if best_fit[i] == 1:
                    binimage_fZ[inbin] = Z_err[i]
                else:
                    binimage_fZ[inbin] = Z2_err[i]
        if 'fplindx' in vars_to_map:
            if best_fit[i] == 1:
                binimage_fplindx[inbin] = 0.0
            else:
                binimage_fplindx[inbin] = plindx_err[i]
        if 'fmdot' in vars_to_map:
            if best_fit[i] == 1:
                binimage_fmkcnorm[inbin] = 0.0
            else:
                binimage_fmkcnorm[inbin] = mkcnorm_err[i]
        if 'fnH' in vars_to_map:
            if best_fit[i] == 1:
                binimage_fnH[inbin] = nH_err[i]
            else:
                binimage_fnH[inbin] = nH2_err[i]
        if 'fnorm' in vars_to_map:
            if second_comp == 'mekal' or second_comp == 'apec':
                if data2['norm1'][i] <= data2['norm2'][i]:
                    binimage_fnorm1[inbin] = norm1_err[i]
                    binimage_fnorm2[inbin] = norm2_err[i]
                else:
                    binimage_fnorm1[inbin] = norm2_err[i]
                    binimage_fnorm2[inbin] = norm1_err[i]
                if best_fit[i] == 1: # single-temp model preferred
                    binimage_fnorm1[inbin] = norm_err[i]
                    binimage_fnorm2[inbin] = norm_err[i]
            else:
                if best_fit[i] == 1:
                    binimage_fnorm[inbin] = norm_err[i]
                else:
                    binimage_fnorm[inbin] = norm2_err[i]
        if 'chi2' in vars_to_map:
            if best_fit[i] == 1:
                binimage_chi2[inbin] = data1['chi2'][i]
            else:
                binimage_chi2[inbin] = data2['chi2'][i]
        if Fprob != None:
            binimage_Fprob[inbin] = Fprob[i]

    if root == None:
        root = 'output'
    if 'kT' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binmap[0].data = binimage_kT1
            binmap.writeto(root+'_kT1_map.fits', clobber=clobber)
            binmap[0].data = binimage_kT2
            binmap.writeto(root+'_kT2_map.fits', clobber=clobber)
        else:
            binmap[0].data = binimage_kT
            binmap.writeto(root+'_kT_map.fits', clobber=clobber)
    if 'Z' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binmap[0].data = binimage_Z1
            binmap.writeto(root+'_Z1_map.fits', clobber=clobber)
            binmap[0].data = binimage_Z2
            binmap.writeto(root+'_Z2_map.fits', clobber=clobber)
        else:
            binmap[0].data = binimage_Z
            binmap.writeto(root+'_Z_map.fits', clobber=clobber)
    if 'plindx' in vars_to_map:
        binmap[0].data = binimage_plindx
        binmap.writeto(root+'_plindx_map.fits', clobber=clobber)
    if 'mdot' in vars_to_map:
        binmap[0].data = binimage_mkcnorm
        binmap.writeto(root+'_mdot_map.fits', clobber=clobber)
    if 'nH' in vars_to_map:
        binmap[0].data = binimage_nH
        binmap.writeto(root+'_nH_map.fits', clobber=clobber)
    if 'norm' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binmap[0].data = binimage_norm1
            binmap.writeto(root+'_norm1_map.fits', clobber=clobber)
            binmap[0].data = binimage_norm2
            binmap.writeto(root+'_norm2_map.fits', clobber=clobber)
        else:
            binmap[0].data = binimage_norm
            binmap.writeto(root+'_norm_map.fits', clobber=clobber)
    if 'fkT' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binmap[0].data = binimage_fkT1
            binmap.writeto(root+'_fkT1_map.fits', clobber=clobber)
            binmap[0].data = binimage_fkT2
            binmap.writeto(root+'_fkT2_map.fits', clobber=clobber)
        else:
            binmap[0].data = binimage_fkT
            binmap.writeto(root+'_fkT_map.fits', clobber=clobber)
    if 'fZ' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binmap[0].data = binimage_fZ1
            binmap.writeto(root+'_fZ1_map.fits', clobber=clobber)
            binmap[0].data = binimage_fZ2
            binmap.writeto(root+'_fZ2_map.fits', clobber=clobber)
        else:
            binmap[0].data = binimage_fZ
            binmap.writeto(root+'_fZ_map.fits', clobber=clobber)
    if 'fnorm' in vars_to_map:
        if second_comp == 'mekal' or second_comp == 'apec':
            binmap[0].data = binimage_fnorm1
            binmap.writeto(root+'_fnorm1_map.fits', clobber=clobber)
            binmap[0].data = binimage_fnorm2
            binmap.writeto(root+'_fnorm2_map.fits', clobber=clobber)
        else:
            binmap[0].data = binimage_fnorm
            binmap.writeto(root+'_fnorm_map.fits', clobber=clobber)
    if 'fplindx' in vars_to_map:
        binmap[0].data = binimage_fplindx
        binmap.writeto(root+'_fplindx_map.fits', clobber=clobber)
    if 'fmdot' in vars_to_map:
        binmap[0].data = binimage_fmkcnorm
        binmap.writeto(root+'_fmdot_map.fits', clobber=clobber)
    if 'fnH' in vars_to_map:
        binmap[0].data = binimage_fnH
        binmap.writeto(root+'_fnH_map.fits', clobber=clobber)
    if 'chi2' in vars_to_map:
        binmap[0].data = binimage_chi2
        binmap.writeto(root+'_chi2_map.fits', clobber=clobber)
    if Fprob != None:
        binmap[0].data = binimage_Fprob
        binmap.writeto(root+'_Ftest_map.fits', clobber=clobber)


def read_fit_results(filename, second_comp=None):
    """
    Reads in the results output by fit_spectra.py and returns a dictionary.

    Inputs:  filename - output of fit_spectra.py with fit results.
             second_comp - name of second component use in fit ('mekal', 'apec', 'pow', 'mkcflow')

    Outputs: Returns dictionary of fit results.

    """
    # Read in the fit results file, given the format implied by second_comp
    if second_comp == None:
        dtype = {'names': ('reg_id', 'kT', 'kT_lo', 'kT_hi', 'Z', 'Z_lo', 'Z_hi', 'norm', 'norm_lo', 'norm_hi', 'nH', 'nH_lo', 'nH_hi', 'chi2', 'totcnts', 'nbins'), 'formats': ('i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4')}
    else:
        if second_comp == 'mekal' or second_comp == 'apec':
            dtype = {'names': ('reg_id', 'kT1', 'kT1_lo', 'kT1_hi', 'Z1', 'Z1_lo', 'Z1_hi', 'norm1', 'norm1_lo', 'norm1_hi', 'kT2', 'kT2_lo', 'kT2_hi', 'Z2', 'Z2_lo', 'Z2_hi', 'norm2', 'norm2_lo', 'norm2_hi', 'nH', 'nH_lo', 'nH_hi', 'chi2', 'totcnts', 'nbins'), 'formats': ('i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4')}
        if second_comp == 'pow':
            dtype = {'names': ('reg_id', 'kT', 'kT_lo', 'kT_hi', 'Z', 'Z_lo', 'Z_hi', 'norm', 'norm_lo', 'norm_hi', 'plindx', 'plindx_lo', 'plindx_hi', 'plnorm', 'plnorm_lo', 'plnorm_hi', 'nH', 'nH_lo', 'nH_hi', 'chi2', 'totcnts', 'nbins'), 'formats': ('i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4')}
        if second_comp == 'mkcflow':
            dtype = {'names': ('reg_id', 'kT', 'kT_lo', 'kT_hi', 'Z', 'Z_lo', 'Z_hi', 'norm', 'norm_lo', 'norm_hi', 'mkcnorm', 'mkcnorm_lo', 'mkcnorm_hi', 'nH', 'nH_lo', 'nH_hi', 'chi2', 'totcnts', 'nbins'), 'formats': ('i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4')}

    data = numpy.loadtxt(filename, dtype=dtype)
    return data


def apply_binning_to_image(binmap_file, image_file, root=None, clobber=False):
    """
    Applies the binning scheme from the binmap to an image of the same shape.

    Inputs:  binmap_file - fits file of map of bins (pixel values = bin numbers)
             image_file - fits file of image (must have the same shape as the binmap)
             root - root name of ouput map; defaults to image_file_root + "_binned"
             clobber - if True, overwrite any existing files

    Outputs: Binned version of the input image with each bin set to the mean value
             inside that bin, named "root.fits".

    """
    if root == None:
        root = os.path.splitext(image_file)[0] + '_binned'

    # Check if min bin is negative or starts or ends on the image boundary.
    # If so, assume it is not wanted (e.g., for wvt bin maps).
    import pyfits
    binmap = pyfits.open(binmap_file)
    binimage = binmap[0].data
    minbin = int(binimage.min())
    maxbin = int(binimage.max())
    if minbin < 0:
        minbin = 0
    inbin = numpy.where(binimage == minbin)
    if 0 in inbin[0] or numpy.size(binimage,0)-1 in inbin[0]:
        minbin += 1
    nbins = maxbin - minbin + 1
    image = pyfits.open(image_file)
    im = image[0].data

    # Check that the binmap and image have the same shape
    if binimage.shape != im.shape:
        sys.exit('ERROR: Input binmap and image must have the same shape.')

    # make copy of the binmap
    binimage_out = binimage.astype(float)

    for i in range(nbins):
        inbin = numpy.where(binimage == i + minbin)
        binimage_out[inbin] = numpy.mean(im[inbin])

    binmap[0].data = binimage_out
    binmap.writeto(root+'.fits', clobber=clobber)


def compare_fits(fit1_file, fit2_file, second_comp, null_prob=0.05):
    """
    Performs an F-test on two fits and returns index of better fit (1 or 2).

    Inputs:  fit1_file - output of fit_spectra for fit 1.
             fit2_file - output of fit_spectra for fit 2.
             second_comp - name of component added to model 1 to get model 2
                           (pow, mekal/apec, mkcflow).
             null_prob - probability of null hypothesis below which model 2 is
                         prefered.

    Outputs: Returns 1 or 2 for each region in the input fit files, where 1
             indicates that the second fit is NOT significantly better than the
             first one and 2 indicates that it is.

    """
    # Read in chi2 and num_bins for each fit
    data1 = read_fit_results(fit1_file)
    data2 = read_fit_results(fit2_file, second_comp=second_comp)
    p1 = 3
    if second_comp == 'mekal' or second_comp == 'apec':
        p2 = 6
    if second_comp == 'pow':
        p2 = 5
    if second_comp == 'mkcflow':
        p2 = 4

    nreg = len(data1)
    index_of_pref = []
    Fprob_list = []
    for i in range(nreg):
        dof1 = data1['nbins'][i] - p1
        dof2 = data2['nbins'][i] - p2
        stat1 =  data1['chi2'][i] * dof1 # chi^2
        stat2 =  data2['chi2'][i] * dof2
        Fprob = calc_ftest(dof1, stat1, dof2, stat2)
        Fprob_list.append(Fprob)
        # print 'fit #: '+str(i+1)+',  F-test prob: '+str(Fprob)
        if Fprob <= null_prob:
            index_of_pref.append(2)
        else:
            index_of_pref.append(1)
    return index_of_pref, Fprob_list

if __name__=='__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <binmap> <evt2_file> <bg_file> <pbk_file> <asol_file> <msk_file> <bpix_file> <redshift> <nH_Gal> <root>\n\nArguments:\n  <binmap>      map of bins, with values equal to bin number\n  <evt2_file>   events file or file of list of events files (e.g., @evt2.list)\n  <bg_file>     background file or file of list of background files (e.g., @bg.list)\n  <pbk_file>    pbk0 file or file of list of pbk0 files (e.g., @pbk0.list)\n  <asol_file>   asol1 file or file of list of asol1 files (e.g., @asol1.list).\n                If there are more than one asol1 files for an observation,\n                list them on one line, separated with a comma and ordered by time\n  <msk_file>    msk1 file or file of list of msk1 files (e.g., @msk1.list)\n  <bpix_file>    bad pixel file or file of list of bad pixel files (e.g., @bpix.list)\n  <redshift>    redshift of source\n  <nH_Gal>      Galactic N_H (10^22 cm^-2)\n  <root>        root of output map(s)', version="%prog 0.6")
    parser.add_option('--vars_to_map', dest='vars_to_map', help='Variable(s) to map (kT, Z, nH, norm, plindx, mdot, fkT, fZ, fnH, fnorm, fplindx, fmdot, chi2); default = "kT, Z"', metavar='VAL', default='kT, Z')
    parser.add_option('--add_comp', dest='second_comp', help='Add a component to the default single-temperature model (pow, mekal/apec, mkcflow); default = None', metavar='VAL', default=None)
    parser.add_option('--bin', dest='binning', help='Binning for spectra (counts); default = 25', metavar='VAL', default='25')
    parser.add_option('--kT', dest='kT_guess', help='Initial guess for kT (keV); default = 3', metavar='VAL', default='3.0')
    parser.add_option('--Ab', dest='Ab_guess', help='Initial guess for abundance (solar); default = 0.3', metavar='VAL', default='0.3')
    parser.add_option('--plindx', dest='plindx_guess', help='Initial guess for power-law index; default = 1.8', metavar='VAL', default='1.8')
    parser.add_option('--lo_energy', dest='lo_energy', help='Low energy bound for spectra (keV); default = 0.5', metavar='VAL', default='0.5')
    parser.add_option('--hi_energy', dest='hi_energy', help='Upper energy bound for spectra (keV); default = 7', metavar='VAL', default='7.0')
    parser.add_option('--plasma_model', dest='plasma_model', help='plasma model to use in fit (mekal or apec); default = mekal', metavar='STR', default='mekal')
    parser.add_option('--fix_nh', dest='fix_nH', help='Freeze nH (yes/no); default = yes', metavar='VAL', default='yes')
    parser.add_option('--fix_abund', dest='fix_abund', help='Freeze abundance (yes/no); default = no', metavar='VAL', default='no')
    parser.add_option('--binmap_bin', dest='binmap_binning', help='Binning for binmap (pixels); default = None', metavar='VAL', default=None)
    parser.add_option('--min_rate', dest='min_cnt_rate_ratio', help='Minimum count rate ratio (relative to maximum count rate in region) below which observations are rejected; default = 0.3', metavar='VAL', default='0.3')
    parser.add_option('--binmap_minx', dest='binmap_minx', help='Minimum sky x for binmap (sky coords); default = None', metavar='VAL', default=None)
    parser.add_option('--binmap_miny', dest='binmap_miny', help='Minimum sky y for binmap (sky coords); default = None', metavar='VAL', default=None)
    parser.add_option('-e', action='store_true', dest='skip_extract', help='Skip extraction', default=False)
    parser.add_option('-f', action='store_true', dest='skip_fit', help='Skip fitting', default=False)
    parser.add_option('-b', action='store_true', dest='bin_in_extract', help='Do binning during extraction instead of during fitting', default=False)
    parser.add_option('-v', action='store_true', dest='verbose', help='Enable verbose mode', default=False)
    parser.add_option('-p', action='store_true', dest='plot', help='Enable plotting of spectral fits', default=False)
    parser.add_option('-a', action='store_true', dest='add_spectra', help='Add spectra together for fitting (instead of fitting the spectra simultaneously)', default=False)
    parser.add_option('-c', action='store_true', dest='clobber', help='Clobber any existing files', default=False)
    (options, args) = parser.parse_args()
    if len(args) == 10:
        binmap = args[0]
        evt2_file = args[1]
        bg_file = args[2]
        pbk_file = args[3]
        asol_file = args[4]
        msk_file = args[5]
        bpix_file = args[6]
        redshift = args[7]
        nH_Gal = args[8]
        root = args[9]
        v_to_map = options.vars_to_map
        if ',' in v_to_map:
            vars_to_map = v_to_map.split(',')
            for i in range(len(vars_to_map)):
                vars_to_map[i] = vars_to_map[i].strip()
        else:
            vars_to_map = v_to_map.split()
        second_comp = options.second_comp
        allowed_vars_to_map = ['kT', 'Z', 'nH', 'norm', 'plindx', 'mdot', 'fkT', 'fZ', 'fnH', 'fnorm', 'fplindx', 'fmdot', 'chi2']
        for var in vars_to_map:
            if var in allowed_vars_to_map:
                if var == 'plindx' and second_comp != 'pow':
                    sys.exit('ERROR: Map variable "'+var+'" allowed only with --add_comp=pow.')
                if var == 'mdot' and second_comp != 'mkcflow':
                    sys.exit('ERROR: Map variable "'+var+'" allowed only with --add_comp=mkcflow.')
            else:
                sys.exit('ERROR: Map variable "'+var+'" not allowed.')
        binning = options.binning
        kT_guess = options.kT_guess
        Ab_guess = options.Ab_guess
        plindx_guess = options.plindx_guess
        lo_energy = options.lo_energy
        hi_energy = options.hi_energy
        plasma_model = options.plasma_model
        fix_nH = options.fix_nH
        if fix_nH == 'yes':
            fix_nH_Gal = True
        else:
            fix_nH_Gal = False
        fix_abundance = options.fix_abund
        if fix_abundance == 'yes':
            fix_abund = True
        else:
            fix_abund = False
        binmap_bin = options.binmap_binning
        if binmap_bin != None:
            binmap_bin = int(binmap_bin)
        binmap_minx = options.binmap_minx
        if binmap_minx != None:
            binmap_minx = float(binmap_minx)
        binmap_miny = options.binmap_miny
        if binmap_miny != None:
            binmap_miny = float(binmap_miny)
        min_cnt_rate_ratio = options.min_cnt_rate_ratio
        if min_cnt_rate_ratio != None:
            min_cnt_rate_ratio = float(min_cnt_rate_ratio)
        skip_extract = options.skip_extract
        skip_fit = options.skip_fit
        bin_in_extract = options.bin_in_extract
        if bin_in_extract:
            if add_spectra:
                sys.exit('Binning must be done during fitting if spectra are combined')
            binning_extract = int(binning)
            binning_fit = None
        else:
            binning_extract = None
            binning_fit = int(binning)
        verbose = options.verbose
        make_plots = options.plot
        add_spectra = options.add_spectra
        quiet = not verbose
        clobber = options.clobber

        # Read evt2 file names from the evt2_file if evt2_file begins with '@'
        # Do the same for bg_file, pbk_file, asol_file, and msk_file. For the
        # the asol_file, which may be a list itself (e.g., "asol1, asol2"), we
        # need to split using "," as the separator and make it a list of lists.
        evt2_list = stack_to_list(evt2_file, adjust_path=True)
        bg_list = stack_to_list(bg_file, adjust_path=True)
        pbk_list = stack_to_list(pbk_file, adjust_path=True)
        asol_list = stack_to_list(asol_file, adjust_path=True, stack_of_stacks=True)
        msk_list = stack_to_list(msk_file, adjust_path=True)
        bpix_list = stack_to_list(bpix_file, adjust_path=True)

        # Check that each observation has an evt2, bg, pbk, asol, and msk file. Also,
        # if the files do not have absolute paths, prepend "../" to each file
        # so that their relative paths work from the spectra subdirectory.
        nobs = len(evt2_list)
        if nobs != len(bg_list) or nobs != len(pbk_list) or nobs != len(msk_list) or nobs != len(asol_list) or nobs != len(bpix_list):
            sys.exit('ERROR: You must give the same number of evt2, bg, pbk, and msk files and at least as many asol files as evt2 files.')

        # Check, if clobber = False, that output map files do not already exist
        if clobber == False:
            for var in vars_to_map:
                if os.path.isfile(root+'_'+var+'_map.fits'):
                    sys.exit('ERROR: Output file '+root+'_'+var+'_map.fits'+' exists and clobber = False.')

        # Make regions from the input binmap
        if not skip_extract:
            print('\nDetermining regions from binmap...')
            region_list_binmap = make_regions_from_binmap(binmap, root+'_spectra', minx=binmap_minx, miny=binmap_miny, bin=binmap_bin, skip_dmimglasso=False, clobber=clobber)
        if skip_extract and not skip_fit:
            print('\nDetermining regions from binmap...')
            region_list_binmap = make_regions_from_binmap(binmap, root+'_spectra', minx=binmap_minx, miny=binmap_miny, bin=binmap_bin, skip_dmimglasso=True, clobber=clobber)

        # Extract the spectra and responses
        if not skip_extract:
            os.chdir(root+'_spectra')
            hdr_obs1 = pyfits.getheader(evt2_list[0], 1)
            for j in range(nobs):
                if j == 0:
                    region_list = transform_regions(region_list_binmap, hdr_obs1, hdr_obs1, 'obs'+str(j+1)+'_', clobber=clobber)
                else:
                    hdr = pyfits.getheader(evt2_list[j], 1)
                    region_list = transform_regions(region_list_binmap, hdr_obs1, hdr, 'obs'+str(j+1)+'_', clobber=clobber)

                wextract(region_list, evt2_list[j], pbk_list[j], asol_list[j], msk_list[j], bg_file=bg_list[j], bpix_file=bpix_list[j], binning=binning_extract, quiet=quiet, clobber=clobber)
            os.chdir('..')

        # Fit the spectra
        if not skip_fit:
            if bin_in_extract:
                sp_file_suffix = '_grp.pi'
            else:
                sp_file_suffix = '.pi'
            for j in range(nobs):
                spectra_list_append = []
                for i in range(len(region_list_binmap)):
                    pi_file = 'obs'+str(j+1)+'_' + os.path.splitext(region_list_binmap[i])[0] + '_sou' + sp_file_suffix
                    spectra_list_append.append(pi_file)
                if j==0:
                    spectra_list = [spectra_list_append]
                else:
                    spectra_list.append(spectra_list_append)

            if add_spectra:
                print('\Combining spectra in each region...')
                os.chdir(root+'_spectra')
                spectra_list = combine_spectra(spectra_list, 'combined', clobber=clobber)
                os.chdir('..')

            os.chdir(root+'_spectra')
            if 'fkT' in vars_to_map or 'fZ' in vars_to_map or 'fnH' in vars_to_map or 'fnorm' in vars_to_map or 'fplindx' in vars_to_map or 'fmdot' in vars_to_map:
                find_errors = True
            else:
                find_errors = False

            if 'nH' in vars_to_map or 'fnH' in vars_to_map:
                fix_nH_Gal = False
            if 'Z' in vars_to_map or 'fZ' in vars_to_map:
                fix_abund = False

            first_reg_num = int(os.path.splitext(region_list_binmap[0])[0][-1:])
            call_sherpa_1T(spectra_list, redshift, nH_Gal, kT_guess, Ab_guess, root, lo_energy=lo_energy, hi_energy=hi_energy, plasma_model=plasma_model, fix_nH_Gal=fix_nH_Gal, find_errors=find_errors, binning=binning_fit, fix_abund=fix_abund, reg_num_to_start=first_reg_num, clobber=clobber, make_plots=make_plots, min_cnt_rate_ratio=min_cnt_rate_ratio)

            # If a second component is specified, fit again with it
            # included and compare chi2 for each region. If the fit
            # for a given region improves significantly with the
            # second component included, use the two-component fit
            # results in the maps; otherwise, use the single-
            # temperature results.
            if second_comp != None:
                if second_comp == 'pow':
                    call_sherpa_1T_plus_pow(spectra_list, redshift, nH_Gal, kT_guess, Ab_guess, plindx_guess, root, lo_energy=lo_energy, hi_energy=hi_energy, plasma_model=plasma_model, fix_nH_Gal=fix_nH_Gal, find_errors=find_errors, binning=binning_fit, reg_num_to_start=first_reg_num, clobber=clobber, make_plots=make_plots, fix_abund=fix_abund)
                    fit2_file = root + '_wabs_' + plasma_model + '_pow.dat'
                if second_comp == 'mekal' or second_comp == 'apec':
                    call_sherpa_2T(spectra_list, redshift, nH_Gal, kT_guess, Ab_guess, root, lo_energy=lo_energy, hi_energy=hi_energy, plasma_model=plasma_model, fix_nH_Gal=fix_nH_Gal, find_errors=find_errors, binning=binning_fit, reg_num_to_start=first_reg_num, clobber=clobber, make_plots=make_plots, fix_abund=fix_abund)
                    fit2_file = root + '_wabs_2' + plasma_model + '.dat'
                if second_comp == 'mkcflow':
                    call_sherpa_1T_plus_mkcflow(spectra_list, redshift, nH_Gal, kT_guess, Ab_guess, root, lo_energy=lo_energy, hi_energy=hi_energy, plasma_model=plasma_model, fix_nH_Gal=fix_nH_Gal, find_errors=find_errors, binning=binning_fit, reg_num_to_start=first_reg_num, clobber=clobber, make_plots=make_plots, fix_abund=fix_abund)
                    fit2_file = root + '_wabs_' + plasma_model + '_mkcflow.dat'

            os.chdir('..')

        # Make the map
        fit_file = root + '_wabs_' + plasma_model + '.dat'
        if os.path.isfile(root+'_spectra/'+fit_file) == False:
            sys.exit('ERROR: you must perform spectral fitting before mapping.')
        if second_comp != None:
            if second_comp == 'pow': fit2_file = root + '_wabs_' + plasma_model + '_pow.dat'
            if second_comp == 'mekal' or second_comp == 'apec': fit2_file = root + '_wabs_2' + plasma_model + '.dat'
            if second_comp == 'mkcflow': fit2_file = root + '_wabs_' + plasma_model + '_mkcflow.dat'
            fit2_file = os.path.join(root+"_spectra", fit2_file)
            if os.path.isfile(fit2_file) == False:
                sys.exit('ERROR: you must perform spectral fitting before mapping.')

            # Compare chi2 and make list which defines which fit
            # to use in mapping
            best_fit, Fprob = compare_fits(root+'_spectra/'+fit_file, fit2_file, second_comp)
        else:
            fit2_file = None
            best_fit = None
            Fprob = None
        print('\nPainting the maps...')
        paint_map(binmap, root+'_spectra/'+fit_file, vars_to_map, root=root, fit2_file=fit2_file, second_comp=second_comp, best_fit=best_fit, Fprob=Fprob, clobber=clobber)
        print('...done. \n\nOutput maps are named:')
        if second_comp == None:
            for var in vars_to_map:
                print('  '+root+'_'+var+'_map.fits')
        else:
            for var in vars_to_map:
                if var == 'nH' or var == 'chi2' or var == 'plindx' or var == 'mdot':
                    print('  '+root+'_'+var+'_map.fits')
                else:
                    print('  '+root+'_'+var+'1_map.fits')
                    print('  '+root+'_'+var+'2_map.fits')
            print('  '+root+'_Ftest_map.fits')

    else:
        parser.print_help()
