"""
Fits a model to Chandra spectra.

This module can use the output of extract_spectra.py to fit
a model to the spectra. If run as a script, multiple spectra
may be specfied by giving a text file that lists the
spectra as an argument with '@' prepended.

Version 0.6: 3/5/2011 - Updated for Ciao 4.3
Version 0.7: 1/11/2013 - Updated for Ciao 4.5
Version 0.8: 29/5/2014 - Updated for Ciao 4.6

"""
import os
import sys
import subprocess
import numpy
import multiprocessing
import threading
# Check that CIAO was initialized
if os.environ.get("CALDB") == None:
    sys.exit('Please initalize CIAO before running this script.')
from sherpa.astro.ui import *
from pychips import *
import gc


def call_sherpa_1T(spectra, redshift, nH_Gal, kT_guess, Ab_guess, root, lo_energy='0.5', hi_energy='7.0', plasma_model='mekal', min_counts=100, binning=None, reg_num_to_start=0, fix_nH_Gal=True, fix_abund=False, find_errors=False, make_plots=False, min_cnt_rate_ratio=0.3, clobber=False):
    """
    Calls Sherpa to fit a single-temperature model to one or more spectra.

    Inputs:  spectra - list of input PI files. Can be a list of
                       lists if there is more than one observation:
                       e.g., [spectra_obs1, spectra_obs2], where
                       spectra_obs1 = ['reg1.pi, 'reg2.pi', ...]
             redshift - redshift of source
             nH_Gal - Galactic N_H (10^22 cm^-2)
             kT_guess - initial guess for temperature (keV)
             Ab_guess - initial guess for abundance
             root - root of output file with fit results
             lo_energy - lower bound of fitted energy range (keV)
             hi_energy - upper bound of fitted energy range (keV)
             plasma_model - specifies the plasma model (mekal or apec)
             min_counts - minimum number of total counts required for
                          fitting
             binning - number of counts per bin for fitting
             reg_num_to_start - number to start from when numbering the
                                fit results by region
             fix_nH_Gal - if True, freezes nH_Gal
             fix_abund - if True, freezes abundance
             find_errors - if True, calculates errors
             make_plots - if True, make a plot of fit for each region
             min_cnt_rate_ratio - min ratio (relative to max cnt rate
                                  in region) below which to discard
                                  observations
             clobber - if True, overwrite existing file

    Outputs: The fits results are saved to the file:
                 root+'_wabs_'+plasma_model+'.dat'

    """
    if isinstance(spectra, str): spectra = [spectra]
    if isinstance(spectra[0], str):
        nreg = 1 # number of regions
    else:
        nreg = len(spectra[0]) # number of regions
    nobs = len(spectra) # number of observations
    fit_results_file = root + '_wabs_' + plasma_model + '.dat'

    ccolor = ["blue", "gold", "cyan", "forest", "darkred", "red", "gray", "green", "magenta", "orange", "black", "yellow", "turquoise", "firebrick", "brown", "azure", "honeydew", "lime", "mistyrose", "navy"]
    if type(redshift) == str: redshift = float(redshift)
    if type(nH_Gal) == str: nH_Gal = float(nH_Gal)
    if type(kT_guess) == str: kT_guess = float(kT_guess)
    if type(Ab_guess) == str: Ab_guess = float(Ab_guess)
    if type(lo_energy) == str: lo_energy = float(lo_energy)
    if type(hi_energy) == str: hi_energy = float(hi_energy)
    if os.path.isfile(fit_results_file) == False or clobber == True:
        results_file = open(fit_results_file, "w")
        results_file.write('# Fit results for wabs*'+plasma_model+' (zeros indicate that no fitting was performed)\n')
        results_file.write('# Reg_no.  kT  kT_loerr kT_hierr   Z    Z_loerr  Z_hierr  norm    norm_loerr norm_hierr nH_Gal  nH_loerr nH_hierr red_chisq total_counts num_bins\n')

        for i in range(nreg):
            print('\n')
            clean() # reset everything
            gc.collect() # collect garbage every step to avoid memory problems when fitting a large number (>10) of observations simultaneously
            nobs_current_reg = 0 # number of valid spectra for this region
            if nobs > 1:
                cnts = numpy.zeros(nobs) # array to store counts
                max_rate = numpy.zeros(nobs) # max count rate [counts/s/keV]
                src_id = 0 # index of source id
                good_src_ids = numpy.zeros(nobs, dtype=int) - 1
                for j in range(nobs):
                    pi_file = spectra[j][i]
                    pi_root = os.path.splitext(pi_file)[0]
                    if pi_root[-3:] == 'grp': # check if grouped or not
                        pi_root = pi_root[:-4]
                    bgd_file = pi_root[:-3] + 'bgd.pi'
                    rmf_file = pi_root + '.rmf'
                    arf_file = pi_root + '.arf'
                    pi_file_exists = os.path.isfile(pi_file)
                    bgd_file_exists = os.path.isfile(bgd_file)
                    rmf_file_exists = os.path.isfile(rmf_file)
                    arf_file_exists = os.path.isfile(arf_file)

                    if pi_file_exists and bgd_file_exists and rmf_file_exists and arf_file_exists: # make sure all required files exist before trying to load data
                        nobs_current_reg += 1
                        load_pha(src_id, pi_file)
                        if binning != None:
                            print('Grouping to ' + str(binning) + ' counts...')
                            group_counts(src_id, binning)
                        ignore_id(src_id, 0.0, lo_energy)
                        ignore_id(src_id, hi_energy, None)
                        cnts[j] = calc_data_sum(lo_energy, hi_energy, src_id) # get counts in filtered dataset
                        print('Counts for obs '+str(j+1)+': '+str(int(cnts[j])))
                        cnt_rate = get_rate(src_id, filter=True)
                        if len(cnt_rate) == 0: # when few counts (<50), get_rate can return zero-length array
                            max_rate[j] = 0.0
                        else:
                            max_rate[j] = numpy.max(cnt_rate)
                        subtract(src_id) # subtract the background
                        if src_id == 0:
                            if plasma_model == 'mekal':
                                set_source(src_id, xswabs.abs1 * xsmekal.plsm1)
                            if plasma_model == 'apec':
                                set_source(src_id, xswabs.abs1 * xsapec.plsm1)
                        else:
                            set_source(src_id, abs1 * plsm1)
                        good_src_ids[j] = src_id
                        src_id += 1
                # Filter out ignored observations
                good_src_ids_indx = numpy.where(good_src_ids >= 0)
                good_src_ids = good_src_ids[good_src_ids_indx]
                max_rate = max_rate[good_src_ids_indx]
                cnts = cnts[good_src_ids_indx]

                # If min_cnt_rate_ratio is specified, check whether the count rate
                # of any observation falls below the limit.
                if min_cnt_rate_ratio != None:
                    max_rate_overall = numpy.max(max_rate)
                    max_rate_ratios = max_rate / max_rate_overall
                    lowcr_src_ids_indx = numpy.where(max_rate_ratios < min_cnt_rate_ratio)
                    highcr_src_ids_indx = numpy.where(max_rate_ratios >= min_cnt_rate_ratio)
                    if len(lowcr_src_ids_indx) > 0:
                        lowcr_src_ids = good_src_ids[lowcr_src_ids_indx]
                        good_src_ids = good_src_ids[highcr_src_ids_indx]
                        cnts = cnts[highcr_src_ids_indx]
                        for b in range(len(lowcr_src_ids)):
                            print('Removing observation '+str(lowcr_src_ids[b]+1)+' (dataset '+str(lowcr_src_ids[b])+') for low count rate.')
                            delete_data(lowcr_src_ids[b])
                            nobs_current_reg -= 1

            if nobs == 1:
                pi_file = spectra[0][i]
                pi_root = os.path.splitext(pi_file)[0]
                if pi_root[-3:] == 'grp': # check if grouped or not
                    pi_root = pi_root[:-4]
                bgd_file = pi_root[:-3] + 'bgd.pi'
                rmf_file = pi_root + '.rmf'
                arf_file = pi_root + '.arf'
                pi_file_exists = os.path.isfile(pi_file)
                bgd_file_exists = os.path.isfile(bgd_file)
                rmf_file_exists = os.path.isfile(rmf_file)
                arf_file_exists = os.path.isfile(arf_file)

                if pi_file_exists and bgd_file_exists and rmf_file_exists and arf_file_exists: # make sure all required files exist before trying to load data
                    nobs_current_reg += 1
                    valid_obs_nums.append(1)
                    load_pha(pi_file)
                    if binning != None:
                        group_counts(src_id, binning)
                    if plasma_model == 'mekal':
                        set_source(xswabs.abs1 * xsmekal.plsm1)
                    if plasma_model == 'apec':
                        set_source(xswabs.abs1 * xsapec.plsm1)
                    ignore(0.0, lo_energy)
                    ignore(hi_energy, None)
                    cnts[0] = calc_data_sum(lo_energy, hi_energy) # get counts in filtered dataset
                    subtract()

            # Check whether total counts >= min_counts.
            # If so, fit; if not, skip the fit
            totcnts = numpy.sum(cnts)
            if totcnts >= min_counts:
                if nobs_current_reg > 1:
                    print('\nFitting '+str(nobs_current_reg)+' spectra in region '+str(i + reg_num_to_start)+' ('+str(int(totcnts))+' counts total)...')
                else:
                    print('\nFitting 1 spectrum in region '+str(i + reg_num_to_start)+' ('+str(int(totcnts))+' counts total)...')
                abs1.nH = nH_Gal
                abs1.cache = 0
                if fix_nH_Gal:
                    freeze(abs1.nH)
                else:
                    thaw(abs1.nH)
                plsm1.kt = kT_guess
                thaw(plsm1.kt)
                plsm1.abundanc = Ab_guess
                if fix_abund:
                    freeze(plsm1.abundanc)
                else:
                    thaw(plsm1.abundanc)
                plsm1.redshift = redshift
                freeze(plsm1.redshift)
                plsm1.cache = 0

                fit()
                fit_result = get_fit_results()
                red_chi2 = fit_result.rstat
                num_bins = fit_result.numpoints
                if fix_nH_Gal:
                    nH = nH_Gal
                    kT = fit_result.parvals[0]
                    if fix_abund:
                        Z = Ab_guess
                        norm = fit_result.parvals[1]
                    else:
                        Z = fit_result.parvals[1]
                        norm = fit_result.parvals[2]
                else:
                    nH = fit_result.parvals[0]
                    kT = fit_result.parvals[1]
                    if fix_abund:
                        Z = Ab_guess
                        norm = fit_result.parvals[2]
                    else:
                        Z = fit_result.parvals[2]
                        norm = fit_result.parvals[3]
                del fit_result

                if make_plots:
                    if len(good_src_ids) > 10:
                        nplots = numpy.ceil(len(good_src_ids) / 10.0)
                    else:
                        nplots = 1
                    for plot_num in range(nplots):
                        start_indx = 0 + numpy.floor(len(good_src_ids) / nplots) * plot_num
                        if plot_num == nplots - 1:
                            end_indx = len(good_src_ids)
                        else:
                            end_indx = numpy.floor(len(good_src_ids) / nplots) * (plot_num + 1.0)

                        clear() # delete any plot windows
                        add_window()
                        #add_window(["display", False])
                        set_preference("frame.transparency", "true")
                        #set_preference("window.display", "false")
                        label_ypos = 0.2
                        ccolor_indx = 0
                        for p in good_src_ids[start_indx:end_indx]:
                            plot_fit_delchi(p, overplot=True, clearwindow=False)
                            log_scale(X_AXIS)
                            set_current_plot("plot1")
                            limits(X_AXIS, 0.1, 10) # keV
                            limits(Y_AXIS, 1E-7, 0.2)
                            log_scale(Y_AXIS)
                            set_curve("crv1", ["symbol.color", ccolor[ccolor_indx], "err.color", ccolor[ccolor_indx]])
                            set_curve("crv2", ["line.color", ccolor[ccolor_indx], "err.color", ccolor[ccolor_indx]])
                            label_ypos = label_ypos / 2.0
                            add_label(4.0, label_ypos, "obs"+str(p+1))
                            set_label(["color", ccolor[ccolor_indx]])
                            set_current_plot("plot2")
                            limits(X_AXIS, 0.1, 10) # keV
                            limits(Y_AXIS, -5.0, 5.0)
                            set_curve(["symbol.color", ccolor[ccolor_indx], "err.color", ccolor[ccolor_indx]])
                            set_plot_xlabel("Energy (keV)")
                            set_plot_ylabel("Sigma")
                            ccolor_indx += 1
                        set_current_plot("plot1")
                        set_plot_title("Fit for region "+str(i+reg_num_to_start))
                        set_plot_ylabel("Counts s^{-1} keV^{-1}")
                        if nplots > 1:
                            fit_plot_file = "reg" + str(i+reg_num_to_start) + "_plot" + str(plot_num) + "_1T_fit.pdf"
                        else:
                            fit_plot_file = "reg" + str(i+reg_num_to_start) + "_1T_fit.pdf"
                        print_window(fit_plot_file, ["orientation", "portrait", "clobber", True])

                if find_errors:
                    covar()
                    covar_result = get_covar_results()
                    if fix_nH_Gal:
                        nH_loerr = 0.0
                        nH_hierr = 0.0
                        kT_loerr = covar_result.parmins[0]
                        kT_hierr = covar_result.parmaxes[0]
                        if fix_abund:
                            Z_loerr = 0.0
                            Z_hierr = 0.0
                            norm_loerr = covar_result.parmins[1]
                            norm_hierr = covar_result.parmaxes[1]
                        else:
                            Z_loerr = covar_result.parmins[1]
                            Z_hierr = covar_result.parmaxes[1]
                            norm_loerr = covar_result.parmins[2]
                            norm_hierr = covar_result.parmaxes[2]
                    else:
                        nH_loerr =covar_result.parmins[0]
                        kT_loerr = covar_result.parmins[1]
                        nH_hierr = covar_result.parmaxes[0]
                        kT_hierr = covar_result.parmaxes[1]
                        if fix_abund:
                            Z_loerr = 0.0
                            Z_hierr = 0.0
                            norm_loerr = covar_result.parmins[2]
                            norm_hierr = covar_result.parmaxes[2]
                        else:
                            Z_loerr = covar_result.parmins[2]
                            Z_hierr = covar_result.parmaxes[2]
                            norm_loerr = covar_result.parmins[3]
                            norm_hierr = covar_result.parmaxes[3]
                    del covar_result

                    # Check for failed errors (= None) and set them to +/- best-fit value
                    if not fix_nH_Gal:
                        if nH_loerr == None: nH_loerr = -nH
                        if nH_hierr == None: nH_hierr = nH
                    if kT_loerr == None: kT_loerr = -kT
                    if kT_hierr == None: kT_hierr = kT
                    if not fix_abund:
                        if Z_loerr == None: Z_loerr = -Z
                        if Z_hierr == None: Z_hierr = Z
                    if norm_loerr == None: norm_loerr = -norm
                    if norm_hierr == None: norm_hierr = norm
                else:
                    kT_loerr = 0.0
                    Z_loerr = 0.0
                    nH_loerr = 0.0
                    norm_loerr = 0.0
                    kT_hierr = 0.0
                    Z_hierr = 0.0
                    nH_hierr = 0.0
                    norm_hierr = 0.0

            else: # if total counts < min_counts, just write zeros
                print('\n  Warning: no fit performed for for region '+str(i + reg_num_to_start)+':')
                if nobs > 1:
                    print('  Spectra have insufficient counts after filtering or do not exist.')
                else:
                    print('  Spectrum has insufficient counts after filtering or does not exist.')
                print('  --> All parameters for this region set to 0.0.')
                kT = 0.0
                Z = 0.0
                nH = 0.0
                norm = 0.0
                kT_loerr = 0.0
                Z_loerr = 0.0
                nH_loerr = 0.0
                norm_loerr = 0.0
                kT_hierr = 0.0
                Z_hierr = 0.0
                nH_hierr = 0.0
                norm_hierr = 0.0
                red_chi2 = 0.0
                num_bins = 0

            results_file.write('%7r %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.4e %6.4e %6.4e %7.4f %7.4f %7.4f %7.4f %8.1f %8r\n' % (i+reg_num_to_start, kT, kT_loerr, kT_hierr, Z, Z_loerr, Z_hierr, norm, norm_loerr, norm_hierr, nH, nH_loerr, nH_hierr, red_chi2, totcnts, num_bins) )

        results_file.close()

        # Finally, make sure all regions have an entry in the fit results file
        # (this shouldn't be needed, but just in case...)
        dtype = {'names': ('reg_id', 'kT', 'kT_lo', 'kT_hi', 'Z', 'Z_lo', 'Z_hi', 'norm', 'norm_lo', 'norm_hi', 'nH', 'nH_lo', 'nH_hi', 'chi2', 'totcnts', 'nbins'), 'formats': ('i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4')}
        data = numpy.loadtxt(fit_results_file, dtype=dtype)
        reg_num = data["reg_id"]
        missing_regs = []
        n_missing = 0
        for i in range(len(reg_num)):
            if int(reg_num[i]) != i + reg_num_to_start + n_missing:
                missing_regs.append(i + reg_num_to_start + n_missing)
                n_missing += 1
        if n_missing > 0:
            results_file = open(fit_results_file, "a") # append missing regions
            for i in range(n_missing):
                results_file.write('%7r %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.4e %6.4e %6.4e %7.4f %7.4f %7.4f %7.4f %8.1f %8r\n' % (missing_regs[i]+reg_num_to_start, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0) )

    else:
        print('\n  Output file ('+fit_results_file+') exists and clobber = False.')


def call_sherpa_2T(spectra, redshift, nH_Gal, kT_guess, Ab_guess, root, lo_energy='0.5', hi_energy='7.0', plasma_model='mekal', min_counts=100, binning=None, reg_num_to_start=0, fix_nH_Gal=True, fix_abund=False, find_errors=False, make_plots=False, min_cnt_rate_ratio=0.3, clobber=False):
    """
    Calls Sherpa to fit a two-temperature model to one or more spectra.

    Inputs:  spectra - list of input PI files. Can be a list of
                       lists if there is more than one observation:
                       e.g., [spectra_obs1, spectra_obs2], where
                       spectra_obs1 = ['reg1.pi, 'reg2.pi', ...]
             redshift - redshift of source
             nH_Gal - Galactic N_H (10^22 cm^-2)
             kT_guess - initial guess for temperature (keV)
             Ab_guess - initial guess for abundance
             root - root of output file with fit results
             lo_energy - lower bound of fitted energy range (keV)
             hi_energy - upper bound of fitted energy range (keV)
             plasma_model - specifies the plasma model (mekal or apec)
             min_counts - minimum number of total counts required for
                          fitting
             binning - number of counts per bin for fitting
             reg_num_to_start - number to start from when numbering the
                                fit results by region
             fix_nH_Gal - if True, freezes nH_Gal
             fix_abund - if True, freezes abundance
             find_errors - if True, calculates errors
             make_plots - if True, make a plot of fit for each region
             min_cnt_rate_ratio - min ratio (relative to max cnt rate
                                  in region) below which to discard
                                  observations
             clobber - if True, overwrite existing file

    Outputs: The fits results are saved to the file:
                 root+'_wabs_2'+plasma_model+'.dat'

    """
    if isinstance(spectra, str): spectra = [spectra]
    if isinstance(spectra[0], str):
        nreg = 1 # number of regions
    else:
        nreg = len(spectra[0]) # number of regions
    nobs = len(spectra) # number of observations
    fit_results_file = root + '_wabs_2' + plasma_model + '.dat'

    ccolor = ["blue", "gold", "cyan", "forest", "darkred", "red", "gray", "green", "magenta", "orange", "black", "yellow", "turquoise", "firebrick", "brown", "azure", "honeydew", "lime", "mistyrose", "navy"]
    if type(redshift) == str: redshift = float(redshift)
    if type(nH_Gal) == str: nH_Gal = float(nH_Gal)
    if type(kT_guess) == str: kT_guess = float(kT_guess)
    if type(Ab_guess) == str: Ab_guess = float(Ab_guess)
    if type(lo_energy) == str: lo_energy = float(lo_energy)
    if type(hi_energy) == str: hi_energy = float(hi_energy)
    if os.path.isfile(fit_results_file) == False or clobber == True:
        results_file = open(fit_results_file, "w")
        results_file.write('# Fit results for wabs*'+plasma_model+' (zeros indicate that no fitting was performed)\n')
        results_file.write('# Reg_no.  kT1 kT1_loerr kT1_hierr Z1   Z1_loerr Z1_hierr norm1  norm1_loerr norm1_hierr kT2 kT2_loerr kT2_hierr  Z2   Z2_loerr Z2_hierr  norm2 norm2_loerr norm2_hierr nH_Gal  nH_loerr nH_hierr red_chisq total_counts num_bins\n')
        if make_plots:
            add_window(["display", False])

        for i in range(nreg):
            print('\n')
            clean() # reset everything
            gc.collect() # collect garbage every step to avoid memory problems when fitting a large number (>10) of observations simultaneously
            nobs_current_reg = 0 # number of valid spectra for this region
            if nobs > 1:
                cnts = numpy.zeros(nobs) # array to store counts
                max_rate = numpy.zeros(nobs) # max count rate [counts/s/keV]
                src_id = 0 # index of source id
                good_src_ids = numpy.zeros(nobs, dtype=int) - 1
                for j in range(nobs):
                    pi_file = spectra[j][i]
                    pi_root = os.path.splitext(pi_file)[0]
                    if pi_root[-3:] == 'grp': # check if grouped or not
                        pi_root = pi_root[:-4]
                    bgd_file = pi_root[:-3] + 'bgd.pi'
                    rmf_file = pi_root + '.rmf'
                    arf_file = pi_root + '.arf'
                    pi_file_exists = os.path.isfile(pi_file)
                    bgd_file_exists = os.path.isfile(bgd_file)
                    rmf_file_exists = os.path.isfile(rmf_file)
                    arf_file_exists = os.path.isfile(arf_file)

                    if pi_file_exists and bgd_file_exists and rmf_file_exists and arf_file_exists: # make sure all required files exist before trying to load data
                        nobs_current_reg += 1
                        load_pha(src_id, pi_file)
                        if binning != None:
                            print('Grouping to ' + str(binning) + ' counts...')
                            group_counts(src_id, binning)
                        ignore_id(src_id, 0.0, lo_energy)
                        ignore_id(src_id, hi_energy, None)
                        cnts[j] = calc_data_sum(lo_energy, hi_energy, src_id) # get counts in filtered dataset
                        print('Counts for obs '+str(j+1)+': '+str(int(cnts[j])))
                        cnt_rate = get_rate(src_id, filter=True)
                        if len(cnt_rate) == 0: # when few counts (<50), get_rate can return zero-length array
                            max_rate[j] = 0.0
                        else:
                            max_rate[j] = numpy.max(cnt_rate)
                        subtract(src_id) # subtract the background
                        if src_id == 0:
                            if plasma_model == 'mekal':
                                set_source(src_id, xswabs.abs1 * (xsmekal.plsm1 + xsmekal.plsm2))
                            if plasma_model == 'apec':
                                set_source(src_id, xswabs.abs1 * (xsapec.plsm1 + xsapec.plsm2))
                        else:
                            set_source(src_id, abs1 * (plsm1 + plsm2))
                        good_src_ids[j] = src_id
                        src_id += 1
                # Filter out ignored observations
                good_src_ids_indx = numpy.where(good_src_ids >= 0)
                good_src_ids = good_src_ids[good_src_ids_indx]
                max_rate = max_rate[good_src_ids_indx]
                cnts = cnts[good_src_ids_indx]

                # If min_cnt_rate_ratio is specified, check whether the count rate
                # of any observation falls below the limit.
                if min_cnt_rate_ratio != None:
                    max_rate_overall = numpy.max(max_rate)
                    max_rate_ratios = max_rate / max_rate_overall
                    lowcr_src_ids_indx = numpy.where(max_rate_ratios < min_cnt_rate_ratio)
                    highcr_src_ids_indx = numpy.where(max_rate_ratios >= min_cnt_rate_ratio)
                    if len(lowcr_src_ids_indx) > 0:
                        lowcr_src_ids = good_src_ids[lowcr_src_ids_indx]
                        good_src_ids = good_src_ids[highcr_src_ids_indx]
                        cnts = cnts[highcr_src_ids_indx]
                        for b in range(len(lowcr_src_ids)):
                            print('Removing observation '+str(lowcr_src_ids[b]+1)+' (dataset '+str(lowcr_src_ids[b])+') for low count rate.')
                            delete_data(lowcr_src_ids[b])
                            nobs_current_reg -= 1

            if nobs == 1:
                pi_file = spectra[i]
                pi_root = os.path.splitext(pi_file)[0]
                if pi_root[-3:] == 'grp': # check if grouped or not
                    pi_root = pi_root[:-4]
                bgd_file = pi_root[:-3] + 'bgd.pi'
                rmf_file = pi_root + '.rmf'
                arf_file = pi_root + '.arf'
                pi_file_exists = os.path.isfile(pi_file)
                bgd_file_exists = os.path.isfile(bgd_file)
                rmf_file_exists = os.path.isfile(rmf_file)
                arf_file_exists = os.path.isfile(arf_file)

                if pi_file_exists and bgd_file_exists and rmf_file_exists and arf_file_exists: # make sure all required files exist before trying to load data
                    nobs_current_reg += 1
                    load_pha(pi_file)
                    if binning != None:
                        group_counts(src_id, binning)
                    if plasma_model == 'mekal':
                        set_source(xswabs.abs1 * (xsmekal.plsm1 + xsmekal.plsm2))
                    if plasma_model == 'apec':
                        set_source(xswabs.abs1 * (xsapec.plsm1 + xsapec.plsm2))
                    ignore(0.0, lo_energy)
                    ignore(hi_energy, None)
                    cnts[0] = calc_data_sum(lo_energy, hi_energy) # get counts in filtered dataset
                    subtract()

            # Check whether total counts >= min_counts.
            # If so, fit; if not, skip the fit
            totcnts = numpy.sum(cnts)
            if totcnts >= min_counts:
                if nobs_current_reg > 1:
                    print('\nFitting '+str(nobs_current_reg)+' spectra in region '+str(i + reg_num_to_start)+' ('+str(totcnts)+' counts total)...')
                else:
                    print('\nFitting 1 spectrum in region '+str(i + reg_num_to_start)+' ('+str(numpy.sum(cnts))+' counts total)...')
                abs1.nH = nH_Gal
                abs1.cache = 0
                if fix_nH_Gal:
                    freeze(abs1.nH)
                else:
                    thaw(abs1.nH)
                plsm1.kt = kT_guess
                thaw(plsm1.kt)
                plsm1.abundanc = Ab_guess
                if fix_abund:
                    freeze(plsm1.abundanc)
                else:
                    thaw(plsm1.abundanc)
                plsm1.redshift = redshift
                freeze(plsm1.redshift)
                plsm1.cache = 0
                plsm2.kt = kT_guess
                thaw(plsm2.kt)
                plsm2.abundanc = Ab_guess
                if fix_abund:
                    freeze(plsm2.abundanc)
                else:
                    thaw(plsm2.abundanc)
                link(plsm1.abundanc, plsm2.abundanc)
                plsm2.redshift = redshift
                freeze(plsm2.redshift)
                plsm2.cache = 0

                set_method("moncar")
                fit()
                fit_result = get_fit_results()
                red_chi2 = fit_result.rstat
                num_bins = fit_result.numpoints
                if fix_nH_Gal:
                    nH = nH_Gal
                    kT1 = fit_result.parvals[0]
                    if fix_abund:
                        Z1 = Ab_guess
                        norm1 = fit_result.parvals[1]
                        kT2 = fit_result.parvals[2]
                        norm2 = fit_result.parvals[3]
                    else:
                        Z1 = fit_result.parvals[1]
                        norm1 = fit_result.parvals[2]
                        kT2 = fit_result.parvals[3]
                        norm2 = fit_result.parvals[4]
                else:
                    nH = fit_result.parvals[0]
                    kT1 = fit_result.parvals[1]
                    if fix_abund:
                        Z1 = Ab_guess
                        norm1 = fit_result.parvals[2]
                        kT2 = fit_result.parvals[3]
                        norm2 = fit_result.parvals[4]
                    else:
                        Z1 = fit_result.parvals[2]
                        norm1 = fit_result.parvals[3]
                        kT2 = fit_result.parvals[4]
                        norm2 = fit_result.parvals[5]
                del fit_result

                if find_errors:
                    covar()
                    covar_result = get_covar_results()
                    if fix_nH_Gal:
                        nH_loerr = 0.0
                        nH_hierr = 0.0
                        kT1_loerr = covar_result.parmins[0]
                        kT1_hierr = covar_result.parmaxes[0]
                        if fix_abund:
                            Z1_loerr = 0.0
                            Z1_hierr = 0.0
                            norm1_loerr = covar_result.parmins[1]
                            norm1_hierr = covar_result.parmaxes[1]
                            kT2_loerr = covar_result.parmins[2]
                            norm2_loerr = covar_result.parmins[3]
                            kT2_hierr = covar_result.parmaxes[2]
                            norm2_hierr = covar_result.parmaxes[3]
                        else:
                            Z1_loerr = covar_result.parmins[1]
                            norm1_loerr = covar_result.parmins[2]
                            kT1_hierr = covar_result.parmaxes[0]
                            Z1_hierr = covar_result.parmaxes[1]
                            norm1_hierr = covar_result.parmaxes[2]
                            kT2_loerr = covar_result.parmins[3]
                            norm2_loerr = covar_result.parmins[4]
                            kT2_hierr = covar_result.parmaxes[3]
                            norm2_hierr = covar_result.parmaxes[4]
                    else:
                        nH_loerr =covar_result.parmins[0]
                        nH_hierr = covar_result.parmaxes[0]
                        kT1_loerr = covar_result.parmins[1]
                        kT1_hierr = covar_result.parmaxes[1]
                        if fix_abund:
                            Z1_loerr = 0.0
                            Z1_hierr = 0.0
                            norm1_loerr = covar_result.parmins[2]
                            norm1_hierr = covar_result.parmaxes[2]
                            kT2_loerr = covar_result.parmins[3]
                            norm2_loerr = covar_result.parmins[4]
                            kT2_hierr = covar_result.parmaxes[3]
                            norm2_hierr = covar_result.parmaxes[4]
                        else:
                            Z1_loerr = covar_result.parmins[2]
                            norm1_loerr = covar_result.parmins[3]
                            Z1_hierr = covar_result.parmaxes[2]
                            norm1_hierr = covar_result.parmaxes[3]
                            kT2_loerr = covar_result.parmins[4]
                            norm2_loerr = covar_result.parmins[5]
                            kT2_hierr = covar_result.parmaxes[4]
                            norm2_hierr = covar_result.parmaxes[5]
                    del covar_result

                    # Check for failed errors (= None) and set them to +/- best-fit value
                    if not fix_nH_Gal:
                        if nH_loerr == None: nH_loerr = -nH
                        if nH_hierr == None: nH_hierr = nH
                    if kT1_loerr == None: kT1_loerr = -kT1
                    if kT1_hierr == None: kT1_hierr = kT1
                    if not fix_abund:
                        if Z1_loerr == None: Z1_loerr = -Z1
                        if Z1_hierr == None: Z1_hierr = Z1
                    if norm1_loerr == None: norm1_loerr = -norm1
                    if norm1_hierr == None: norm1_hierr = norm1
                    if kT2_loerr == None: kT2_loerr = -kT2
                    if kT2_hierr == None: kT2_hierr = kT2
                    if norm2_loerr == None: norm2_loerr = -norm2
                    if norm2_hierr == None: norm2_hierr = norm2
                else:
                    kT1_loerr = 0.0
                    Z1_loerr = 0.0
                    nH_loerr = 0.0
                    norm1_loerr = 0.0
                    kT1_hierr = 0.0
                    Z1_hierr = 0.0
                    nH_hierr = 0.0
                    norm1_hierr = 0.0
                    kT2_loerr = 0.0
                    norm2_loerr = 0.0
                    kT2_hierr = 0.0
                    norm2_hierr = 0.0

            else: # if total counts < min_counts, just write zeros
                print('\n  Warning: no fit performed for for region '+str(i + reg_num_to_start)+':')
                if nobs > 1:
                    print('  Spectra have insufficient counts after filtering or do not exist.')
                else:
                    print('  Spectrum has insufficient counts after filtering or does not exist.')
                print('  --> All parameters for this region set to 0.0.')
                kT1 = 0.0
                Z1 = 0.0
                nH = 0.0
                norm1 = 0.0
                kT2 = 0.0
                norm2 = 0.0
                kT1_loerr = 0.0
                Z1_loerr = 0.0
                nH_loerr = 0.0
                norm1_loerr = 0.0
                kT1_hierr = 0.0
                Z1_hierr = 0.0
                nH_hierr = 0.0
                norm1_hierr = 0.0
                kT2_loerr = 0.0
                norm2_loerr = 0.0
                kT2_hierr = 0.0
                norm2_hierr = 0.0
                red_chi2 = 0.0
                num_bins = 0

            results_file.write('%7r %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.4e %6.4e %6.4e %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.4e %6.4e %6.4e %7.4f %7.4f %7.4f %7.4f %8.1f %8r\n' % (i+reg_num_to_start, kT1, kT1_loerr, kT1_hierr, Z1, Z1_loerr, Z1_hierr, norm1, norm1_loerr, norm1_hierr, kT2, kT2_loerr, kT2_hierr, Z1, Z1_loerr, Z1_hierr, norm2, norm2_loerr, norm2_hierr, nH, nH_loerr, nH_hierr, red_chi2, totcnts, num_bins) )

        results_file.close()
    else:
        print('\n  Output file ('+fit_results_file+') exists and clobber = False.')


def call_sherpa_1T_plus_pow(spectra, redshift, nH_Gal, kT_guess, Ab_guess, plindx_guess, root, lo_energy='0.5', hi_energy='7.0', plasma_model='mekal', min_counts=100, binning=None, reg_num_to_start=0, fix_nH_Gal=True, fix_abund=False, find_errors=False, min_cnt_rate_ratio=0.3, make_plots=False, clobber=False):
    """
    Calls Sherpa to fit a single-temperature-plus-power-law model to one or more spectra.

    Inputs:  spectra - list of input PI files. Can be a list of
                       lists if there is more than one observation:
                       e.g., [spectra_obs1, spectra_obs2], where
                       spectra_obs1 = ['reg1.pi, 'reg2.pi', ...]
             redshift - redshift of source
             nH_Gal - Galactic N_H (10^22 cm^-2)
             kT_guess - initial guess for temperature (keV)
             Ab_guess - initial guess for abundance
             plindx_guess - intial guess for the power-law index
             root - root of output file with fit results
             lo_energy - lower bound of fitted energy range (keV)
             hi_energy - upper bound of fitted energy range (keV)
             plasma_model - specifies the plasma model (mekal or apec)
             min_counts - minimum number of total counts required for
                          fitting
             binning - number of counts per bin for fitting
             reg_num_to_start - number to start from when numbering the
                                fit results by region
             fix_nH_Gal - if True, freezes nH_Gal
             fix_abund - if True, freezes abundance
             find_errors - if True, calculates errors
             make_plots - if True, make a plot of fit for each region
             min_cnt_rate_ratio - min ratio (relative to max cnt rate
                                  in region) below which to discard
                                  observations
             clobber - if True, overwrite existing file

    Outputs: The fits results are saved to the file:
                 root+'_wabs_'+plasma_model+'_pow.dat'

    """
    if isinstance(spectra, str): spectra = [spectra]
    if isinstance(spectra[0], str):
        nreg = 1 # number of regions
    else:
        nreg = len(spectra[0]) # number of regions
    nobs = len(spectra) # number of observations
    fit_results_file = root + '_wabs_' + plasma_model + '_pow.dat'

    ccolor = ["blue", "gold", "cyan", "forest", "darkred", "red", "gray", "green", "magenta", "orange", "black", "yellow", "turquoise", "firebrick", "brown", "azure", "honeydew", "lime", "mistyrose", "navy"]
    if type(redshift) == str: redshift = float(redshift)
    if type(nH_Gal) == str: nH_Gal = float(nH_Gal)
    if type(kT_guess) == str: kT_guess = float(kT_guess)
    if type(Ab_guess) == str: Ab_guess = float(Ab_guess)
    if type(plindx_guess) == str: plindx_guess = float(plindx_guess)
    if type(lo_energy) == str: lo_energy = float(lo_energy)
    if type(hi_energy) == str: hi_energy = float(hi_energy)
    if os.path.isfile(fit_results_file) == False or clobber == True:
        results_file = open(fit_results_file, "w")
        results_file.write('# Fit results for wabs*'+plasma_model+' (zeros indicate that no fitting was performed)\n')
        results_file.write('# Reg_no.  kT  kT_loerr kT_hierr   Z    Z_loerr  Z_hierr  norm   norm_loerr norm_hierr plindx  indx_loerr indx_hierr plnorm plnorm_loerr plnorm_hierr nH_Gal  nH_loerr nH_hierr red_chisq total_counts num_bins\n')

        for i in range(nreg):
            print('\n')
            clean() # reset everything
            gc.collect() # collect garbage every step to avoid memory problems when fitting a large number (>10) of observations simultaneously
            nobs_current_reg = 0 # number of valid spectra for this region
            good_src_ids = numpy.zeros(nobs, dtype=int) - 1
            if nobs > 1:
                cnts = numpy.zeros(nobs) # array to store counts
                max_rate = numpy.zeros(nobs) # max count rate [counts/s/keV]
                src_id = 0 # index of source id
                for j in range(nobs):
                    pi_file = spectra[j][i]
                    pi_root = os.path.splitext(pi_file)[0]
                    if pi_root[-3:] == 'grp': # check if grouped or not
                        pi_root = pi_root[:-4]
                    bgd_file = pi_root[:-3] + 'bgd.pi'
                    rmf_file = pi_root + '.rmf'
                    arf_file = pi_root + '.arf'
                    pi_file_exists = os.path.isfile(pi_file)
                    bgd_file_exists = os.path.isfile(bgd_file)
                    rmf_file_exists = os.path.isfile(rmf_file)
                    arf_file_exists = os.path.isfile(arf_file)

                    if pi_file_exists and bgd_file_exists and rmf_file_exists and arf_file_exists: # make sure all required files exist before trying to load data
                        nobs_current_reg += 1
                        load_pha(src_id, pi_file)
                        if binning != None:
                            print('Grouping to ' + str(binning) + ' counts...')
                            group_counts(src_id, binning)
                        ignore_id(src_id, 0.0, lo_energy)
                        ignore_id(src_id, hi_energy, None)
                        cnts[j] = calc_data_sum(lo_energy, hi_energy, src_id) # get counts in filtered dataset
                        print('Counts for obs '+str(j+1)+': '+str(int(cnts[j])))
                        cnt_rate = get_rate(src_id, filter=True)
                        if len(cnt_rate) == 0: # when few counts (<50), get_rate can return zero-length array
                            max_rate[j] = 0.0
                        else:
                            max_rate[j] = numpy.max(cnt_rate)
                        subtract(src_id) # subtract the background
                        if src_id == 0:
                            if plasma_model == 'mekal':
                                set_source(src_id, xswabs.abs1 * (xsmekal.plsm1 + xspowerlaw.pow1))
                            if plasma_model == 'apec':
                                set_source(src_id, xswabs.abs1 * (xsapec.plsm1 + xspowerlaw.pow1))
                        else:
                            set_source(src_id, abs1 * (plsm1 + pow1))
                        good_src_ids[j] = src_id
                        src_id += 1
                # Filter out ignored observations
                good_src_ids_indx = numpy.where(good_src_ids >= 0)
                good_src_ids = good_src_ids[good_src_ids_indx]
                max_rate = max_rate[good_src_ids_indx]
                cnts = cnts[good_src_ids_indx]

                # If min_cnt_rate_ratio is specified, check whether the count rate
                # of any observation falls below the limit.
                if min_cnt_rate_ratio != None:
                    max_rate_overall = numpy.max(max_rate)
                    max_rate_ratios = max_rate / max_rate_overall
                    lowcr_src_ids_indx = numpy.where(max_rate_ratios < min_cnt_rate_ratio)
                    highcr_src_ids_indx = numpy.where(max_rate_ratios >= min_cnt_rate_ratio)
                    if len(lowcr_src_ids_indx) > 0:
                        lowcr_src_ids = good_src_ids[lowcr_src_ids_indx]
                        good_src_ids = good_src_ids[highcr_src_ids_indx]
                        cnts = cnts[highcr_src_ids_indx]
                        for b in range(len(lowcr_src_ids)):
                            print('Removing observation '+str(lowcr_src_ids[b]+1)+' (dataset '+str(lowcr_src_ids[b])+') for low count rate.')
                            delete_data(lowcr_src_ids[b])
                            nobs_current_reg -= 1

            if nobs == 1:
                pi_file = spectra[i]
                pi_root = os.path.splitext(pi_file)[0]
                if pi_root[-3:] == 'grp': # check if grouped or not
                    pi_root = pi_root[:-4]
                bgd_file = pi_root[:-3] + 'bgd.pi'
                rmf_file = pi_root + '.rmf'
                arf_file = pi_root + '.arf'
                pi_file_exists = os.path.isfile(pi_file)
                bgd_file_exists = os.path.isfile(bgd_file)
                rmf_file_exists = os.path.isfile(rmf_file)
                arf_file_exists = os.path.isfile(arf_file)

                if pi_file_exists and bgd_file_exists and rmf_file_exists and arf_file_exists: # make sure all required files exist before trying to load data
                    nobs_current_reg += 1
                    load_pha(pi_file)
                    if binning != None:
                        group_counts(src_id, binning)
                    if plasma_model == 'mekal':
                        set_source(xswabs.abs1 * (xsmekal.plsm1 + xspowerlaw.pow1))
                    if plasma_model == 'apec':
                        set_source(xswabs.abs1 * (xsapec.plsm1 + xspowerlaw.pow1))
                    ignore(0.0, lo_energy)
                    ignore(hi_energy, None)
                    cnts[0] = calc_data_sum(lo_energy, hi_energy) # get counts in filtered dataset
                    subtract()

            # Check whether total counts >= min_counts.
            # If so, fit; if not, skip the fit
            totcnts = numpy.sum(cnts)
            if totcnts >= min_counts:
                if nobs_current_reg > 1:
                    print('\nFitting '+str(nobs_current_reg)+' spectra in region '+str(i + reg_num_to_start)+' ('+str(totcnts)+' counts total)...')
                else:
                    print('\nFitting 1 spectrum in region '+str(i + reg_num_to_start)+' ('+str(numpy.sum(cnts))+' counts total)...')
                abs1.nH = nH_Gal
                if fix_nH_Gal:
                    freeze(abs1.nH)
                else:
                    thaw(abs1.nH)
                plsm1.kt = kT_guess
                thaw(plsm1.kt)
                plsm1.abundanc = Ab_guess
                thaw(plsm1.abundanc)
                plsm1.redshift = redshift
                freeze(plsm1.redshift)
                pow1.PhoIndex = plindx_guess

                set_method("moncar")
                fit()
                fit_result = get_fit_results()
                red_chi2 = fit_result.rstat
                num_bins = fit_result.numpoints
                if fix_nH_Gal:
                    nH = nH_Gal
                    kT = fit_result.parvals[0]
                    Z = fit_result.parvals[1]
                    norm = fit_result.parvals[2]
                    powindx = fit_result.parvals[3]
                    pownorm = fit_result.parvals[4]
                else:
                    nH = fit_result.parvals[0]
                    kT = fit_result.parvals[1]
                    Z = fit_result.parvals[2]
                    norm = fit_result.parvals[3]
                    powindx = fit_result.parvals[4]
                    pownorm = fit_result.parvals[5]
                del fit_result

                if make_plots:
                    if len(good_src_ids) > 10:
                        nplots = numpy.ceil(len(good_src_ids) / 10.0)
                    else:
                        nplots = 1
                    for plot_num in range(nplots):
                        start_indx = 0 + numpy.floor(len(good_src_ids) / nplots) * plot_num
                        if plot_num == nplots - 1:
                            end_indx = len(good_src_ids)
                        else:
                            end_indx = numpy.floor(len(good_src_ids) / nplots) * (plot_num + 1.0)

                        clear() # delete any plot windows
                        add_window()
                        #add_window(["display", False])
                        set_preference("frame.transparency", "true")
                        #set_preference("window.display", "false")
                        label_ypos = 0.2
                        ccolor_indx = 0
                        for p in good_src_ids[start_indx:end_indx]:
                            plot_fit_delchi(p, overplot=True, clearwindow=False)
                            log_scale(X_AXIS)
                            set_current_plot("plot1")
                            limits(X_AXIS, 0.1, 10) # keV
                            limits(Y_AXIS, 1E-7, 0.2)
                            log_scale(Y_AXIS)
                            set_curve("crv1", ["symbol.color", ccolor[ccolor_indx], "err.color", ccolor[ccolor_indx]])
                            set_curve("crv2", ["line.color", ccolor[ccolor_indx], "err.color", ccolor[ccolor_indx]])
                            label_ypos = label_ypos / 2.0
                            add_label(4.0, label_ypos, "obs"+str(p+1))
                            set_label(["color", ccolor[ccolor_indx]])
                            set_current_plot("plot2")
                            limits(X_AXIS, 0.1, 10) # keV
                            limits(Y_AXIS, -5.0, 5.0)
                            set_curve(["symbol.color", ccolor[ccolor_indx], "err.color", ccolor[ccolor_indx]])
                            set_plot_xlabel("Energy (keV)")
                            set_plot_ylabel("Sigma")
                            ccolor_indx += 1
                        set_current_plot("plot1")
                        set_plot_title("Fit for region "+str(i+reg_num_to_start))
                        set_plot_ylabel("Counts s^{-1} keV^{-1}")
                        if nplots > 1:
                            fit_plot_file = "reg" + str(i+reg_num_to_start) + "_plot" + str(plot_num) + "_1T+pow_fit.pdf"
                        else:
                            fit_plot_file = "reg" + str(i+reg_num_to_start) + "_1T+pow_fit.pdf"
                        print_window(fit_plot_file, ["orientation", "portrait", "clobber", True])

                if find_errors:
                    covar()
                    covar_result = get_covar_results()
                    if fix_nH_Gal:
                        nH_loerr = 0.0
                        nH_hierr = 0.0
                        kT_loerr = covar_result.parmins[0]
                        Z_loerr = covar_result.parmins[1]
                        norm_loerr = covar_result.parmins[2]
                        powindx_loerr = covar_result.parmins[3]
                        pownorm_loerr = covar_result.parmins[4]
                        kT_hierr = covar_result.parmaxes[0]
                        Z_hierr = covar_result.parmaxes[1]
                        norm_hierr = covar_result.parmaxes[2]
                        powindx_hierr = covar_result.parmaxes[3]
                        pownorm_hierr = covar_result.parmaxes[4]
                    else:
                        nH_loerr =covar_result.parmins[0]
                        kT_loerr = covar_result.parmins[1]
                        Z_loerr = covar_result.parmins[2]
                        norm_loerr = covar_result.parmins[3]
                        powindx_loerr = covar_result.parmins[4]
                        pownorm_loerr = covar_result.parmins[5]
                        nH_hierr = covar_result.parmaxes[0]
                        kT_hierr = covar_result.parmaxes[1]
                        Z_hierr = covar_result.parmaxes[2]
                        norm_hierr = covar_result.parmaxes[3]
                        powindx_hierr = covar_result.parmaxes[4]
                        pownorm_hierr = covar_result.parmaxes[5]
                    del covar_result

                    # Check for failed errors (= None) and set them to +/- best-fit value
                    if not fix_nH_Gal:
                        if nH_loerr == None: nH_loerr = -nH
                        if nH_hierr == None: nH_hierr = nH
                    if kT_loerr == None: kT_loerr = -kT
                    if kT_hierr == None: kT_hierr = kT
                    if Z_loerr == None: Z_loerr = -Z
                    if Z_hierr == None: Z_hierr = Z
                    if norm_loerr == None: norm_loerr = -norm
                    if norm_hierr == None: norm_hierr = norm
                    if powindx_loerr == None: powindx_loerr = -powindx
                    if powindx_hierr == None: powindx_hierr = powindx
                    if pownorm_loerr == None: pownorm_loerr = -pownorm
                    if pownorm_hierr == None: pownorm_hierr = pownorm

                else:
                    kT_loerr = 0.0
                    Z_loerr = 0.0
                    nH_loerr = 0.0
                    norm_loerr = 0.0
                    powindx_loerr = 0.0
                    pownorm_loerr = 0.0
                    kT_hierr = 0.0
                    Z_hierr = 0.0
                    nH_hierr = 0.0
                    norm_hierr = 0.0
                    powindx_hierr = 0.0
                    pownorm_hierr = 0.0

            else: # if total counts < min_counts, just write zeros
                print('\n  Warning: no fit performed for for region '+str(i + reg_num_to_start)+':')
                if nobs > 1:
                    print('  Spectra have insufficient counts after filtering or do not exist.')
                else:
                    print('  Spectrum has insufficient counts after filtering or does not exist.')
                print('  --> All parameters for this region set to 0.0.')
                kT = 0.0
                Z = 0.0
                nH = 0.0
                norm = 0.0
                kT_loerr = 0.0
                Z_loerr = 0.0
                nH_loerr = 0.0
                norm_loerr = 0.0
                powindx_loerr = 0.0
                pownorm_loerr = 0.0
                kT_hierr = 0.0
                Z_hierr = 0.0
                nH_hierr = 0.0
                norm_hierr = 0.0
                powindx_hierr = 0.0
                pownorm_hierr = 0.0
                red_chi2 = 0.0

            results_file.write('%7r %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.4e %6.4e %6.4e %7.4f %7.4f %7.4f %6.4e %6.4e %6.4e %7.4f %7.4f %7.4f %7.4f %8.1f %8r\n' % (i+reg_num_to_start, kT, kT_loerr, kT_hierr, Z, Z_loerr, Z_hierr, norm, norm_loerr, norm_hierr, powindx, powindx_loerr, powindx_hierr, pownorm, pownorm_loerr, pownorm_hierr, nH, nH_loerr, nH_hierr, red_chi2, totcnts, num_bins) )
        results_file.close()

        # Finally, make sure all regions have an entry in the fit results file
        # (this shouldn't be needed, but just in case...)
        dtype = {'names': ('reg_id', 'kT', 'kT_lo', 'kT_hi', 'Z', 'Z_lo', 'Z_hi', 'norm', 'norm_lo', 'norm_hi', 'plindx', 'plindx_lo', 'plindx_hi', 'plnorm', 'plnorm_lo', 'plnorm_hi', 'nH', 'nH_lo', 'nH_hi', 'chi2', 'totcnts', 'nbins'), 'formats': ('i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4')}
        data = numpy.loadtxt(fit_results_file, dtype=dtype)
        reg_num = data["reg_id"]
        missing_regs = []
        n_missing = 0
        for i in range(len(reg_num)):
            if int(reg_num[i]) != i + reg_num_to_start + n_missing:
                missing_regs.append(i + reg_num_to_start + n_missing)
                n_missing += 1
        if n_missing > 0:
            results_file = open(fit_results_file, "a") # append missing regions
            for i in range(n_missing):
                results_file.write('%7r %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.4e %6.4e %6.4e %7.4f %7.4f %7.4f %6.4e %6.4e %6.4e %7.4f %7.4f %7.4f %7.4f %8.1f %8r\n' % (missing_regs[i]+reg_num_to_start, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0) )

    else:
        print('\n  Output file ('+fit_results_file+') exists and clobber = False.')


if __name__=='__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <spectra_list> <redshift> <nH_Gal> <kT_guess> <Ab_guess> <root>\n\nArguments:\n  <spectra_list>  input PI file (may be list; if so, prepend filename with "@")\n  <redshift>      redshift of source\n  <nH_Gal>        Galactic N_H (10^22 cm^-2)\n  <kT_guess>      initial guess for temperature\n  <Ab_guess>      initial guess for abundance\n  <root>          root of output file containing fit results', version="%prog 0.55")
    parser.add_option('--lo_energy', dest='lo_energy', help='lower energy bound for fit (keV); default = 0.5', metavar='VAL', default='0.5')
    parser.add_option('--hi_energy', dest='hi_energy', help='upper energy bound for fit (keV); default = 7.0', metavar='VAL', default='7.0')
    parser.add_option('--plasma_model', dest='plasma_model', help='plasma model to use in fit (mekal or apec); default = mekal', metavar='VAL', default='mekal')
    parser.add_option('--fix_nh', dest='fix_nH', help='Freeze nH; default = True', metavar='VAL', default=True)
    parser.add_option('--fix_abund', dest='fix_abund', help='Freeze abundance; default = False', metavar='VAL', default=False)
    parser.add_option('-c', action='store_true', dest='clobber', help='clobber any existing files', default=False)
    (options, args) = parser.parse_args()
    if len(args) == 6:
        spectra_list = args[0]
        redshift = args[1]
        nH_Gal = args[2]
        kT_guess = args[3]
        Ab_guess = args[4]
        root = args[5]
        lo_energy = options.lo_energy
        hi_energy = options.hi_energy
        plasma_model = options.plasma_model
        fix_nH_Gal = options.fix_nH
        fix_abund = options.fix_abund
        clobber = options.clobber

        # Read spectra file names from the spectra_list if it begins with '@'
        if spectra_list[0] == '@':
            spectra_list_file = open(spectra_list[1:], "r")
            spectra_list = spectra_list_file.readlines()
            spectra_list_file.close()
            for i in range(len(spectra_list)): spectra_list[i] = spectra_list[i].rstrip() # trim newlines
        else:
            if len(spectra_list) == 1: spectra_list = [spectra_list]

        # Call Sherpa
        call_sherpa_1T(spectra_list, redshift, nH_Gal, kT_guess, Ab_guess, root, lo_energy=lo_energy, hi_energy=hi_energy, plasma_model=plasma_model, fix_nH_Gal=fix_nH_Gal, fix_abund=fix_abund, clobber=clobber)

    else:
        parser.print_help()
