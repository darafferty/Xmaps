"""
Takes a list of CIAO regions and extracts spectra and responses.

This module automates the extraction of spectra and responses
from Chandra data using the CIAO tool specextract. If run as a
script, multiple region files may be specfied by giving a text
file that lists the regions as an argument with '@' prepended.

The following routines are included:
  wextract - extract spectra and responses
  make_regions_from_binmap - makes CIAO region files defined by
                             a binmap
  transform_regions - transform CIAO region files defined using
                      one observation to another
  ad2xy - transform ra, dec to x, y image coordinates
  xy2ad - transform x, y image coordinates to ra, dec

Version 0.6: 3/5/2011 - Updated for Ciao 4.3
Version 0.7: 1/11/2013 - Updated for Ciao 4.5

"""
import os
import subprocess
import numpy
import multiprocessing
import threading
import sys
from misc_functions import stack_to_list
if os.environ.get("CALDB") == None:
    sys.exit('Please initalize CIAO before running.')

def wextract(region_list, evt2_file, pbk_file, asol_list, msk_file, bpix_file, root=None, bg_file=None, bg_region=None, binning=None, ncores=None, quiet=True, clobber=False):
    """
    Extracts spectra and makes appropriate weighted rmf and arf files.

    Inputs:  region_list - list of region file(s) in CIAO format (e.g. ['reg1.reg', 'reg2.reg'])
             evt2_file - events file
             pbk_file - pbk0 file
             asol_list - list of time-ordered asol1 file(s) (e.g. ['asol1_1', 'asol1_2'])
             msk_file - msk1 file
             bpix_file - bad pixel file
             root - root string for the created files (root_sou.pi, root_sou_grp.pi, etc.); if None,
                    it will be set to the root of each region file (without the '.reg').
             bg_file - background file(s); if None, a local background from the input evt2_file
                       will be made using the bg_region.
             bg_region - region to use for background; if None, the regions defined in the region file are used
             binning - number of counts per bin for output spectra
             ncores - number of cores to use (default = all)
             clobber - if True, overwrite any existing files
             quiet - if True, suppress all screen output from specextract

    """
    # Check that CIAO is initialized
    if os.environ.get("CALDB") == None:
        sys.exit('Please initalize CIAO before running.')

    nreg = len(region_list)
    if isinstance(region_list, str): # if not a list, make into list
        region_list = [region_list]
    if isinstance(asol_list, str): # if not a list, make into list
        asol_list = [asol_list]

    # Run asphist to make an aspect histogram for the input evt2 file
    # (specextract will do this, but if there are mulitple regions, it
    # is more efficient to do it once at the start)
    if len(asol_list) > 1:
        # make a stack
        asol_file = 'asol_stack.lst'
        asol_fileobj = open('asol_stack.lst', 'w')
        for asol in asol_list:
            asol_fileobj.write(asol+'\n')
        asol_fileobj.close()
    else:
        asol_file = asol_list[0]
    if root == None:
        reg_root = os.path.splitext(region_list[0])[0]
        indx_reg = reg_root.find('_reg1')
        if indx_reg > 0:
            asphist_file = reg_root[:indx_reg] + '_asphist.fits'
        else:
            asphist_file = reg_root + '_asphist.fits'
    else:
        asphist_file = root + '_asphist.fits'
    # Check if clobber is set
    if clobber == True:
        clobbertext = 'clobber=yes'
    else:
        clobbertext = ''
    cmd = ['asphist', asol_file, asphist_file, evt2_file, clobbertext]
    print('Making an aspect histogram ('+asphist_file+') for '+evt2_file+'...')
    p = subprocess.call(cmd)

    # Get number of cores if not specified and determine
    # the number of threads to use.
    # First, if quiet = False, use only one core for clarity;
    # otherwise, the output from multiple processes get mixed up.
    if quiet == False:
        ncores = 1
    if ncores == None or ncores > multiprocessing.cpu_count():
        ncores = multiprocessing.cpu_count()
    startindx = 0
    endindx = 0
    indxdel = int(nreg/ncores)
    if indxdel < 1:
        indxdel = 1
        nthreads = nreg
    else:
        nthreads = ncores

    # Now divide the input region list among the threads and do extraction
    threads = []
    # Setup the threads
    for i in range(nthreads):
        startindx = endindx
        endindx = startindx + indxdel
        if i == nthreads-1: endindx = nreg
        region_sublist = region_list[startindx:endindx]
        threads.append(threading.Thread(target=wextract_worker, args=(region_sublist, evt2_file, pbk_file, asphist_file, msk_file, bpix_file), kwargs={"root":root, "bg_file":bg_file, "bg_region":bg_region, "binning":binning, "clobber":clobber, "quiet":quiet, "pfiles_number":i}))

    # Start extraction
    for thr in threads:
        thr.start()

    # Wait until all regions are done before continuing
    for thr in threads:
        thr.join()

    # Check whether all spectra and responses were extracted successfully
    print('\nExtraction complete for '+evt2_file+':')
    made_all = True
    for i,region in enumerate(region_list):
        if root == None:
            pi_file = os.path.splitext(region)[0] + '_sou.pi'
        else:
            if nreg == 1:
                pi_file = root + '_sou.pi'
            else:
                pi_file = root + '_' + str(i) + '_sou.pi'
        if not os.path.isfile(pi_file):
            made_all = False
            print('  No spectra or responses made for region '+region+'.')
    if made_all == True:
        print('  All spectra and responses made successfully.')
    print(' ')


def wextract_worker(region_list, evt2_file, pbk_file, asphist_file, msk_file, bpix_file, root=None, bg_file=None, bg_region=None, binning=None, clobber=False, quiet=False, pfiles_number=None):
    """
    Worker script for wextract that allows multithreading.

    Inputs:  region_list - list of region file(s) in CIAO format (e.g. ['reg1.reg', 'reg2.reg'])
             evt2_file - events file
             pbk_file - pbk0 file
             asphist_file - aspect histogram file
             msk_file - msk1 file
             bpix_file - bad pixel file
             root - root string for the created files (root_sou.pi, root_sou_grp.pi, etc.); if None,
                    it will be set to the root of each region file (without the '.reg').
             bg_file - background file(s); if None, a local background from the input evt2_file
                       will be made using the bg_region.
             bg_region - region to use for background; if None, the regions defined in the region file are used
             binning - number of counts per bin for output spectra
             pfiles_number - PFILES directory number (e.g., 1 = cxcds_param1, etc.); used for multithreading only
             clobber - if True, overwrite any existing files
             quiet - if True, suppress all screen output

    """
    nreg = len(region_list)
    binarfwmap = 1 # set this to 2 or more if there are memory problems

    # Check if clobber is set
    if clobber == True:
        clobbertext = 'clobber=yes'
    else:
        clobbertext = ''

    # Check if binning is set
    if binning != None:
        binningtext1 = 'binspec='+str(binning)
        binningtext2 = 'grouptype=NUM_CTS'
    else:
        binningtext1 = 'binspec=NONE'
        binningtext2 = 'grouptype=NONE'

    # Check if pfiles_number is set. Specifying this ensures that multiple
    # CIAO instances run simultaneously without conflict by setting the PFILES
    # environment variable to a temporary directory (although it seems to
    # work fine even without it since pset is not currently used in this script).
    env = os.environ.copy()
    if pfiles_number != None:
        pfile_dir = './cxcds_param'+str(pfiles_number)
        if not os.path.exists(pfile_dir):
            os.mkdir(pfile_dir)
        # Set PFILES to a temporary directory in the current directory (note use of ; and :)
        env["PFILES"] = pfile_dir+';'+env["ASCDS_CONTRIB"]+'/param'+':'+env["ASCDS_INSTALL"]+'/param'

    # Now loop over regions and do extraction
    for i,region in enumerate(region_list):
        if os.path.isfile(region):
            print('Extracting spectra/responses from region ' + region + '...')

            # Define output file names (follow acisspec conventions)
            if root == None:
                pi_root = os.path.splitext(region)[0] + '_sou'
                pi_file = pi_root + '.pi'
                bg_pi_file = os.path.splitext(region)[0] + '_bgd.pi'
            else:
                if nreg == 1:
                    pi_root = root + '_sou'
                    pi_file = pi_root + '.pi'
                    bg_pi_file = root + '_bgd.pi'
                else:
                    pi_root = root + '_' + str(i) + '_sou'
                    pi_file = pi_root + '.pi'
                    bg_pi_file = root + str(i) + '_bgd.pi'

            # Extract source spectrum and make responses using SPECEXTRACT
            evt2_filter = evt2_file + '[sky=region(' + region + ')]'
            cmd = ['punlearn', 'specextract']
            p = subprocess.call(cmd, env=env)
            if bg_file == None and bg_region != None:
                bgtext = 'bkgfile=' + evt2_file + '[sky=region(' + bg_region + ')]'
            else:
                bgtext = 'bkgfile=NONE'

            cmd = ['specextract', evt2_filter, 'outroot='+pi_root, 'asp='+asphist_file, 'mskfile='+msk_file, 'badpixfile='+bpix_file, 'weight=yes', 'bkgresp=no', 'correct=no', 'combine=no', 'dafile=CALDB', bgtext, binningtext1, binningtext2, clobbertext, 'binarfwmap='+str(binarfwmap)]

            if quiet:
                p = subprocess.call(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            else:
                p = subprocess.call(cmd, env=env)

            # Extract background spectrum if a bg_file is given and a source spectrum was made
            if bg_file != None and os.path.isfile(pi_file):
                cmd = ['punlearn', 'dmextract']
                p = subprocess.call(cmd, env=env)
                bg_filter = bg_file + '[sky=region(' + region + ')][bin PI]'
                cmd = ['dmextract', 'infile=' + bg_filter, 'outfile=' + bg_pi_file, 'opt=pha1', 'error=gaussian', 'bkgerror=gaussian', 'wmap=[energy=300:2000][bin det=8]', clobbertext]
                if quiet:
                    p = subprocess.call(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # if quiet==True, simply redirect stdout and stderr to pipe
                else:
                    p = subprocess.call(cmd, env=env)

            #  Update header of source PI file with bg_file
            if os.path.isfile(pi_file):
                cmd = ['punlearn', 'dmhedit']
                p = subprocess.call(cmd, env=env)
                cmd = ['dmhedit', 'infile='+pi_file, 'filelist=', 'operation=add', 'key=BACKFILE', 'value='+bg_pi_file]
                if quiet:
                    p = subprocess.call(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                else:
                    p = subprocess.call(cmd, env=env)
                if binning != None:
                    cmd = ['dmhedit', 'infile='+pi_root+'_grp.pi', 'filelist=', 'operation=add', 'key=BACKFILE', 'value='+bg_pi_file]
                    if quiet:
                        p = subprocess.call(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    else:
                        p = subprocess.call(cmd, env=env)

        else:
            print('Region file ' + region + ' not found! No extraction done for this region.')


def make_regions_from_binmap(binmap_file, output_dir, reg_format='fits', minx=None, miny=None, bin=None, skip_dmimglasso=False, clobber=False):
    """
    Make CIAO region files from an input binmap.

    Inputs:  binmap - fits file of map of bins (pixel values = bin numbers)
             output_dir - name of output directory where region
                          files will be written
             reg_format - format of ouput region files ('fits' or 'ascii');
                          default = 'fits'
             minx - minimum x sky coordinate of binmap; default = None
             miny - minimum y sky coordinate of binmap; default = None
             bin - pixel binning used for binmap; default = None
             clobber - if True, overwrite any existing files

    Outputs: CIAO region files for each bin, suitable for spectral
             extraction. Output files are named "output_dir/reg*.fits"
             (for fits format) or "output_dir/reg*.reg" (for ascii format).

    Returns: A list of the region files created (without the output_dir
             prepended).

    Uses "dmimglasso" from CIAO to build the polygonal regions.

    """
    import astropy.io.fits as pyfits

    # Check if we need to find minx, miny, and binning of the binmap
    if minx == None or miny == None or bin == None:
        # Get pixel scale parameters out of the binmap header as needed
        hdr_binmap = pyfits.getheader(binmap_file)
        if minx == None:
            if 'CRVAL1P' in hdr_binmap:
                minx = hdr_binmap['CRVAL1P']
            else:
                sys.exit('ERROR: The binmap header does not have pixel coordinate information. Please specify the minx, miny, and binning for the binmap.')
        if miny == None:
            if 'CRVAL2P' in hdr_binmap:
                miny = hdr_binmap['CRVAL2P']
            else:
                sys.exit('ERROR: The binmap header does not have pixel coordinate information. Please specify the minx, miny, and binning for the binmap.')
        if bin == None:
            if 'CDELT1P' in hdr_binmap:
                binx = int(hdr_binmap['CDELT1P'])
            else:
                sys.exit('ERROR: The binmap header does not have pixel coordinate information. Please specify the minx, miny, and binning for the binmap.')
            if 'CDELT2P' in hdr_binmap:
                biny = int(hdr_binmap['CDELT2P'])
            else:
                biny = binx
    if bin != None:
        binx = bin
        biny = bin
    if not os.path.exists(output_dir):
        p = subprocess.call(['mkdir', output_dir])
    if clobber:
        p = subprocess.call(['rm', '-f', output_dir+'/bin_*.reg'])

    # Check if min bin is negative or starts or ends on the image boundary.
    # If so, assume it is not wanted (e.g., for wvt bin maps).
    binmap = pyfits.open(binmap_file, mode="readonly")
    binimage = binmap[0].data
    minbin = int(binimage.min())
    maxbin = int(binimage.max())
    if minbin < 0:
        minbin = 0
    inbin = numpy.where(binimage == minbin)
    if 0 in inbin[0] or numpy.size(binimage,0)-1 in inbin[0]:
        minbin += 1
    print('  Using minbin='+str(minbin)+', maxbin='+str(maxbin)+', minx='+str(minx)+', miny='+str(miny)+', binx='+str(binx)+', biny='+str(biny))

    # For each bin, construct region using CIAO's "dmimglasso"
    if not skip_dmimglasso:
        region_comment = '# Region file format: CIAO version 1.0\n'
        if clobber:
            clb_txt = 'yes'
        else:
            clb_txt = 'no'
        for i in range(minbin, maxbin+1):
            out_region_fits_file = output_dir+'/reg'+str(i)+'.fits'
            inbin = numpy.where(binimage == i)
            if len(inbin[0]) == 0:
                continue
            for j in range(len(inbin[0])):
                xpos = min(minx + binx * (inbin[1][j] + 1), minx + binx * binimage.shape[1])
                ypos = min(miny + biny * (inbin[0][j] + 1), miny + biny * binimage.shape[0])
                cmd = ['dmimglasso', binmap_file, out_region_fits_file, str(xpos), str(ypos), '0.1', '0.1', 'value=delta',     'maxdepth=1000000', 'clobber='+clb_txt]
                p = subprocess.call(cmd)
                if p == 0:
                    break

            #
            # Check for failure condition
            #
            # If dmimglasso fails, use ascii conversion work-around
            #
            if p != 0:
                cmd2 = ['dmmakereg', 'region('+output_dir+'/regions/xaf_'+str(i)+'.reg)', out_region_fits_file]
                q = subprocess.call(cmd2)

                print(i, xpos, ypos, out_region_fits_file, p)
                print(cmd)
                print(cmd2)
                print(q)


            if reg_format == 'ascii':
                if os.path.isfile(out_region_fits_file):
                    reg = pyfits.open(out_region_fits_file)
                    reg_params=reg[1].data.tolist()
                    xvertices = reg_params[0][0] # array of x coordinates of vertices
                    yvertices = reg_params[0][1] # array of y coordinates of vertices

                    out_region_file = open(output_dir+'/reg'+str(i)+'.reg', "w")
                    out_region_file.write(region_comment)

                    for j in range(len(xvertices)):
                        if j == 0:
                            region_text = 'polygon(%7.2f,%7.2f' % (xvertices[j], yvertices[j])
                        else:
                            region_text = region_text + ',%7.2f,%7.2f' % (xvertices[j], yvertices[j])
                    region_text = region_text + ')\n'
                    out_region_file.writelines(region_text)
                    out_region_file.close()
                    reg.close()
                    p = subprocess.call(['rm', '-f', out_region_fits_file])

    # Build region list
    bin_region_list = []
    for i in range(minbin, maxbin+1):
        if reg_format == 'ascii':
            # Check that each region file exists before adding it to the list
            if os.path.isfile(output_dir+'/reg'+str(i)+'.reg'):
                bin_region_list.append('reg'+str(i)+'.reg')
        else:
            # Check that each region file exists before adding it to the list
            filename = "reg%d.fits" % i
            path = os.path.join(output_dir, filename)
            if os.path.isfile(path) and pyfits.open(path)[1].data is not None:
                bin_region_list.append(filename)
            else:
                print("Warning: not using %s" % filename)
    return bin_region_list


def transform_regions(region_list, hdr_in, hdr_out, preroot, reg_format = 'fits', clobber = False):
    """
    Transforms CIAO region files for use with a different observation.

    Inputs:  region_list - list of region file(s) in CIAO format
             hdr_in - pyfits header of for input regions
             hdr_out - pyfits header of output regions
             preroot - prepend string for output region files
             reg_format - format of input/ouput region files ('fits' or 'ascii');
                          default = 'fits'
             clobber - if True, overwrite any existing files

    Ouputs:  New region files, prepended with the preroot.

    Returns: A list of the adjusted regions.

    """
    # Check if input and output headers are the same. If so,
    # just copy region files to new names.
    if hdr_in == hdr_out:
        new_region_list = []
        for region_file in region_list:
            if os.path.isfile(preroot + region_file) == False or clobber == True:
                if os.path.isfile(preroot + region_file):
                    p = subprocess.call(['rm', '-f', preroot + region_file])
                p = subprocess.call(['cp', region_file, preroot + region_file])
                new_region_list.append(preroot + region_file)
        return new_region_list

    # Get pixel scale parameters out of the input header
    # if it is the header of an image
    if 'CRVAL1P' in hdr_in:
        CRVAL1P_in = hdr_in['CRVAL1P']
        CRPIX1P_in = hdr_in['CRPIX1P']
        CDELT1P_in = hdr_in['CDELT1P']
        CRVAL2P_in = hdr_in['CRVAL2P']
        CRPIX2P_in = hdr_in['CRPIX2P']
        CDELT2P_in = hdr_in['CDELT2P']
        CRVAL1P_out = hdr_out['CRVAL1P']
        CRPIX1P_out = hdr_out['CRPIX1P']
        CDELT1P_out = hdr_out['CDELT1P']
        CRVAL2P_out = hdr_out['CRVAL2P']
        CRPIX2P_out = hdr_out['CRPIX2P']
        CDELT2P_out = hdr_out['CDELT2P']

    for region_file in region_list:
        # Read in regions
        if reg_format == 'ascii':
            cur_region_file = open(region_file, "r")
            region = cur_region_file.readlines()
            if region[0][0] == '#': # remove first line if it's a comment
                has_comment = True
                region_comment = region[0]
                region = region[1:]
            else:
                has_comment = False
            cur_region_file.close()
            nreg = len(region)
            for i in range(nreg): region[i] = region[i].rstrip() # trim newlines
        else:
            import astropy.io.fits as pyfits
            reg_orig = pyfits.open(region_file)
            region = [reg_orig[1].data] # region parameters are stored in second exten

        for r, reg in enumerate(region):
            # Find region type
            if reg_format == 'ascii':
                post0 = 0
                post1 = reg.find("(")
                reg_type = reg[post0:post1]
            else:
                reg_type = reg.SHAPE[0]

            if reg_type == 'polygon' or reg_type == 'Polygon':
                # Find coordinates of the vertices
                if reg_format == 'ascii':
                    coords = reg[post1+1:-1].split(',')
                    num_ver = len(coords)/2
                else:
                    num_ver = numpy.size(reg.X,1)

                for i in range(num_ver):
                    if reg_format == 'ascii':
                        xphys_in = float(coords[i*2])
                        yphys_in = float(coords[i*2+1])
                    else:
                        xphys_in = reg.X[0][i]
                        yphys_in = reg.Y[0][i]

                    # Change physical coordinates to image coordinates if needed
                    if 'CRVAL1P' in hdr_in:
                        xim_in = (xphys_in - CRVAL1P_in - CRPIX1P_in) / CDELT1P_in
                        yim_in = (yphys_in - CRVAL2P_in - CRPIX2P_in) / CDELT2P_in
                    else:
                        xim_in = xphys_in
                        yim_in = yphys_in

                    # Change image coordinates to ra, dec
                    a, d = xy2ad(xim_in, yim_in, hdr_in)

                    # Change ra, dec to output image coordinates
                    xim_out, yim_out = ad2xy(a, d, hdr_out)

                    # Change output image coordinates to physical coordinates if needed
                    if 'CRVAL1P' in hdr_in:
                        xphys_out=xim_out*CDELT1P_out+CRVAL1P_out+CRPIX1P_out
                        yphys_out=yim_out*CDELT2P_out+CRVAL2P_out+CRPIX2P_out
                    else:
                        xphys_out = xim_out
                        yphys_out = yim_out

                    # Replace physical coordinates in region with new ones
                    if reg_format == 'ascii':
                        if i==0:
                            new_coords_str = '%7.2f,%7.2f' % (xphys_out, yphys_out)
                        else:
                            new_coords_str = new_coords_str + ',%7.2f,%7.2f' % (xphys_out, yphys_out)
                        region[r] = reg[:post1+1] + new_coords_str + ')\n'
                    else:
                        reg.X[0][i] = xphys_out
                        reg.Y[0][i] = yphys_out

            if reg_type == 'rotbox':
                if reg_format == 'fits':
                    sys.exit('ERROR: Sorry, only polygons are currently supported for FITS regions.')

                # Find center(s) in physical coordinates
                posx0 = reg.find("(") + 1
                posx1 = reg.find(",")
                xphys_in = float(reg[posx0:posx1])
                posy0 = posx1 + 1
                posy1 = posy0 + reg[posy0:].find(",")
                yphys_in = float(reg[posy0:posy1])

                # Change physical coordinates to image coordinates if needed
                if 'CRVAL1P' in hdr_in:
                    xim_in = (xphys_in - CRVAL1P_in - CRPIX1P_in) / CDELT1P_in
                    yim_in = (yphys_in - CRVAL2P_in - CRPIX2P_in) / CDELT2P_in
                else:
                    xim_in = xphys_in
                    yim_in = yphys_in

                # Change image coordinates to ra, dec
                a, d = xy2ad(xim_in, yim_in, hdr_in)

                # Change ra, dec to output image coordinates
                xim_out, yim_out = ad2xy(a, d, hdr_out)

                # Change output image coordinates to physical coordinates if needed
                if 'CRVAL1P' in hdr_in:
                    xphys_out=xim_out*CDELT1P_out+CRVAL1P_out+CRPIX1P_out
                    yphys_out=yim_out*CDELT2P_out+CRVAL2P_out+CRPIX2P_out
                else:
                    xphys_out = xim_out
                    yphys_out = yim_out

                # Replace physical coordinates in region with new ones
                new_coords_str = '%7.2f,%7.2f' % (xphys_out, yphys_out)
                region[r] = reg[:posx0] + new_coords_str + reg[posy1:] + '\n'

        # Write new region file
        if reg_format == 'ascii':
            if os.path.isfile(preroot + region_file) == False or clobber == True:
                out_region_file = open(preroot + region_file, "w")
                if has_comment:
                    out_region_file.write(region_comment)
                out_region_file.writelines(region)
                out_region_file.close()
        else:
            if os.path.isfile(preroot + region_file) == False or clobber == True:
                reg_orig.writeto(preroot + region_file, clobber=True)

    new_region_list = []
    for i in range(len(region_list)):
        new_region_list.append(preroot + region_list[i])
    return new_region_list


def ad2xy(a, d, hdr):
    """
    Transforms ra, dec to image coordinates.

    Inputs:  a - ra in degrees
             d - dec in degrees
             hdr - pyfits header containing WCS info

    Returns: x, y image coordinates corresponding to a, d

    """
    if 'CRVAL1P' in hdr: # -> header of an image
        CRVAL1 = hdr['CRVAL1P']
        CRPIX1 = hdr['CRPIX1P']
        CDELT1 = hdr['CDELT1P']
        CRVAL2 = hdr['CRVAL2P']
        CRPIX2 = hdr['CRPIX2P']
        CDELT2 = hdr['CDELT2P']
    else: # -> header of an events file
        CRVAL1 = hdr['TCRVL11']
        CRPIX1 = hdr['TCRPX11']
        CDELT1 = hdr['TCDLT11']
        CRVAL2 = hdr['TCRVL12']
        CRPIX2 = hdr['TCRPX12']
        CDELT2 = hdr['TCDLT12']

    cd_mat = numpy.matrix([[1.0, 0.0], [0.0, 1.0]]) # assume no rotation for Chandra images
    cd_mat[0,0] = cd_mat[0,0] * CDELT1
    cd_mat[0,1] = cd_mat[0,1] * CDELT1
    cd_mat[1,1] = cd_mat[1,1] * CDELT2
    cd_mat[1,0] = cd_mat[1,0] * CDELT2

    radeg = 180.0 / numpy.pi
    h = numpy.sin(d / radeg) * numpy.sin(CRVAL2 / radeg) + numpy.cos(d / radeg) * numpy.cos(CRVAL2 / radeg) * numpy.cos(a / radeg - CRVAL1 / radeg)
    xsi = numpy.cos(d / radeg) * numpy.sin(a / radeg - CRVAL1 / radeg) / h * radeg
    eta = (numpy.sin(d / radeg) * numpy.cos(CRVAL2 / radeg) - numpy.cos(d / radeg) * numpy.sin(CRVAL2 / radeg) * numpy.cos(a / radeg - CRVAL1 / radeg)) / h * radeg
    cdinv = cd_mat.I
    xdif = cdinv[0,0] * xsi + cdinv[0,1] * eta
    ydif = cdinv[1,0] * xsi + cdinv[1,1] * eta
    x = xdif + CRPIX1 - 1
    y = ydif + CRPIX2 - 1

    return x, y


def xy2ad(x, y, hdr):
    """
    Transforms image coordinates to ra, dec.

    Inputs:  x - x image coordinate
             y - y image coordinate
             hdr - pyfits header containing WCS info

    Returns: ra, dec corresponding to x, y

    """
    if 'CRVAL1P' in hdr: # -> header of an image
        CRVAL1 = hdr['CRVAL1P']
        CRPIX1 = hdr['CRPIX1P']
        CDELT1 = hdr['CDELT1P']
        CRVAL2 = hdr['CRVAL2P']
        CRPIX2 = hdr['CRPIX2P']
        CDELT2 = hdr['CDELT2P']
    else: # -> header of an events file
        CRVAL1 = hdr['TCRVL11']
        CRPIX1 = hdr['TCRPX11']
        CDELT1 = hdr['TCDLT11']
        CRVAL2 = hdr['TCRVL12']
        CRPIX2 = hdr['TCRPX12']
        CDELT2 = hdr['TCDLT12']

    cd_mat = numpy.matrix([[1.0, 0.0], [0.0, 1.0]]) # assume no rotation for Chandra images
    cd_mat[0,0] = cd_mat[0,0] * CDELT1
    cd_mat[0,1] = cd_mat[0,1] * CDELT1
    cd_mat[1,1] = cd_mat[1,1] * CDELT2
    cd_mat[1,0] = cd_mat[1,0] * CDELT2

    xdif = x - (CRPIX1 - 1)
    ydif = y - (CRPIX2 - 1)
    xsi = cd_mat[0,0] * xdif + cd_mat[0,1] * ydif
    eta = cd_mat[1,0] * xdif + cd_mat[1,1] * ydif
    radeg = 180.0 / numpy.pi
    beta = numpy.cos(CRVAL2 / radeg) - eta / radeg * numpy.sin(CRVAL2 / radeg)
    a = (numpy.arctan2(xsi / radeg, beta) + CRVAL1 / radeg) * radeg
    gamma = numpy.sqrt((xsi / radeg)**2 + beta**2)
    d = numpy.arctan2(eta / radeg * numpy.cos(CRVAL2 / radeg) + numpy.sin(CRVAL2 / radeg), gamma) * radeg

    return a, d


if __name__=='__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <region_file> <evt2_file> <pbk_file> <asol_file> <msk_file>\n\nArguments:\n  <region_file>  region file in CIAO format (may be list; if so, prepend filename with "@")\n  <evt2_file>    events file\n  <pbk_file>     pbk0 file\n  <asol_file>    asol1 file (may be list; if so, prepend filename with "@")\n  <msk_file>     msk1 file', version="%prog 0.5")
    parser.add_option('--root', dest='root', help='root for output files; default = None', metavar='VAL', default=None)
    parser.add_option('--bg_region', dest='bg_region', help='background region file; default = None', metavar='FILE', default=None)
    parser.add_option('--bg_file', dest='bg_file', help='background fits file; default = None', metavar='FILE', default=None)
    parser.add_option('--bin', dest='bin_cnts', help='number of counts per spectral bin; default = None', metavar='VAL', default=None)
    parser.add_option('--ncores', dest='ncores', help='number of cores to use; default = all', metavar='VAL', default=None)
    parser.add_option('-c', action='store_true', dest='clobber', help='clobber any existing files', default=False)
    (options, args) = parser.parse_args()
    if len(args) == 5:
        region_file = args[0]
        evt2_file = args[1]
        pbk_file = args[2]
        asol_file = args[3]
        msk_file = args[4]
        root = options.root
        bg_file = options.bg_file
        bg_region = options.bg_region
        clobber = options.clobber
        binning = options.bin_cnts
        ncores = options.ncores

        # Read region file names from the region_file if region_file begins with '@'
        if region_file[0] == '@':
            region_list_file = open(region_file[1:], "r")
            region_list = region_list_file.readlines()
            region_list_file.close()
            for i in range(len(region_list)): region_list[i] = region_list[i].rstrip() # trim newlines
        else:
            region_list = [region_file]

        asol_list = stack_to_list(asol_file)

        # Do extraction using threads
        wextract(region_list, evt2_file, pbk_file, asol_list, msk_file, root=root, bg_file=bg_file, bg_region=bg_region, binning=binning, ncores=ncores, clobber=clobber)
    else:
        parser.print_help()
