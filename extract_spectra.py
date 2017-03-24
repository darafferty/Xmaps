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

Version 0.6: 3/5/2011 - Updated for Ciao 4.3
Version 0.7: 1/11/2013 - Updated for Ciao 4.5

"""
import os
import subprocess
import numpy
import multiprocessing
import threading
import sys
import shutil

import pycrates
import pytransform
import stk

from misc_functions import stack_to_list

if os.environ.get("CALDB") is None:
    sys.exit('Please initalize CIAO before running.')


def wextract(region_list, evt2_file, pbk_file, asol_list, msk_file,
             bpix_file, root=None, bg_file=None, bg_region=None,
             binning=None, ncores=None, quiet=True, clobber=False):
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
    if os.environ.get("CALDB") is None:
        sys.exit('Please initalize CIAO before running.')

    nreg = len(region_list)
    if isinstance(region_list, str):  # if not a list, make into list
        region_list = [region_list]
    if isinstance(asol_list, str):  # if not a list, make into list
        asol_list = [asol_list]

    # Run asphist to make an aspect histogram for the input evt2 file
    # (specextract will do this, but if there are mulitple regions, it
    # is more efficient to do it once at the start)
    if len(asol_list) > 1:
        # make a stack
        asol_file = 'asol_stack.lst'
        asol_fileobj = open('asol_stack.lst', 'w')
        for asol in asol_list:
            asol_fileobj.write(asol + '\n')
        asol_fileobj.close()
    else:
        asol_file = asol_list[0]
    if root is None:
        reg_root = os.path.splitext(region_list[0])[0]
        indx_reg = reg_root.find('_reg1')
        if indx_reg > 0:
            asphist_file = reg_root[:indx_reg] + '_asphist.fits'
        else:
            asphist_file = reg_root + '_asphist.fits'
    else:
        asphist_file = root + '_asphist.fits'
    # Check if clobber is set
    if clobber:
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
    if not quiet:
        ncores = 1
    if ncores is None or ncores > multiprocessing.cpu_count():
        ncores = multiprocessing.cpu_count()
    startindx = 0
    endindx = 0
    indxdel = int(nreg / ncores)
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
        if i == nthreads - 1:
            endindx = nreg
        region_sublist = region_list[startindx:endindx]
        threads.append(threading.Thread(target=wextract_worker, args=(region_sublist, evt2_file, pbk_file, asphist_file, msk_file, bpix_file), kwargs={"root":root, "bg_file":bg_file, "bg_region":bg_region, "binning":binning, "clobber":clobber, "quiet":quiet, "pfiles_number":i}))

    # Start extraction
    for thr in threads:
        thr.start()

    # Wait until all regions are done before continuing
    for thr in threads:
        thr.join()

    # Check whether all spectra and responses were extracted successfully
    print('\nExtraction complete for ' + evt2_file + ':')
    made_all = True
    for i, region in enumerate(region_list):
        if root is None:
            pi_file = os.path.splitext(region)[0] + '_sou.pi'
        else:
            if nreg == 1:
                pi_file = root + '_sou.pi'
            else:
                pi_file = root + '_' + str(i) + '_sou.pi'
        if not os.path.isfile(pi_file):
            made_all = False
            print('  No spectra or responses made for region ' +
                  region + '.')
    if made_all:
        print('  All spectra and responses made successfully.')
    print(' ')


def wextract_worker(region_list, evt2_file, pbk_file, asphist_file,
                    msk_file, bpix_file, root=None, bg_file=None,
                    bg_region=None, binning=None, clobber=False,
                    quiet=False, pfiles_number=None):
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
    binarfwmap = 1  # set this to 2 or more if there are memory problems

    # Check if clobber is set
    if clobber:
        clobbertext = 'clobber=yes'
    else:
        clobbertext = ''

    # Check if binning is set
    if binning is not None:
        binningtext1 = 'binspec=' + str(binning)
        binningtext2 = 'grouptype=NUM_CTS'
    else:
        binningtext1 = 'binspec=NONE'
        binningtext2 = 'grouptype=NONE'

    # Check if pfiles_number is set. Specifying this ensures that multiple
    # CIAO instances run simultaneously without conflict by setting the PFILES
    # environment variable to a temporary directory (although it seems to
    # work fine even without it since pset is not currently used in this script).
    env = os.environ.copy()
    if pfiles_number is not None:
        pfile_dir = './cxcds_param' + str(pfiles_number)
        if not os.path.exists(pfile_dir):
            os.mkdir(pfile_dir)
        # Set PFILES to a temporary directory in the current directory (note use of ; and :)
        env["PFILES"] = pfile_dir+';'+env["ASCDS_CONTRIB"]+'/param'+':'+env["ASCDS_INSTALL"]+'/param'

    # Now loop over regions and do extraction
    for i, region in enumerate(region_list):
        if os.path.isfile(region):
            print('Extracting spectra/responses from region ' + region + '...')

            # Define output file names (follow acisspec conventions)
            if root is None:
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
            if bg_file is None and bg_region is not None:
                bgtext = 'bkgfile=' + evt2_file + \
                         '[sky=region(' + bg_region + ')]'
            else:
                bgtext = 'bkgfile=NONE'

            cmd = ['specextract', evt2_filter, 'outroot='+pi_root, 'asp='+asphist_file, 'mskfile='+msk_file, 'badpixfile='+bpix_file, 'weight=yes', 'bkgresp=no', 'correct=no', 'combine=no', 'dafile=CALDB', bgtext, binningtext1, binningtext2, clobbertext, 'binarfwmap='+str(binarfwmap)]

            if quiet:
                p = subprocess.call(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            else:
                p = subprocess.call(cmd, env=env)

            # Extract background spectrum if a bg_file is given and a source spectrum was made
            if bg_file is not None and os.path.isfile(pi_file):
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
                if binning is not None:
                    cmd = ['dmhedit', 'infile='+pi_root+'_grp.pi', 'filelist=', 'operation=add', 'key=BACKFILE', 'value='+bg_pi_file]
                    if quiet:
                        p = subprocess.call(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    else:
                        p = subprocess.call(cmd, env=env)

        else:
            print('Region file ' + region + ' not found! No extraction done for this region.')


def find_sky_transform(cr):
    """Return the SKY transform for the input image crate.

    Parameters
    ----------
    cr
        A IMAGECrate, e.g. as returned by pycrates.read_file().
        It could also support tables, but that is trickier to do,
        so currently restricted to an image.

    Returns
    -------
    tr : None or a LINEAR2DTransform
        Returns None if no transform can be found, otherwise
        it returns a 2D transform.

    Notes
    -----
    The mapping from FITS logical coordinates - what DS9 calls
    "Image" coordinates - to the Chandra SKY system - which DS9
    calls "Physical" coordinates - can be stored as two 1D
    transforms (since Chandra data has no rotation) or a 2D
    transform. This routine tries to "hide" this difference
    and just return a 2D transform, whatever is in the input file.

    If the files all matched Chandra data products, then it would
    be relatively easy, but without knowing the axis names it is
    rather messy. An alternative would be to manually add keywords
    to the image files to make them look more like Chandra images,
    and avoid the complexity here. This would include keyword/values
    like

        MTYPE1 = 'SKY'
        MFORM1 = 'X,Y'
        MTYPE2 = 'EQPOS'
        MFORM2 = 'RA,DEC'

    as well as the CRVAL1/2P, CRPIX1/2P, and CDELT1/2P keywords
    (for the SKY coordinate transform) and CRVAL1/2, CRPIX1/2, and
    CDELT1/2 for the celestial coordinates.

    I expect this routine to fail horribly once it is used on
    actual data.

    """

    assert isinstance(cr, pycrates.IMAGECrate)

    axes = cr.get_axisnames()
    if axes == []:
        # Assume that if there's no axis names then there's no
        # transform data.
        return None

    # The simple case is if there's a "vector" column (e.g. SKY)
    # with a 2D transform.
    #
    try:
        tr = cr.get_transform(axes[0])
    except KeyError:
        # For now return None, but could try something like
        # iterating through the other axis names, if there
        # are any.
        return None

    if isinstance(tr, pytransform.LINEAR2DTransform):
        return tr

    elif not isinstance(tr, pytransform.LINEARTransform):
        # For now return None, but could try something like
        # iterating through the other axis names, if there
        # are any.
        return None

    # Assume that the second component is the second
    # axis.
    xtr = tr
    try:
        ytr = cr.get_transform(axes[1])
    except KeyError:
        return None

    # Create a 2D transform based on the two 1D transforms.
    #
    trs = [xtr, ytr]
    scales = [itr.get_parameter('SCALE').get_value()
              for itr in trs]
    offsets = [itr.get_parameter('OFFSET').get_value()
               for itr in trs]
    out = pytransform.LINEAR2DTransform()
    out.get_parameter('ROTATION').set_value(0)
    out.get_parameter('SCALE').set_value(scales)
    out.get_parameter('OFFSET').set_value(offsets)
    return out


def make_regions_from_binmap(binmap_file, output_dir,
                             reg_format='fits',
                             minx=None, miny=None, bin=None,
                             skip_dmimglasso=False,
                             clobber=False):
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

    cr = pycrates.read_file(binmap_file)
    file_sky_transform = find_sky_transform(cr)

    # Create the mapping from logical (image or pixel) coordinates
    # to the SKY system. If none of minx, miny, or bin are given
    # then the transform from the file is used, otherwise a new
    # transform is created using the user-supplied information.
    #
    # This is not quite the same as the original version. In part,
    # the original code was a little-more general, in that it
    # could read in "partial" data on the SKY coordinate
    # transformation - e.g. if the file only had CDELT1P and CDELT2P
    # keywords then the original version would pick this up. Using
    # the crates approach is more of an all-or-nothing: you either
    # have all the keywords or none of them.
    #
    if minx is None and miny is None and bin is None:
        sky_transform = file_sky_transform

    else:
        if file_sky_transform is None:
            sys.exit('ERROR: The binmap header does not have ' +
                     'pixel coordinate information. Please specify ' +
                     'the minx, miny, and binning for the binmap.')

        scales = file_sky_transform.get_parameter('SCALE').get_value()
        offsets = file_sky_transform.get_parameter('OFFSET').get_value()

        # There is a slight difference here in how the OFFSET values
        # are processed (these correspond to the CRVAL1/2P values).
        # The transform offsets are defined for the logical coordinate
        # (0, 0), whereas the file values are defined for whatever
        # the CRPIX1/2P values are. Normally these are (0.5,0.5),
        # and I believe the original code rounded the CRVALXP values
        # down, which effectively matches things up.
        #
        if minx is not None:
            offsets[0] = minx * 1.0
        if miny is not None:
            offsets[1] = miny * 1.0

        if bin is not None:
            scales = [bin * 1.0, bin * 1.0]

        sky_transform = pytransform.LINEAR2DTransform()
        sky_transform.get_parameter('ROTATION').set_value(0)
        sky_transform.get_parameter('SCALE').set_value(scales)
        sky_transform.get_parameter('OFFSET').set_value(offsets)

    # This logic could be moved into the if statement above to
    # avoid unneeded work, but it's clearer here.
    scales = sky_transform.get_parameter('SCALE').get_value()
    offsets = sky_transform.get_parameter('SCALE').get_value()
    binx, biny = scales
    minx, miny = offsets
    file_sky_transform = None

    if not os.path.exists(output_dir):
        p = subprocess.call(['mkdir', output_dir])

    if clobber:
        p = subprocess.call(['rm', '-f', output_dir + '/bin_*.reg'])

    # Check if min bin is negative or starts or ends on the image
    # boundary. If so, assume it is not wanted (e.g., for wvt bin
    # maps).
    binimage = cr.get_image().values
    minbin = int(binimage.min())
    maxbin = int(binimage.max())
    if minbin < 0:
        minbin = 0

    inbin = numpy.where(binimage == minbin)
    if 0 in inbin[0] or numpy.size(binimage, 0) - 1 in inbin[0]:
        minbin += 1

    print('  Using minbin=' + str(minbin) +
          ', maxbin=' + str(maxbin) +
          ', minx=' + str(minx) +
          ', miny=' + str(miny) +
          ', binx=' + str(binx) +
          ', biny=' + str(biny))

    # For each bin, construct region using CIAO's "dmimglasso"
    #
    # The coordinates returned by numpy.where are 0 based, but
    # the FITS logical/image coordinate system is 1 based, so
    # a conversion is needed when passing to sky_transform.
    #
    if not skip_dmimglasso:
        region_comment = '# Region file format: CIAO version 1.0\n'
        if clobber:
            clb_txt = 'yes'
        else:
            clb_txt = 'no'

        for i in range(minbin, maxbin + 1):
            out_region_fits_file = output_dir + '/reg' + \
                str(i) + '.fits'
            inybin, inxbin = numpy.where(binimage == i)
            if len(inybin) == 0:
                continue

            # Convert the j,i values from where into the FITS
            # logical coordinate system lx,ly and then convert
            # to SKY values. It is important that lcoords is
            # a floating-point value, and not an integer one,
            # to ensure that no truncation of the result happens.
            #
            # An alternative would be to set coord=logical
            # when running dmimglasso, so that
            # inxbin+1, inybin+1 could be used (i.e. no need
            # for the coordinate conversion.
            #
            lcoords = numpy.vstack((inxbin + 1.0, inybin + 1.0)).T
            scoords = sky_transform.apply(lcoords)

            for xpos, ypos in scoords:
                # Is this restriction needed?
                xpos = min(xpos, minx + binx * binimage.shape[1])
                ypos = min(ypos, miny + biny * binimage.shape[0])
                cmd = ['dmimglasso', binmap_file,
                       out_region_fits_file,
                       str(xpos), str(ypos), '0.1', '0.1',
                       'value=delta', 'maxdepth=1000000',
                       'clobber=' + clb_txt]
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
                    reg = pycrates.read_file(out_region_fits_file)
                    vertices = reg.get_column(0).values
                    xvertices = vertices[0, 0]
                    yvertices = vertices[0, 1]

                    out_region_file = open(output_dir + '/reg' +
                                           str(i) + '.reg', "w")
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
    for i in range(minbin, maxbin + 1):
        if reg_format == 'ascii':
            # Check that each region file exists before adding it to the list
            rname = 'reg' + str(i) + '.reg'
            if os.path.isfile(output_dir + '/' + rname):
                bin_region_list.append(rname)
        else:
            # Check that each region file exists before adding it to the list
            filename = "reg%d.fits" % i
            path = os.path.join(output_dir, filename)
            # It is not 100% clear what the equivalent to checking
            # `pyfits.open(path)[1].data is not None`. I am going
            # to use a check for a non-zero number of rows, forcing
            # the second block (CXC Data model starts counting at
            # a block number of 1, not 0).
            if os.path.isfile(path) and \
                    pycrates.read_file(path + "[2]").get_nrows() > 0:
                bin_region_list.append(filename)
            else:
                print("Warning: not using %s" % filename)
    return bin_region_list


def copy_regions(region_list, preroot, clobber=False):
    """Copy regions to new names.

    This is for when the input and output coordinate coordinate
    systems are the same.
    """

    new_region_list = []
    for region_file in region_list:
        newname = preroot + region_file
        if clobber or not os.path.isfile(newname):
            if os.path.isfile(newname):
                subprocess.check_call(['rm', '-f', newname])

            shutil.copy(region_file, newname)
            new_region_list.append(newname)

    return new_region_list


def transform_region_ascii(infile, outfile, wcs_in, wcs_out):
    """Transform coordinates in ASCII input file.

    This routine is called assuming the output file is to
    be clobbered if it already exists.
    """

    with open(infile, 'r') as fh:
        regions = fh.readlines()

    with open(outfile, 'w') as ofh:
        for region in regions:
            if region.startswith('#'):
                ofh.write(region + '\n')
                continue

            region = region.rstrip()
            post0 = 0
            post1 = region.find("(")
            reg_type = region[post0:post1]

            if reg_type in ['polygon', 'Polygon']:
                # convert from a 1D array into a 2D one
                coords_in = [float(f)
                             for f in region[post1 + 1:-1].split(',')]

                assert coords_in.size % 2 == 0
                # Use integer division here
                coords_in.resize(2, coords_in.size // 2)

                # The conversion can be applied to all the
                # pairs at once, but it requires the data be
                # in the "right" shape.
                #
                coords_cel = wcs_in.apply(coords_in.T)
                coords_out = wcs_out.invert(coords_cel)

                # The coords_out array is not transposed (to
                # match the input) since it makes it easier
                # to convert back to a string.
                coords_str = ",".join(["{:7.2f}".format(c)
                                       for c in coords_out])

                out = reg_type + '(' + coords_str + ')'

            elif reg_type == 'rotbox':

                # Just need to convert the center of the box, since
                # the assumption is that the pixel scale is the
                # same in both the input and output systems.
                #
                toks = region[post1 + 1:].split(",")
                assert len(toks) > 2

                xphys_in = float(toks[0])
                yphys_in = float(toks[1])

                # The handling of nD arrays by the apply and invert
                # methods of transform objects is, at best, strange
                # to describe.
                #
                coords_cel = wcs_in.apply([[xphys_in, yphys_in]])
                coords_out = wcs_out.invert(coords_cel)

                xphys_out = coords_out[0][0]
                yphys_out = coords_out[0][1]
                coords_str = '{:7.2f},{:7.2f},'.format(xphys_out,
                                                       yphys_out)

                # Hopefully this re-creates the remainded of the
                # string (i.e. after the center of the box).
                #
                out = reg_type + '(' + coords_str + ",".join(toks[2:])

            else:
                # copy over the line
                out = region

            ofh.write(out + '\n')


def transform_region_fits(infile, outfile, wcs_in, wcs_out):
    """Transform coordinates in FITS input file.

    This routine is called assuming the output file is to
    be clobbered if it already exists.
    """

    # Read in the whole file in case there are any other interesting
    # blocks in the file that should be retained.
    #
    ds = pycrates.CrateDataset(infile, mode='r')
    cr = ds.get_crate(2)
    assert isinstance(cr, pycrates.TABLECrate)

    # NOTE: the EQPOS column could be read directly, which would
    # avoid the need to convert the input file from physical to
    # celestial coordinates.
    #
    shapes = cr.get_column('SHAPE').values
    pos = cr.get_column('POS').values

    # Since the shapes are being processed on a row-by-row bases,
    # do not try and convert all the coordinates in one go (which
    # is supported by the apply and invert methods).
    #
    for i in xrange(0, cr.get_nrows()):

        shape = shapes[i].upper()
        if shape == 'POLYGON':

            coords_cel = wcs_in.apply(pos[i].T)
            coords_out = wcs_out.invert(coords_cel).T

            # Overwrite this row
            pos[i] = coords_out

        elif shape == 'ROTBOX':
            # It should be possible to encode ROTBOX in CXC FITS
            # region files.
            #
            sys.exit('ERROR: Sorry, only polygons are currently supported for FITS regions.')

    cr.get_column('POS').values = pos
    ds.write(outfile, clobber=True)


def transform_regions(region_list, wcs_in, wcs_out, preroot,
                      reg_format='fits', clobber=False):
    """
    Transforms CIAO region files for use with a different observation.

    Inputs:  region_list - list of region file(s) in CIAO format
             wcs_in - EQPOS transform for input coordinate system
                      (physical to celestial)
             wcs_out - EQPOS transform for output coordinate system
                       (physical to celestial)
             preroot - prepend string for output region files
             reg_format - format of input/ouput region files
                          ('fits' or 'ascii'); default = 'fits'
             clobber - if True, overwrite any existing files

    Ouputs:  New region files, prepended with the preroot.

    Returns: A list of the adjusted regions.

    It is assumed that sky_in and sky_out are different; if they
    are the same use copy_regions().

    """

    if reg_format == 'ascii':
        func = transform_region_ascii
    elif reg_format == 'fits':
        func = transform_region_fits
    else:
        raise ValueError("Unsupported reg_format=" + reg_format)

    outfiles = []
    for region_file in region_list:

        # Do not do anything if the output file already exists and
        # clobber is not set.
        #
        outfile = preroot + region_file
        outfiles.append(outfile)
        if not clobber and os.path.exists(outfile):
            continue

        func(region_file, outfile, wcs_in, wcs_out)

    return outfiles


if __name__ == '__main__':

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
        # stk.build reads in a stack - e.g. a single value, or a
        # comma-separated list of names, or a filename with a leading
        # '@' and returns a list of values.
        region_list = stk.build(args[0])
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

        # TODO: stack_to_list could be simplified by using stk.build
        #       but leave that for a later date.
        asol_list = stack_to_list(asol_file)

        # Do extraction using threads
        wextract(region_list, evt2_file, pbk_file, asol_list, msk_file,
                 root=root, bg_file=bg_file, bg_region=bg_region,
                 binning=binning, ncores=ncores, clobber=clobber)
    else:
        parser.print_help()
