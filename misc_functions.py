# Miscellaneous helper functions
#
# Version 0.7: 1/11/2013 - Updated for Ciao 4.5


def stack_to_list(infile, stack_of_stacks=False, adjust_path=False):
    """
    Takes a input text file ('@' at the start denotes a stack) and returns a list
    """
    if infile[0] != '@':
        if ',' in infile:
            stack_list = infile.split(',')
            for i in range(len(stack_list)):
                stack_list[i] = stack_list[i].strip() # trim spaces
        else:
            stack_list = [infile]
    else:
        stack_file = open(infile[1:], "r")
        stack_list = stack_file.readlines()
        stack_file.close()
        for i, stack_entry in enumerate(stack_list):
            stack_list[i] = stack_entry.rstrip() # trim newlines
            if stack_of_stacks and ',' in stack_entry:
                stack_list[i] = stack_entry.split(',')
                for j, substack in enumerate(stack_list[i]):
                    stack_list[i][j] = substack.strip() # trim spaces
        blank_lines=True
        while blank_lines==True:
            try:
                stack_list.remove('')
            except ValueError:
                blank_lines=False

    # If the files do not have absolute paths, prepend "../" to each file
    # so that their relative paths work from the spectra subdirectory.
    if adjust_path:
        for i in range(len(stack_list)):
            if isinstance(stack_list[i], list):
                for j in range(len(stack_list[i])):
                    if stack_list[i][j][0] != '/': stack_list[i][j] = '../' + stack_list[i][j]
            else:
                if stack_list[i][0] != '/': stack_list[i] = '../' + stack_list[i]
    return stack_list


def combine_spectra(spectra_list, outroot, method='sum', bscale_method='asca',
    quiet=False):
    """
    Combines spectra for fitting using the CIAO tool combine_spectra

    The resulting combined spectrum is named {outroot}_src.pi

    Returns a list with the name of the resulting combined spectrum.
    """
    import subprocess

    cmd = ['punlearn', 'combine_spectra']
    p = subprocess.call(cmd)

    spectra_list_txt = ','.join(spectra_list)

    cmd = ['combine_spectra', spectra_list_txt, outroot, 'method='+method,
        'bscale_method='+bscale_method]

    if quiet:
        p = subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        p = subprocess.call(cmd)

    return [outroot + '_src.pi']


