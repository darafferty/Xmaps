Xmaps
=====

Python scripts to make maps of spectral parameters with Chandra data.

The make_maps.py script will extract and fit spectra from regions defined
in an input binmap. The binmap is simply a FITS image with the values that
lie in each bin equal to the bin's number. Binmaps can be generated in a
number of ways, such as with the Contbin algorithm of Sanders (2006) or
the WVT algorithm of Diehl & Statler (2007).

CIAO must be initialized before starting the script.

For usage information, run the make\_maps.py script without arguments or
as "make_maps.py --help".

Available models
----------------

By default, the script fits a WABS * MEKAL model to each spectrum. A second
component can be added using the "--add\_comp=component" argument, where
component is one of the following:

*   mekal/apec: add a 2nd single-temperature component (WABS*[MEKAL + MEKAL]).
    The abundances of the two single-temperature components are tied together.
*   pow: add a power-law component (WABS*[MEKAL + POW])
*   mkcflow: add a MKCFLOW component (WABS*[MEKAL + MKCFLOW])

If a second component is specified, the fit is performed twice: once for the
single-temperature model and again with the second component included. The chi2
values for the two fits are compared in each region. If the fit for a given
region improves significantly with the second component included (from an
F-test with a null probability of 0.05), the two-component fit is used in the
final maps; otherwise, the single-temperature results are used (if the mapped
parameter is not one of the single-temperature model, its value is set to 0.0
when the single-temperature model is preferred).

In the case of the two-temperature model, two maps are made: a low-temperature
and a high-temperature map. If the single-temperature fit is preferred in this
case, its value of the temperature is used in both maps. In all cases when a
second component is specified, an F-test probability map is made.

Spectral parameters that can be mapped
--------------------------------------

The following parameters can be mapped (by passing them to make\_maps.py with
the --vars\_to_map argument; e.g., "--vars\_to\_map="kT, Z, fkT"):

*   kT / fkT: MEKAL or APEC temperature / error (keV)
*   Z / fZ: MEKAL or APEC abundance / error (Solar)
*   nH / fnH: WABS column density / error (10^22 cm^-2)
*   norm / fnorm: MEKAL or APEC normalization / error
*   chi2: chi-squared value of the fit
*   plindx / fplindx: power-law index / error; requires "--add\_comp=pow"
*   mdot / fmdot: MKCFLOW cooling rate / error (solar mass/yr); requires "--add\_comp=mkcflow"

