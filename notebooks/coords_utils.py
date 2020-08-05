"""
 Taken from Erin Sheldon's esutil package at:
 https://github.com/esheldon/esutil/blob/master/esutil/coords.py
 2020
    NAME
        coords
    PURPOSE
        A set of astronomical utilities for dealing with coordinates and
        coordinate transformations.
    COORDINATE TRANSFORMATIONS
        euler:
            A generic routine for transforming between Galactic, Celestial,
            and ecliptic coords.  The following wrapper routines are also
            supplied for convenience:
        l,b = eq2gal(ra, dec, b1950=False, dtype='f8')
            Convert equatorial to glactic coordinates.
        # The following use the same interface:
        gal2eq
            Convert galactic to equatorial coordinates.
        eq2ec
            Convert equatorial to ecliptic coordinates.
        ec2eq
            Convert ecliptic to equatorial coordinates.
        ec2gal
            Convert ecliptic to galactic coordinates.
        gal2ec
            Convert galactic to ecliptic coordinates.
        # These SDSS specific functions do not use euler
        eq2sdss
            Convert between equatorial and corrected SDSS survey coords.
        sdss2eq
            Convert between corrected SDSS survey and equatorial coords.
        eq2xyz: Convert equatorial to x,y,z on the sphere according to
            the following transform:
                    x = sin(pi/2-dec)*cos(ra)
                    y = sin(pi/2-dec)*sin(ra)
                    z = cos(pi/2-dec)
        xyz2eq:
            inverse of eq2xyz
        sphdist:
            Calculate the arc length between two sets of points on the sphere.
            Currently only takes ra,dec.
        shiftlon:
            shift the input longitude.  By default wrap the coordinate to
            -180,180.  If a shift is entered, return the new value
            lon-shift such that the range is still [0,360)
        shiftra:
            shift right ascension.  This just calls shiftlon
        radec2aitoff:
            Convert ra,dec to aitoff coordinates.
        dec_parse(decstring)
            parse a colon separated string representing declination ito
            degrees.
        ra_parse(decstring)
            parse a colon separated string representing right ascension ito
            degrees.
        randsphere(numrand, system='eq', ra_range=[0,360], dec_range=[-90,90]):
            Generate random points on the sphere.  By default ra,dec are
            returned.  If system='xyz' then x,y,z are returned.
        randcap(nrand,ra,dec,rad,get_radius=False):
            Create random points in a cap, or disc, centered at the
            input ra,dec location and with radius rad.
        rect_area(lon_min, lon_max, lat_min, lat_max)
            Calculate the area of a rectangle on the sphere.
"""
import math
import numpy as np

try:
    import numpy
    from numpy import where, sin, cos, arccos, arcsin, arctan2, sqrt, rad2deg, deg2rad, zeros
    have_numpy = True
except:
    have_numpy = False


_license = """
  Copyright (C) 2009  Erin Sheldon
    This program is free software; you can redistribute it and/or modify it
    under the terms of version 2 of the GNU General Public License as
    published by the Free Software Foundation.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

PI = math.pi
HALFPI = PI / 2.0
D2R = PI / 180.0
R2D = 1.0 / D2R

_centra = 185.
_centdec = 32.5
_sdsspar = {'center_ra': _centra,
            'center_dec': _centdec,
            'node': (_centra - 90.0) * D2R,
            'etapole': _centdec * D2R,
            'etaoffset': 91.25,
            'doc': """
    A set of transformation functions for use with SDSS coordinate systems.
    eq2sdss(): Convert between equatorial and corrected SDSS survey coords.
    sdss2eq(): Convert between corrected SDSS survey and equatorial coords.
    Don't use these unless you have to, as these are the old coordinates
        eq2survey(): Convert between equatorial and uncorrected survey coords.
        survey2eq(): Convert between uncorrected survey and equatorial coords.
    Adapted from astrotools
        Erin Sheldon, NYU, 2006-03-11
    Force data type and allow selection of dtype through keyword.
        Erin Sheldon, NYU, 2007-05-23
"""
            }

def euler(ai_in, bi_in, select, b1950=False, dtype='f8'):
    """
    NAME:
        euler
    PURPOSE:
        Transform between Galactic, celestial, and ecliptic coordinates.
    CALLING SEQUENCE:
        long_out, lat_out =
            euler(long_in, lat_in, type, b1950=False, dtype='f8')
    INPUTS:
       long_in - Input Longitude in DEGREES, scalar or vector.
       lat_in  - Input Latitude in DEGREES
       select  - Integer (1-6) specifying type of coordinate transformation.
      select   From          To        |   select      From            To
       1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec
       2     Galactic       RA-DEC     |     5       Ecliptic      Galactic
       3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic
      Celestial coordinates (RA, Dec) should be given in equinox J2000
      unless the b1950=True keyword is set.
    OUTPUTS:
       long_out - Output Longitude in DEGREES
       lat_out  - Output Latitude in DEGREES
    INPUT KEYWORD:
       b1950 - If this keyword is true then input and output
             celestial and ecliptic coordinates should be given in equinox
             B1950.
    REVISION HISTORY:
       Written W. Landsman,  February 1987
       Adapted from Fortran by Daryl Yentis NRL
       Converted to IDL V5.0   W. Landsman   September 1997
       Made J2000 the default, added /FK4 keyword  W. Landsman December 1998
       Add option to specify SELECT as a keyword W. Landsman March 2003
       Converted from IDL to numerical Python: Erin Sheldon, NYU, 2008-07-02
    """

    # Make a copy as an array. ndmin=1 to avoid messed up scalar arrays
    ai = numpy.array(ai_in, ndmin=1, copy=True, dtype=dtype)
    bi = numpy.array(bi_in, ndmin=1, copy=True, dtype=dtype)

    twopi = 2.0 * PI
    fourpi = 4.0 * PI

    #   J2000 coordinate conversions are based on the following constants
    #   (see the Hipparcos explanatory supplement).
    #  eps = 23.4392911111d           Obliquity of the ecliptic
    #  alphaG = 192.85948d            Right Ascension of Galactic North Pole
    #  deltaG = 27.12825d             Declination of Galactic North Pole
    #  lomega = 32.93192d             Galactic longitude of celestial equator
    #  alphaE = 180.02322d            Ecliptic longitude of Galactic North Pole
    #  deltaE = 29.811438523d         Ecliptic latitude of Galactic North Pole
    #  Eomega  = 6.3839743d           Galactic longitude of ecliptic equator
    # Parameters for all the different conversions
    if b1950:
        # equinox = '(B1950)'
        psi = numpy.array([0.57595865315, 4.9261918136,
                           0.00000000000, 0.0000000000,
                           0.11129056012, 4.7005372834], dtype=dtype)
        stheta = numpy.array([0.88781538514, -0.88781538514,
                              0.39788119938, -0.39788119938,
                              0.86766174755, -0.86766174755], dtype=dtype)
        ctheta = numpy.array([0.46019978478, 0.46019978478,
                              0.91743694670, 0.91743694670,
                              0.49715499774, 0.49715499774], dtype=dtype)
        phi = numpy.array([4.9261918136, 0.57595865315,
                           0.0000000000, 0.00000000000,
                           4.7005372834, 0.11129056012], dtype=dtype)
    else:
        # equinox = '(J2000)'

        psi = numpy.array([0.57477043300, 4.9368292465,
                           0.00000000000, 0.0000000000,
                           0.11142137093, 4.71279419371], dtype=dtype)
        stheta = numpy.array([0.88998808748, -0.88998808748,
                              0.39777715593, -0.39777715593,
                              0.86766622025, -0.86766622025], dtype=dtype)
        ctheta = numpy.array([0.45598377618, 0.45598377618,
                              0.91748206207, 0.91748206207,
                              0.49714719172, 0.49714719172], dtype=dtype)
        phi = numpy.array([4.9368292465, 0.57477043300,
                           0.0000000000, 0.00000000000,
                           4.71279419371, 0.11142137093], dtype=dtype)

    # zero offset
    i = select - 1
    a = (ai * D2R) - phi[i]

    b = bi * D2R
    sb = sin(b)
    cb = cos(b)
    cbsa = cb * sin(a)
    b = -stheta[i] * cbsa + ctheta[i] * sb
    w, = numpy.where(b > 1.0)
    if w.size > 0:
        b[w] = 1.0
    bo = arcsin(b) * R2D

    a = arctan2(ctheta[i] * cbsa + stheta[i] * sb, cb * cos(a))

    ao = ((a + psi[i] + fourpi) % twopi) * R2D

    return ao, bo

#
# Some clearer shortcut functions which call Euler
#
def eq2gal(ra, dec, b1950=False, dtype='f8'):
    """
    NAME
        eq2gal
    PURPOSE
        Convert from equatorial to galactic coordinates in units of degrees.
    CALLING SEQUENCE
        l,b = eq2gal(ra, dec, b1950=False, dtype='f8')
    INPUTS
        ra, dec: Equatorial coordinates.  May be Numpy arrays, sequences, or
            scalars as long as they are all the same length.  They must be
            convertible to a Numpy array with the specified datatype.
    KEYWORDS
        b1950:  If True, use b1950 coordiates.  By default j2000 are used.
        dtype:  The datatype of the output arrays.  Default is f8
    OUTPUTS
        l, b:  Galactic longitude and latitude.  The returned value is always
            a Numpy array with the specified dtype
    REVISION HISTORY
        Created Erin Sheldon, NYU, 2008-07-02
    """
    return euler(ra, dec, 1, b1950=b1950, dtype=dtype)


def gal2eq(l, b, b1950=False, dtype='f8'):
    """
    NAME
        gal2eq
    PURPOSE
        Convert from galactice to equatorial coordinates in units of degrees.
    CALLING SEQUENCE
        ra,dec = gal2eq(l, b, b1950=False, dtype='f8')
    INPUTS
        l, b: Galactic coordinates.  May be Numpy arrays, sequences, or
            scalars as long as they are all the same length.  They must be
            convertible to a Numpy array with the specified datatype.
    KEYWORDS
        b1950:  If True, use b1950 coordiates.  By default j2000 are used.
        dtype:  The datatype of the output arrays.  Default is f8
    OUTPUTS
        ra, dec:  Equatorial longitude and latitude.  The returned value is
            always a Numpy array with the specified dtype
    REVISION HISTORY
        Created Erin Sheldon, NYU, 2008-07-02
    """

    return euler(l, b, 2, b1950=b1950, dtype=dtype)

def eq2ec(ra, dec, b1950=False, dtype='f8'):
    """
    NAME
        eq2ec
    PURPOSE
        Convert from equatorial to ecliptic coordinates in units of degrees.
    CALLING SEQUENCE
        lam,beta = eq2ec(ra, dec, b1950=False, dtype='f8')
    INPUTS
        ra, dec: Equatorial coordinates.  May be Numpy arrays, sequences, or
            scalars as long as they are all the same length.  They must be
            convertible to a Numpy array with the specified datatype.
    KEYWORDS
        b1950:  If True, use b1950 coordiates.  By default j2000 are used.
        dtype:  The datatype of the output arrays.  Default is f8
    OUTPUTS
        lam, beta:  Ecliptic longitude and latitude.  The returned value is
            always a Numpy array with the specified dtype
    REVISION HISTORY
        Created Erin Sheldon, NYU, 2008-07-02
    """

    return euler(ra, dec, 3, b1950=b1950, dtype=dtype)

def ec2eq(lam, beta, b1950=False, dtype='f8'):
    """
    NAME
        ec2eq
    PURPOSE
        Convert from ecliptic to equatorial coordinates in units of degrees.
    CALLING SEQUENCE
        ra,dec = eq2gal(lam, beta, b1950=False, dtype='f8')
    INPUTS
        lam,beta: Ecliptic coordinates.  May be Numpy arrays, sequences, or
            scalars as long as they are all the same length.  They must be
            convertible to a Numpy array with the specified datatype.
    KEYWORDS
        b1950:  If True, use b1950 coordiates.  By default j2000 are used.
        dtype:  The datatype of the output arrays.  Default is f8
    OUTPUTS
        ra,dec:  Equatorial longitude and latitude.  The returned value is
            always a Numpy array with the specified dtype
    REVISION HISTORY
        Created Erin Sheldon, NYU, 2008-07-02
    """

    return euler(lam, beta, 4, b1950=b1950, dtype=dtype)

def ec2gal(lam, beta, b1950=False, dtype='f8'):
    """
    NAME
        ec2gal
    PURPOSE
        Convert from ecliptic to galactic coordinates in units of degrees.
    CALLING SEQUENCE
        l,b = eq2gal(lam, beta, b1950=False, dtype='f8')
    INPUTS
        lam, beta: Ecliptic coordinates.  May be Numpy arrays, sequences, or
            scalars as long as they are all the same length.  They must be
            convertible to a Numpy array with the specified datatype.
    KEYWORDS
        b1950:  If True, use b1950 coordiates.  By default j2000 are used.
        dtype:  The datatype of the output arrays.  Default is f8
    OUTPUTS
        l, b:  Galactic longitude and latitude.  The returned value is always
            a Numpy array with the specified dtype
    REVISION HISTORY
        Created Erin Sheldon, NYU, 2008-07-02
    """

    galact_long, galact_lat = euler(lam, beta, 5, b1950=b1950, dtype=dtype)
    sort_idxs = np.argsort(galact_long)
    return galact_long[sort_idxs], galact_lat[sort_idxs]

def gal2ec(l, b, b1950=False, dtype='f8'):
    """
    NAME
        gal2ec
    PURPOSE
        Convert from Galactic to Ecliptic coordinates in units of degrees.
    CALLING SEQUENCE
        lam,beta = eq2gal(l, b, b1950=False, dtype='f8')
    INPUTS
        l, b: Galactic coordinates.  May be Numpy arrays, sequences, or
            scalars as long as they are all the same length.  They must be
            convertible to a Numpy array with the specified datatype.
    KEYWORDS
        b1950:  If True, use b1950 coordiates.  By default j2000 are used.
        dtype:  The datatype of the output arrays.  Default is f8
    OUTPUTS
        lam,beta:  Ecliptic longitude and latitude.  The returned value is
            always a Numpy array with the specified dtype
    REVISION HISTORY
        Created Erin Sheldon, NYU, 2008-07-02
    """

    return euler(l, b, 6, b1950=b1950, dtype=dtype)


def _thetaphi2xyz(theta, phi):
    """
    theta and phi in radians relative to the SDSS node at ra=95 degrees
    """
    x = cos(theta) * cos(phi)
    y = sin(theta) * cos(phi)
    z = sin(phi)

    return x, y, z

def _xyz2thetaphi(x, y, z):
    """
    returns theta, phi in radians relative to the SDSS node at ra=95 degrees
    """
    phi = arcsin(z)
    theta = arctan2(y, x)

    return theta, phi

def eq2xyz(ra, dec, dtype='f8', units='deg'):
    """
    Convert equatorial coordinates RA and DEC to x,y,z on the unit sphere
    parameters
    ----------
    ra: scalar or array
        Right ascension. Can be an array
    dec: scalar or array
        Declination. Can be an array
    units: string, optional
        'deg' if the input is degrees, 'rad' if input
        is in radians.  Default is degrees.
    Notes:
        This follows the same convention as the STOMP package.
    """

    theta = numpy.array(ra, ndmin=1, copy=True, dtype=dtype)
    phi = numpy.array(dec, ndmin=1, copy=True, dtype=dtype)

    # in place is more efficient
    if units == 'deg':
        numpy.deg2rad(theta, theta)
        numpy.deg2rad(phi, phi)

    theta -= _sdsspar['node']

    return _thetaphi2xyz(theta, phi)

def xyz2eq(xin, yin, zin, units='deg'):
    """
    Convert x,y,z on the unit sphere to RA DEC.
    parameters
    ----------
    x,y,z:
        scalars or arrays as given by eq2xyz
    units: string, optional
        'deg' if the output is to be degrees, 'rad' if it is to be radians.
        Default is degrees.
    Notes:
        This follows the same convention as the STOMP package.
    """

    x = numpy.array(xin, ndmin=1, copy=False)
    y = numpy.array(yin, ndmin=1, copy=False)
    z = numpy.array(zin, ndmin=1, copy=False)

    theta, phi = _xyz2thetaphi(x, y, z)
    theta += _sdsspar['node']

    if units == 'deg':
        numpy.rad2deg(theta, theta)
        numpy.rad2deg(phi, phi)

    atbound(theta, 0.0, 360.0)

    # theta->ra, phi->dec
    return theta, phi

def sphdist(ra1, dec1, ra2, dec2, units=['deg', 'deg']):
    """
    Get the arc length between two points on the unit sphere
    parameters
    ----------
    ra1,dec1,ra2,dec2: scalar or array
        Coordinates of two points or sets of points.
        Must be the same length.
    units: sequence
        A sequence containing the units of the input and output.  Default
        ['deg',deg'], which means inputs and outputs are in degrees.  Units
        can be 'deg' or 'rad'
    """

    units_in, units_out = units

    # note x,y,z from eq2xyz always returns 8-byte float
    x1, y1, z1 = eq2xyz(ra1, dec1, units=units_in)
    x2, y2, z2 = eq2xyz(ra2, dec2, units=units_in)

    costheta = x1 * x2 + y1 * y2 + z1 * z2
    costheta.clip(-1.0, 1.0, out=costheta)

    theta = arccos(costheta)

    if units_out == 'deg':
        numpy.rad2deg(theta, theta)
    return theta


def gcirc(ra1deg, dec1deg, ra2deg, dec2deg, getangle=False):
    """
    This is currently very inflexible: degrees in, radians out
    """
    ra1 = numpy.array(ra1deg, dtype='f8', ndmin=1)
    dec1 = numpy.array(dec1deg, dtype='f8', ndmin=1)
    ra2 = numpy.array(ra2deg, dtype='f8', ndmin=1)
    dec2 = numpy.array(dec2deg, dtype='f8', ndmin=1)

    deg2rad(ra1, ra1)
    deg2rad(dec1, dec1)
    deg2rad(ra2, ra2)
    deg2rad(dec2, dec2)

    sindec1 = sin(dec1)
    cosdec1 = cos(dec1)

    sindec2 = sin(dec2)
    cosdec2 = cos(dec2)

    radiff = (ra2 - ra1)
    cosradiff = cos(radiff)
    cosdis = sindec1 * sindec2 + cosdec1 * cosdec2 * cosradiff

    cosdis.clip(-1.0, 1.0, out=cosdis)
    dis = arccos(cosdis)

    if getangle:
        theta = arctan2(sin(radiff),
                        (sindec1 * cosradiff - cosdec1 * sindec2 / cosdec2)) - HALFPI
        return dis, theta

    return dis

# utility functions
def atbound(longitude, minval, maxval):
    """
    Wraps the longitude until it is within the bounds. The bounds must define an entire circle (e.g
    180 to 180 or 0 to 360)
    """
    w, = numpy.where(longitude < minval)
    while w.size > 0:
        longitude[w] += 360.0
        w, = numpy.where(longitude < minval)

    w, = numpy.where(longitude > maxval)
    while w.size > 0:
        longitude[w] -= 360.0
        w, = numpy.where(longitude > maxval)


def atbound2(theta, phi):
    """
    Wraps theta and phi to be between -180 and 180 for theta, and 0 and 360 for phi.
    """
    atbound(theta, -180.0, 180.0)

    w, = numpy.where(numpy.abs(theta) > 90.0)
    if w.size > 0:
        theta[w] = 180.0 - theta[w]
        phi[w] += 180.0

    atbound(theta, -180.0, 180.0)
    atbound(phi, 0.0, 360.0)

    w, = numpy.where(numpy.abs(theta) == 90.0)
    if w.size > 0:
        phi[w] = 0.0

#
# SDSS specific conversions
#


def eq2sdss(ra_in, dec_in, dtype='f8'):
    """
    NAME:
      eq2sdss
    PURPOSE:
       Convert from ra, dec to the corrected clambda, ceta
       SDSS survey coordinate system.  It is corrected so that the
       longitude eta ranges from [-180.0, 180.0] and the latitude
       lambda ranges from [-90.0,90.0].  The standard lambda/eta
       both range from [-180.0,180.0] which doesn't make sense.
       NOTE: lambda is often referred to as longitude but this
       is incorrect since it has poles at [-90,90]
    CALLING SEQUENCE:
      from esutil import coords
      (clambda, ceta) = coords.eq2sdss(ra, dec, dtype='f8')
    INPUTS:
      ra: Equatorial latitude in degrees.
      dec: Equatorial longitude in degrees.
    OPTIONAL INPUTS:
        dtype: The data type of output.  Default is 'f8'. See
        numpy.typeDict for a list of possible types.
        dtype: The data type of output.  Default is 'f8'.
    OUTPUTS:
      clambda: Corrected Survey longitude (actually lattitude) in degrees
      ceta: Corrected Survey latitude (actually logitude) in degrees
    REVISION HISTORY:
      Written: 11-March-2006  Converted from IDL program.
    """

    # Make a copy as an array. ndmin=1 to avoid messed up scalar arrays
    ra = numpy.array(ra_in, ndmin=1, copy=True, dtype=dtype)
    dec = numpy.array(dec_in, ndmin=1, copy=True, dtype=dtype)

    if ra.size != dec.size:
        raise ValueError("RA, DEC must be same size")

    # range checking
    if ra.min() < 0.0 or ra.max() > 360.0:
        raise ValueError('RA must we within [0,360]')
    if dec.min() < -90.0 or dec.max() > 90.0:
        raise ValueError('DEC must we within [-90,90]')

    ra *= D2R
    dec *= D2R
    ra -= _sdsspar['node']

    # generate x,y,z on unit sphere, clearing memory as we go
    cdec = cos(dec)

    x = cos(ra) * cdec
    y = sin(ra) * cdec

    ra = 0
    cdec = 0 # mem

    z = numpy.sin(dec)

    dec = 0 # mem

    # generate clambda, ceta
    # do things in place to save memory

    # clambda = -arcsin( x ) (not a copy clambda=x)
    arcsin(x, x)
    clambda = x
    clambda *= -1


    arctan2(z, y, z)
    ceta = z
    ceta -= _sdsspar['etapole']

    clambda *= R2D
    ceta *= R2D

    atbound(ceta, -180.0, 180.0)

    return (clambda, ceta)

def sdss2eq(clambda_in, ceta_in, dtype='f8'):
    """
    NAME:
      sdss2eq
    PURPOSE:
       Convert corrected clambda, ceta SDSS survey coordinate system t
       equatorial coords.
    CALLING SEQUENCE:
      from esutil import coords
      (ra, dec) = coords.sdss2eq(clambda, ceta, dtype='f8')
    INPUTS:
      clambda: Corrected Survey longitude (actually lattitude) in degrees
      ceta: Corrected Survey latitude (actually logitude) in degrees
    OPTIONAL INPUTS:
        dtype: The data type of output.  Default is 'f8'. See
        numpy.typeDict for a list of possible types.
    OUTPUTS:
      ra: Equatorial latitude in degrees.
      dec: Equatorial longitude in degrees.
    REVISION HISTORY:
      Written: 11-March-2006  Converted from IDL program.
    """

    # Make a copy as an array. ndmin=1 to avoid messed up scalar arrays
    clambda = numpy.array(clambda_in, ndmin=1, copy=True, dtype=dtype)
    ceta = numpy.array(ceta_in, ndmin=1, copy=True, dtype=dtype)

    # range checking
    if clambda.min() < -90.0 or clambda.max() > 90.0:
        raise ValueError('CLAMBDA must we within [-90,90]')
    if ceta.min() < -180.0 or ceta.max() > 180.0:
        raise ValueError('CETA must we within [-180,180]')

    clambda *= D2R
    ceta *= D2R

    x = -sin(clambda)
    y = cos(ceta + _sdsspar['etapole']) * cos(clambda)
    z = sin(ceta + _sdsspar['etapole']) * cos(clambda)

    ra = arctan2(y, x) + _sdsspar['node']
    dec = arcsin(z)

    ra *= R2D
    dec *= R2D
    atbound2(dec, ra)

    return (ra, dec)


def _eq2survey(ra_in, dec_in, dtype='f8'):
    """
    NAME:
      _eq2survey
    PURPOSE:
       Convert from ra, dec to the lambda, eta
       SDSS survey coordinate system.  Note this coordinate system is
       not well defined.  Recommend you use csurvey coords.
    CALLING SEQUENCE:
      from esutil import coords
      (lambda, eta) = coords._eq2survey(ra, dec, dtype='f8')
    INPUTS:
      ra: Equatorial latitude in degrees.
      dec: Equatorial longitude in degrees.
    OPTIONAL INPUTS:
        dtype: The data type of output.  Default is 'f8'. See
        numpy.typeDict for a list of possible types.
    OUTPUTS:
      lambda: SDSS Survey longitude (actually lattitude) in degrees
      eta: SDSS Survey latitude (actually logitude) in degrees
    REVISION HISTORY:
      Written: 11-March-2006  Converted from IDL program.
    """

    # Make a copy as an array. ndmin=1 to avoid messed up scalar arrays
    ra = numpy.array(ra_in, ndmin=1, copy=True, dtype=dtype)
    dec = numpy.array(dec_in, ndmin=1, copy=True, dtype=dtype)

    if ra.size != dec.size:
        raise ValueError("RA, DEC must be same size")

    # range checking
    if ra.min() < 0.0 or ra.max() > 360.0:
        raise ValueError('RA must we within [0,360]')
    if dec.min() < -90.0 or dec.max() > 90.0:
        raise ValueError('DEC must we within [-90,90]')



    ra *= D2R
    dec *= D2R
    ra -= _sdsspar['node']

    # generate x,y,z on unit sphere, clearing memory as we go
    cdec = cos(dec)

    x = cos(ra) * cdec
    y = sin(ra) * cdec

    ra = 0
    cdec = 0 # mem

    z = sin(dec)

    dec = 0 # mem

    # generate lam, eta
    # do things in place to save memory

    # lam = -arcsin( x ) (not a copy lam=x)
    arcsin(x, x)
    lam = x
    lam *= -1


    arctan2(z, y, z)
    eta = z
    eta -= _sdsspar['etapole']

    lam *= R2D
    eta *= R2D

    atbound2(lam, eta)
    atbound(eta, -180.0, 180.0)

    w, = numpy.where(eta > (90.0 - _sdsspar['center_dec']))
    if w.size > 0:
        eta[w] -= 180.0
        lam[w] = 180.0 - lam[w]

    atbound(lam, -180.0, 180.0)

    return (lam, eta)


def dec_parse(decstring):
    """
    parse a colon separated string representing declination into
    degrees.
    parameters
    ----------
    decstring: string
        DD:MM:SS.sss the value is specified in degrees, minutes, seconds
        Only the degrees are required. Additional
        precision (minutes, seconds) are optional in the string (i.e. "12" or
        "12:34" or "12:34:56" are all valid input strings)
    Corrections by Paul Ray and Dave Smith, NRL, 2013-03-19
    """
    dec = 0.0
    sign = 1.0

    if decstring.find("-") >= 0:
        sign = -1.0
    ds = decstring.split(':')
    lds = len(ds)
    if lds >= 1:
        deg = abs(float(ds[0]))
        dec += deg
    if lds >= 2:
        minutes = float(ds[1])
        dec += minutes / 60.0
    if lds >= 3:
        sec = float(ds[2])
        dec += sec / 3600.0
    return dec * sign

def ra_parse(rastring, hours=True):
    """
    parse a colon separated string representing right ascension into
    decimal degrees.
    parameters
    ----------
    rastring: string
        "HH:MM:SS.sss" if hours is True and
        "DD:MM:SS.sss" if hours is False (indicating that
            the value is specified in degrees, minutes, seconds)
        In all cases,  only the hours (or degrees) are required. Additional
        precision (minutes, seconds) are optional in the string (i.e. "12" or
        "12:34" or "12:34:56" are all valid input strings)
    Corrections by Paul Ray and Dave Smith, NRL, 2013-03-19
    """
    ra = 0.0

    rs = rastring.split(':')
    lrs = len(rs)
    if lrs >= 1:
        ra += float(rs[0])
    if lrs >= 2:
        minutes = float(rs[1])
        ra += minutes/60.0
    if lrs >= 3:
        sec = float(rs[2])
        ra += sec/3600.0
    if hours:
        ra *= 15
    return ra

def fitsheader2dict(hdr, ext=0):
    """
    Convert a fits header object into a dict.  A dict provides more expected
    interface to the data but cannot be written back to a fits file without
    transformation.
    """

    hdict = {}
    for key in hdr:
        hdict[key.lower()] = hdr[key]

    return hdict

def shiftlon(lon_input, shift=None, wrap=True):
    """
    Name:
        shiftlon
    Calling Sequence:
        newlon = shiftlon(longitude, wrap=True, shift=0.0)
    Purpose:
        Shift the value of a longitude.  By default, the value is "wrapped" to
        be [-180,180] instead of [0,360]
        If the shift keyword is sent, then the longitude is simply shifted by
        the input value and then constrained to be again on the [0,360) range.
    Input:
        A longitude or array of longitudes on the range [0,360)
    Keywords:
        shift:
            If shift is sent, then lon-shift is returned, constrained to still
            be on [0,360).
        wrap:
            If shift is not sent, and wrap is True, wrap the range to
            [-180,180]
    """
    lon = numpy.array(lon_input, ndmin=1, copy=True, dtype='f8')

    if shift is not None:
        negshift = False
        if shift < 0:
            negshift = True

        abs_shift = abs(shift)

        # make sure in range [0,360)
        abs_shift = abs_shift % 360.0

        if negshift:
            lon += abs_shift

            w, = numpy.where(lon > 360.0)
            if w.size > 0:
                lon[w] -= 360.0
        else:
            lon -= abs_shift

            w, = numpy.where(lon < 0.0)
            if w.size > 0:
                lon[w] += 360.0

    elif wrap:
        w, = where(lon > 180)
        if w.size > 0:
            lon[w] -= 360

    return lon



def shiftra(ra, shift=None, wrap=True):
    """
    Name:
        shiftra
    Calling Sequence:
        newra = shiftra(ra, wrap=True, shift=0.0)
    Purpose:
        Shift the value of a longitude RA.  By default, the value is "wrapped"
        to be [-180,180] instead of [0,360]
        If the shift keyword is sent, then the longitude is simply shifted by
        the input value and then constrained to be again on the [0,360) range.
    Input:
        ra or any other longitude on the range [0,360)
    Keywords:
        shift:
            If shift is sent, then ra-shift is returned, constrained to still
            be on [0,360).
        wrap:
            If shift is not sent, and wrap is True, wrap the range to
            [-180,180]
    """
    return shiftlon(ra, shift=shift, wrap=wrap)



def radec2aitoff(ra, dec):
    """
    Take the ra/dec into aitoff coords
    """

    r2 = numpy.sqrt(2.0)
    f = 2. * r2 / PI


    sra = shiftra(ra)

    alpha2 = sra / 2. * D2R
    delta = dec * D2R

    cdec = cos(delta)

    denom = sqrt(1.0 + cdec * cos(alpha2))

    x = cdec * sin(alpha2) * 2. * r2 / denom
    y = sin(delta) * r2 / denom
    x = x * R2D / f
    y = y * R2D / f

    return x, y

def _check_range(rng, allowed):
    if rng is None:
        rng = allowed
    else:
        if not hasattr(rng, '__len__'):
            raise ValueError("range object does not have len() method")

        if rng[0] < allowed[0] or rng[1] > allowed[1]:
            raise ValueError(f"lon_range should be within [{allowed[0]}, {allowed[1]}]")
    return rng

def randsphere(num, ra_range=None, dec_range=None, system='eq'):
    """
    Generate random points on the sphere
    You can limit the range in ra and dec.  To generate on a spherical cap, see
    randcap()
    parameters
    ----------
    num: integer
        The number of randoms to generate
    ra_range: list, optional
        Should be within range [0,360].  Default [0,360]
    dec_range: list, optional
        Should be within range [-90,90].  Default [-90,90]
    system: string
        Default is 'eq' for the ra-dec system.  Can also be 'xyz'.
    output
    ------
        for system == 'eq' the return is a tuple
            ra,dec = randsphere(...)
        for system == 'xyz' the return is a tuple
            x,y,z = randsphere(...)
    examples
    --------
        ra,dec = randsphere(2000, ra_range=[10,35], dec_range=[-25,15])
        x,y,z = randsphere(2000, system='xyz')
    """

    ra_range = _check_range(ra_range, [0.0, 360.0])
    dec_range = _check_range(dec_range, [-90.0, 90.0])

    ra = numpy.random.random(num)
    ra *= (ra_range[1] - ra_range[0])
    if ra_range[0] > 0:
        ra += ra_range[0]

    # number [-1,1)
    cosdec_min = cos(deg2rad(90.0 + dec_range[0]))
    cosdec_max = cos(deg2rad(90.0 + dec_range[1]))
    v = numpy.random.random(num)
    v *= (cosdec_max - cosdec_min)
    v += cosdec_min

    numpy.clip(v, -1.0, 1.0, v)
    # Now this generates on [0,pi)
    dec = numpy.arccos(v)

    # convert to degrees
    rad2deg(dec, dec)
    # now in range [-90,90.0)
    dec -= 90.0

    if system == 'xyz':
        x, y, z = eq2xyz(ra, dec)
        return x, y, z

    return ra, dec

def randcap(nrand, ra, dec, rad, get_radius=False, dorot=False):
    """
    Generate random points in a sherical cap
    parameters
    ----------
    nrand:
        The number of random points
    ra,dec:
        The center of the cap in degrees.  The ra should be within [0,360) and
        dec from [-90,90]
    rad: float
        radius of the cap, same units as ra,dec
    get_radius: bool, optional
        if true, return radius of each point in radians
    dorot: bool
        If dorot is True, generate the points on the equator and rotate them to
        be centered at the desired location.  This is the default when the dec
        is within 0.1 degrees of the pole, to avoid calculation issues
    """
    # generate uniformly in r**2
    if dec >= 89.9 or dec <= -89.9:
        dorot = True

    if dorot:
        tra, tdec = 90.0, 0.0
        rand_ra, rand_dec, rand_r = randcap(nrand, 90.0, 0.0, rad, get_radius=True)
        rand_ra, rand_dec = rotate(0.0, dec - tdec, 0.0, rand_ra, rand_dec)
        rand_ra, rand_dec = rotate(ra - tra, 0.0, 0.0, rand_ra, rand_dec)
    else:
        rand_r = numpy.random.random(nrand)
        rand_r = sqrt(rand_r) * rad

        # put in degrees
        numpy.deg2rad(rand_r, rand_r)

        # generate position angle uniformly 0,2*PI
        rand_posangle = numpy.random.random(nrand) * 2 * PI

        theta = numpy.array(dec, dtype='f8', ndmin=1, copy=True)
        phi = numpy.array(ra, dtype='f8', ndmin=1, copy=True)
        theta += 90

        numpy.deg2rad(theta, theta)
        numpy.deg2rad(phi, phi)

        sintheta = sin(theta)
        costheta = cos(theta)

        sinr = sin(rand_r)
        cosr = cos(rand_r)

        cospsi = cos(rand_posangle)
        costheta2 = costheta * cosr + sintheta * sinr * cospsi

        numpy.clip(costheta2, -1, 1, costheta2)

        # gives [0,pi)
        theta2 = arccos(costheta2)
        sintheta2 = sin(theta2)

        cosDphi = (cosr - costheta * costheta2) / (sintheta * sintheta2)

        numpy.clip(cosDphi, -1, 1, cosDphi)
        Dphi = arccos(cosDphi)

        # note fancy usage of where
        phi2 = numpy.where(rand_posangle > PI, phi + Dphi, phi - Dphi)

        numpy.rad2deg(phi2, phi2)
        numpy.rad2deg(theta2, theta2)
        rand_ra = phi2
        rand_dec = theta2 - 90.0

        atbound(rand_ra, 0.0, 360.0)

    if get_radius:
        #numpy.rad2deg(rand_r, rand_r)
        return rand_ra, rand_dec, rand_r

    return rand_ra, rand_dec

def randcap_brute(nrand, ra, dec, rad, get_radius=False):
    """
    Generate random points in a sherical cap using brute
    force rejection sampling. This is extremely
    slow and is used for testing purposes only.
    parameters
    ----------
    nrand: int
        The number of random points
    ra,dec: float
        The center of the cap in degrees.  The ra should be within [0,360) and
        dec from [-90,90]
    rad: float
        radius of the cap, same units as ra,dec
    get_radius: bool, optional
        if true, return radius of each point in radians
    """

    ora = zeros(nrand)
    odec = zeros(nrand)
    orad = zeros(nrand)

    ngood = 0
    nleft = nrand

    while ngood < nrand:
        tra, tdec = randsphere(nleft)
        d = sphdist(ra, dec, tra, tdec)
        w, = where(d <= rad)
        if w.size > 0:
            ora[ngood:ngood + w.size] = tra[w]
            odec[ngood:ngood + w.size] = tdec[w]
            orad[ngood:ngood + w.size] = d[w]

            ngood += w.size
            nleft -= w.size

    if get_radius:
        return ora, odec, orad
    return ora, odec

def rotate(phi, theta, psi, ra, dec):
    """
    rotation the given positions on the sphere
    The convention is the usual zxz
    Parameters
    ----------
    phi, theta, psi: numbers
        The euler angles in zxz convention
    ra, dec: numbers or arrays
        positions to be rotated
    """

    if hasattr(ra, '__len__'):
        is_scalar = False
    else:
        is_scalar = True

    ra = numpy.array(ra, ndmin=1, copy=False)
    dec = numpy.array(dec, ndmin=1, copy=False)
    if ra.size != dec.size:
        raise ValueError('ra[%d] has different size than dec[%d]' % (ra.size, dec.size))

    twopi = 2.0 * PI
    fourpi = 4.0 * PI

    # use negative; rotating the points is like rotating
    # the coord system in the opposite direction
    phi = deg2rad(-phi)
    theta = deg2rad(-theta)
    psi = deg2rad(-psi)

    sintheta = sin(theta)
    costheta = cos(theta)


    a = deg2rad(ra) - phi
    b = deg2rad(dec)

    sb = sin(b)
    cb = cos(b)
    cbsa = cb * sin(a)

    b = -sintheta * cbsa + costheta * sb

    w, = numpy.where(b > 1.0)
    if w.size > 0:
        b[w] = 1.0

    dec_out = arcsin(b)

    a = arctan2(costheta * cbsa + sintheta * sb, cb * cos(a))
    ra_out = (a + psi + fourpi) % twopi

    rad2deg(ra_out, out=ra_out)
    rad2deg(dec_out, out=dec_out)

    if is_scalar:
        ra_out = ra_out[0]
        dec_out = dec_out[0]

    return ra_out, dec_out


def rect_area(lon_min, lon_max, lat_min, lat_max):
    """
    Calculate the area of a rectangle on the sphere.
    parameters
    ----------
    lon_min, lon_max, lat_min, lat_max:
        Definition of the rectangle, in degrees
    """
    smax = sin(deg2rad(lat_max))
    smin = sin(deg2rad(lat_min))
    area = (smax - smin) * (lon_max - lon_min)
    return numpy.abs(area) * R2D