J/ApJS/217/27       Morphologies of z<0.01 SDSS-DR7 galaxies       (Ann+, 2015)
================================================================================
A catalog of visually classified galaxies in the local (z ~ 0.01) universe.
    Ann H.B., Seo M., Ha D.K.
   <Astrophys. J. Suppl. Ser., 217, 27 (2015)>
   =2015ApJS..217...27A
================================================================================
ADC_Keywords: Galaxies, nearby ; Morphology ; Photometry, SDSS
Keywords: catalogs; galaxies: general; galaxies: statistics; galaxies: structure

Abstract:
    The morphological types of 5836 galaxies were classified by a visual
    inspection of color images using the Sloan Digital Sky Survey Data
    Release 7 to produce a morphology catalog of a representative sample
    of local galaxies with z<0.01. The sample galaxies are almost complete
    for galaxies brighter than r_pet_=17.77. Our classification system is
    basically the same as that of the Third Reference Catalog of Bright
    Galaxies with some simplifications for giant galaxies. On the other
    hand, we distinguish the fine features of dwarf elliptical (dE)-like
    galaxies to classify five subtypes: dE, blue-cored dwarf ellipticals,
    dwarf spheroidals (dSph), blue dwarf ellipticals (dE_blue_), and dwarf
    lenticulars (dS0). In addition, we note the presence of nucleation in
    dE, dSph, and dS0. Elliptical galaxies and lenticular galaxies
    contribute only ~1.5 and 4.9% of local galaxies, respectively, whereas
    spirals and irregulars contribute ~32.1 and ~42.8%, respectively. The
    dE_blue_ galaxies, which are a recently discovered population of
    galaxies, contribute a significant fraction of dwarf galaxies. There
    seem to be structural differences between dSph and dE galaxies. The
    dSph galaxies are fainter and bluer with a shallower surface
    brightness gradient than dE galaxies. They also have a lower fraction
    of galaxies with small axis ratios (b/a<~0.4) than dE galaxies. The
    mean projected distance to the nearest neighbor galaxy is ~260 kpc.
    About 1% of local galaxies have no neighbors with comparable
    luminosity within a projected distance of 2Mpc.

Description:
    This paper presents a catalog of the morphological types of galaxies
    whose redshifts are less than z=0.01. The morphological types are
    determined by a visual inspection of the color images provided by SDSS
    DR7 (II/294). The majority of galaxies in the present sample come from
    the KIAS-VAGC (Choi et al. 2010JKAS...43..191C) which is based on the
    spectroscopic target galaxies of the SDSS DR7 complemented by the
    bright galaxies with known redshifts from various catalogs.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records  Explanations
--------------------------------------------------------------------------------
ReadMe         80        .  This file
table1.dat    121     5840  Catalog of the morphological types of local galaxies
--------------------------------------------------------------------------------

See also:
 II/294  : The SDSS Photometric Catalog, Release 7 (Adelman-McCarthy+, 2009)
 VII/250 : The 2dF Galaxy Redshift Survey (2dFGRS) (2dFGRS Team, 1998-2003)
 VII/98  : Catalogue of 2810 nearby galaxies (Kraan-Korteweg, 1986)
 J/MNRAS/446/3749 : SDSS nearby galaxies morphologies (Yoshino+, 2015)
 J/MNRAS/445/630  : Early-type galaxies in Ursa Major cluster (Pak+, 2014)
 J/MNRAS/435/2835 : Morphological types from Galaxy Zoo 2 (Willett+, 2013)
 J/AJ/146/151     : Morphology catalog of nearby galaxies from SDSS (Oh+, 2013)
 J/AJ/145/101     : Updated nearby galaxy catalog (Karachentsev+, 2013)
 J/MNRAS/429/2264 : Anatomy of Ursa Majoris (Karachentsev+, 2013)
 J/AJ/144/4       : Properties of dwarf galaxies in the LG (McConnachie+, 2012)
 J/ApJ/743/123    : Motions of galaxies in Coma I cloud (Karachentsev+, 2011)
 J/MNRAS/412/2498 : Galaxy groups in the local universe (Makarov+, 2011)
 J/MNRAS/412/727  : SDSS red gal. automated morph. classification (Cheng+, 2011)
 J/ApJ/722/L120   : Central surface brightness of SDSS galaxies (Fathi, 2010)
 J/ApJ/711/361    : Local Group dE galaxies. II. (Geha+, 2010)
 J/ApJS/186/427   : Detailed morphology of SDSS galaxies (Nair+, 2010)
 J/MNRAS/398/1129 : Central galaxies in groups and clusters (Guo+, 2009)
 J/ApJS/182/216   : Surface photometry of Virgo ellipticals (Kormendy+, 2009)
 J/ApJ/660/1186   : VCC galaxies classification with SDSS (Lisker+, 2007)
 J/AJ/128/163     : Galaxy morphological classification (Lotz+, 2004)
 J/ApJS/147/29    : BRHalpha data of blue compact dwarf gal. (Gil De Paz+, 2003)
 J/ApJS/147/1     : Classification of nearby galaxies (Conselice+, 2003)
 J/ApJ/588/218    : i*g* photometry of SDSS EDR galaxies (Abraham+, 2003)
 J/AJ/90/1681     : The Virgo Cluster Catalog (VCC) (Binggeli+, 1985)
 http://www.sdss.org/ : SDSS home page

Byte-by-byte Description of file: table1.dat
--------------------------------------------------------------------------------
   Bytes Format Units  Label  Explanations
--------------------------------------------------------------------------------
   1- 18  A18   ---    ID     SDSS identifier (ObjID or SDSSJHHMMSS.ss+DDM)
  20- 29  F10.6 deg    RAdeg  Right Ascension in decimal degrees (J2000)
  31- 40  F10.6 deg    DEdeg  Declination in decimal degrees (J2000)
  42- 50  F9.6  ---    z      [-0.003/0.01] Heliocentric spectroscopic redshift
  52- 58  A7    ---    Morph  Morphology from this study
  60- 62  I3    ---    T      [-11/13] Numeric morphology index (1)
  64- 68  F5.2  Mpc    Dist   [0/42.2] Distance (2)
  70- 75  F6.2  mag    rMag   [-22.3/0.5] Absolute SDSS DR7 r band magnitude (3)
  77- 81  F5.2  mag    u-r    [-2.3/15.5]?=-9 Extinction corrected model (u-r)
                               color index
  83- 89  F7.1  pix    amaj   ?=-9 Isophotal semi-major axis measured at
                                   {mu}_r_=25mag/arcsec^2^
  91- 95  F5.2  ---    b/a    [0/1]?=-9 Isophotal axis ratio
  97-121  A25   ---    Name   Primary source identifier (4)
--------------------------------------------------------------------------------
Note (1): For practical purposes, the numerical type T borrowed from RC3 is
          used with some simplification (see section 3.1):
  * ellipticals: T=-5
  * lenticulars: T=-3
  * spirals have their Hubble stages from T=0 (S0/a) to T=9 (Sm)
  * Irregulars: T=10 (Im), T=11 (dI), and
    T=12 (irregulars containing small BCD-like components)
  * Blue compact dwarf (BCD): T=13
  * dwarf ellipticals: T=-6 (dE), T=-7 (dEbc), T=-8 (dSph), T=-9 (dEblue),
    T=-10 (dS0, dS0p), and T=-11 for transition type dwarf (dEs/dI dI/dEs),
Note (2): Derived from the radial velocity relative to the Local Group.
Note (3): We use the model magnitude from SDSS DR7 corrected for the galactic
          extinction (Schlegel et al. 1998ApJ...500..525S). The SDSS model
          magnitude is a magnitude derived by best fitting model (exponential
          profile or de Vaucouleurs profile) to the observed images.
Note (4): The first entry in NED is taken. In cases of no NED name, we make a
          name using the coordinates of the galaxy similar to the SDSS
          object name.
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

================================================================================
(End)                 Greg Schwarz [AAS], Emmanuelle Perret [CDS]    11-May-2015
