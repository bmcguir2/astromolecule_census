#!/usr/bin/env python

#############################################################
#						Revision History					#
#############################################################

# generates and writes out a ton of things

# 1.0 - Project start (5/12/2018)
# 1.1 - Cleaned up code
# 1.2 - Added sulfur compounds to DU calculation.
# 1.3 - Update to O2 detection paper.  Now references earlier work by Larsson et al. 2007.
# 1.4 - Adds rotational constants, dipole moments, and asymmetry parameters. (1/17/2019)
# 2.0 - Updates for 2019 Census

#############################################################
#							Preamble						#
#############################################################

'''

I wrote this code primarily to help me easily generate ascii files for making the plots 
found in the main text (2019 Census of Interstellar [...] Molecules).  I have
attempted to comment the code to some degree to make it more readable and potentially
useful for others.  The functions at the end of the code produce either the figures in
the manuscript a python plotting program, or the data needed to make them in another 
program.  

A few critical notes for using:

1) I wrote this code intending for it to be loaded and used in an interactive python
environment such as IPython or Jupyter Notebooks.  I don't know how it will behave 
outside of these, but it should be amenable to scripting.

2) It is in Python 3; it *should* work in Python 2.7, I think, but no promises.

3) Molecules and Sources are their own classes.  Each molecule or class is given a 
variable name that is intuitive to me (largely removing spaces and punctuation), but
may not be intuitive to everyone.  You can just search the document in plain text to find
the variable name you need.

4) I have written a utility function for providing a quick look at bulk of the data for
an individual molecule or source.

	>> summary(y)
	
where y is a molecule tag or source tag.  If a source is a line of sight
source, add LOS to the end of the standard source tag.  For example, Sgr B2 is SgrB2 and
the line of sight to Sgr B2 is SgrB2LOS.

5) There is also a utility function for writing out summaries to ascii text files:

	>> output_summary(y,filename=None)
	
This will write the output of summary(y) to a text file.  For a single molecule, the
default filename is formula.txt, but this can be overridden when the command is issued.
y can also be a list of molecules, including full_list, which was used to generate the 
ascii version of the database uploaded as supplementary information.

6) Many additional properties of molecules are under development, and should not be
considered exhaustive/comprehensive.  For example, detected isotopologues are being added,
but the absence of an isotopologue from the list should not be considered a non-detection.
Similarly, isomers are being compiled.  In other words, there is a lot of infrastructure 
in place for future development.

7) Notification of any discovered typographical, bookkeeping, or content errors to 
bmcguire@nrao.edu is greatly appreciated.  If you write your own function that works with
the code without modification (with the exception of an additional import), and would like
to contribute it for the next update to the Census, that would also be most welcome.

8) I am not a programmer, coder, or developer.  I am completely confident that this code
is shockingly suboptimal, bloated, and un-pythonic.  But it does what I want it to do =).

9) Dipole moments were obtained from the CDMS and JPL database entries for these species,
except where noted.  As I work through updates, I am including references to the original 
work these databases reference, although in many cases these dipole moments were 
calculated by the database curators themselves.

'''

import os, sys, argparse, math
import numpy as np
from operator import itemgetter
from math import ceil
from datetime import date
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import periodictable as pt
import matplotlib.patches as patches
from colour import Color
import seaborn as sns
from scipy.stats import gaussian_kde as gkde
from scipy.interpolate import make_interp_spline, BSpline

matplotlib.rc('text', usetex = True)
matplotlib.rc('text.latex',preamble=r'\usepackage{cmbright}\usepackage[version=4]{mhchem}')

#Python version check

if sys.version_info.major != 3:

	print("Warning: This code is written in Python 3.  It may not function properly in Python 2.7.")
	
version = 2.0

#############################################################
#					    Telescope Class  					#
#############################################################

class Telescope(object):

	def __init__(self,name,shortname,type=None,wavelength=None,latitude=None,longitude=None,diameter=None,built=None,decommissioned=None,notes=None):
	
		self.name = name
		self.shortname = shortname
		self.type = type
		self.wavelength = wavelength
		self.latitude = latitude
		self.longitude = longitude
		self.diameter = diameter
		self.built = built
		self.decommissioned = decommissioned
		self.notes = notes
		self.ndetects = None
		self.mol_list = None
	
		return
		
	def update_stats(self,mol_list):
	
		'''
		Takes a list of molecules to loop over and add to the list of molecules detected with this telescope, as well as the detects counter.
		'''
		
		my_mols = []
	
		for mol in mol_list:
		
			if self in mol.telescopes:
			
				my_mols.append(mol)
				
		self.mol_list = my_mols
				
		self.ndetects = len(my_mols)
					
		return	


#############################################################
#					    Telescopes	    					#
#############################################################

#An exceptional resource for these is 'Observatories and Telescopes of Modern Times' by David Leverington

GBT = Telescope('Green Bank Telescope','GBT',type='Single Dish',wavelength=['cm','mm'],latitude=38.433056,longitude=-79.839722,diameter=100,built=2004)
IRAM30 = Telescope('IRAM 30-m','IRAM',type='Single Dish',wavelength=['mm','sub-mm'],latitude=37.066161,longitude=-3.392719,diameter=30,built=1984)
Spitzer = Telescope('Spizter','Spitzer',type='Space',wavelength=['IR'],diameter=0.85,built=2003,decommissioned=2020)
Hubble = Telescope('Hubble Space Telescope','Hubble',type='Space',wavelength=['IR','Vis', 'UV'],diameter=2.4,built=1990)
ALMA = Telescope('Atacama Large Millimeter/sub-millimeter Array','ALMA',type='Interferometer',wavelength=['mm','sub-mm'],latitude=-23.0193,longitude=-67.7532,built=2011)
ISO = Telescope('Infrared Space Observatory','ISO',type='Space',wavelength=['IR'],diameter=0.6,built=1995,decommissioned=1998)
NRAO140 = Telescope('NRAO 140-ft','NRAO 140-ft',type='Single Dish',wavelength=['cm'],latitude=38.433056,longitude=-79.839722,diameter=43,built=1965,decommissioned=2008,notes='Technically started operations again in 2014, but not for PI science.')
Algonquin46 = Telescope('Algonquin 46-m Telescope','Algonquin 46-m',type='Single Dish',wavelength=['cm','mm'],latitude=45.955503,longitude=-78.073042,diameter=46,built=1966,decommissioned=1987,notes='Technically in operation much longer, but seems to have ceased PI science in 1987.')
NRAOARO12 = Telescope('NRAO/ARO 12-m Telescope','NRAO/ARO 12-m',type='Single Dish',wavelength=['mm'],latitude=31.9533,longitude=-111.615,diameter=12,built=1984,notes='Originally the NRAO 36-ft (11 m) telescope until 1984, when the dish was replaced (these are counted as separate facilities).  The observatory was handed over to the ARO in 2000 and renamed.  In 2013, the antenna was replaced with a 12-m ALMA prototype antenna.')
NRAO36 = Telescope('NRAO 36-ft Telescope','NRAO 36-ft',type='Single Dish',wavelength=['mm'],latitude=31.9533,longitude=-111.615,diameter=11,built=1967,decommissioned=1984,notes='Became the NRAO/ARO 12-m telescope in 1984.')
Nobeyama45 = Telescope('Nobeyama 45-m Telescope','Nobeyama',type='Single Dish',wavelength=['cm','mm'],latitude=35.9417,longitude=138.4758,diameter=45,built=1982)
Effelsberg100 = Telescope('Effelsberg 100-m Telescope','Effelsberg',type='Single Dish',wavelength=['cm'],latitude=50.5247,longitude=-6.8828,diameter=100,built=1972)
Haystack37 = Telescope('Haystack 37-m Telescope','Haystack',type='Single Dish',wavelength=['cm','mm'],latitude=42.6233,longitude=-71.4882,diameter=37,built=1964)
PdBI = Telescope('Plateu de Bure Interferometer','PdBI',type='Interferometer',wavelength=['mm'],latitude=44.63389,longitude=5.90792,built=1988,decommissioned=2016)
NOEMA = Telescope('Northern Extended Millimeter Array','NOEMA',type='Interferometer',wavelength=['mm'],latitude=44.63389,longitude=5.90792,built=2016)
BIMA = Telescope('Berkeley-Illinois-Maryland Array','BIMA',type='Interferometer',wavelength=['mm'],latitude=40.8178,longitude=-121.473,built=1986,decommissioned=2005,notes='Became part of CARMA.')
OVRO = Telescope('Caltech Owens Valley Radio Observatory Millimeter Array','OVRO',type='Interferometer',wavelength=['mm'],latitude=37.2339,longitude=-118.282,built=1984,decommissioned=2005,notes='Became part of CARMA.')
Yebes40 = Telescope('Yebes RT40-m Telescope','Yebes',type='Single Dish',wavelength=['cm','mm'],latitude=40525208,longitude=-3.088725,built=2007)
NRL85 = Telescope('Maryland Point Observatory Naval Research Lab 85-foot Telescope','NRL 85-ft',type='Single Dish',wavelength=['cm'],diameter=26,latitude=38.3741667,longitude=-77.230833,built=1965,decommissioned=1994,notes='Primarily used for VLBI for much of its later years.')
ATCA = Telescope('Australia Telescope Compact Array','ATCA',type='Interferometer',wavelength=['cm'],latitude=-30.312778,longitude=149.550278,built=1988)
Parkes64 = Telescope('Parkes 64-m Telescope','Parkes',type='Single Dish',wavelength=['cm'],diameter=64,latitude=-32.99778,longitude=148.26292,built=1961)
SMT10 = Telescope('ARO 10-m Submillimeter Telescope','SMT',type='Single Dish',wavelength=['mm','sub-mm'],diameter=10,latitude=32.701658,longitude=-109.871391,built=1993)
SEST15 = Telescope('Swedish-ESO 15-m Submillimetre Telescope','SEST',type='Single Dish',wavelength=['mm','sub-mm'],diameter=15,latitude=-29.26,longitude=-70.73,built=1987,decommissioned=2003)
Goldstone70 = Telescope('Goldstone 72-m (DSS-14; "Mars")','Goldstone',type='Single Dish',wavelength=['cm'],diameter=70,latitude=35.426667,longitude=-116.89,built=1966,notes='Originally a 64-m dish; become 70-m in 1988. Conceivably still PI Science Capable?')
Mitaka6 = Telescope('Tokyo Astronomical Observatory Mitaka 6-m','Mitaka 6-m',type='Single Dish',wavelength=['mm'],diameter=6,latitude=35.675217,longitude=139.538083,built=1970,decommissioned=2018,notes='Moved around quite a bit within Japan until returning (and retiring) in 2018.')
McMath = Telescope('McMath-Pierce Solar Telescope','McMath Solar Telescope',type='Optical',wavelength=['IR','Vis','UV'],diameter=1.6,latitude=31.9584,longitude=-111.595,built=1962)
Bell7m = Telescope('AT&T Bell Laboratories 7-m Telescope', 'Bell 7-m', type='Single Dish',wavelength=['cm'],diameter=7,built=1976,decommissioned=1992)
IRTF = Telescope('NASA Infrared Telescope Facility', 'IRTF',type='Optical',wavelength=['IR'],diameter=3,latitude=19.8263,longitude=-155.473,built=1974)
KPNO4m = Telescope('Mayall 4-m Telescope', 'KPNO 4-m',type='Optical',wavelength=['IR'],diameter=4,latitude=31.9583,longitude=-111.5967,built=1973)
Onsala20m = Telescope('Onsala 20-m Telescope', 'Onsala 20-m',type='Single Dish',wavelength=['cm','mm'],diameter=20,latitude=57.393056,longitude=11.917778,built=1976)
FCRAO14m = Telescope('Five College Radio Observatory 14-m Telescope', 'FCRAO 14-m', type='Single Dish',wavelength=['cm','mm'],latitude=42.391925,longitude=-72.344097,built=1976,decommissioned=2005)
APEX = Telescope('Atacama Pathfinder Experiment', 'APEX', type='Single Dish',wavelength=['mm','sub-mm'],diameter=12,latitude=-23.0058,longitude=-67.7592,built=2005)
CSO = Telescope('Caltech Submillimeter Observatory', 'CSO', type='Single Dish', wavelength=['mm','sub-mm'],diameter=10.4,latitude=19.8225,longitude=-155.70694,built=1986,decommissioned=2015)
MWO4m = Telescope('University of Texas Millimeter Wave Observatory 4.9-m Telescope','MWO 4.9-m',type='Single Dish',wavelength=['mm'],latitude=30.3866,longitude=-97.7269,built=1971,decommissioned=1988)
HatCreek = Telescope('Hat Creek Station 20-ft Telescope','Hat Creek 20-ft',type='Single Dish',wavelength=['cm','mm'],latitude=40.8178,longitude=-121.473,built=1965,decommissioned=1983,notes='Best build date found was "mid 1960s", so 1965 is an estimate.  It appears to have either been subsummed into BIMA or decommissioned when BIMA came online.  The decommissioning date is an estimate.')
SMA = Telescope('Submillimeter Array','SMA',type='Interferometer',wavelength=['mm'],latitude=19.8225,longitude=-155.70694,built=2003)
Herschel = Telescope('Herschel Space Telescope','Herschel',type='Space',wavelength=['sub-mm','IR'],built=2009,decommissioned=2013)
UKIRT = Telescope('United Kingdom Infrared Telescope','UKIRT',type='Optical',wavelength=['IR'],diameter=3.8,latitude=19.8225,longitude=-155.70694,built=1979)
SOFIA = Telescope('Stratospheric Observatory for Infrared Astronomy','SOFIA',type='Airborne',wavelength=['sub-mm','IR'],diameter=2.5,built=2010)
Odin = Telescope('Odin','Odin',type='Space',wavelength=['sub-mm'],diameter=1.1,built=2001)
FUSE = Telescope('Far Ultraviolet Spectroscopic Explorer','FUSE',type='Space',wavelength=['UV'],built=1999,decommissioned=2007)
Kuiper = Telescope('Kuiper Airborne Observatory','KAO',type='Airborne',wavelength=['sub-mm','IR'],built=1974,decommissioned=1995)
MtHopkins = Telescope('Tillinghast 60 inch','Mt. Hopkins 60-in',type='Optical',diameter=1.5,wavelength=['IR'],built=1969,latitude=31.6811,longitude=-110.878)
Aerobee = Telescope('Aerobee-150 Rocket','Aerobee-150 Rocket',type='Airborne',wavelength=['UV'],built=1970,decommissioned=1970,notes='They literally put a spectrometer on a rocket and shot it into the sky, then used the same spectrometer to measure H2 in the laboratory.')
Millstone = Telescope('Lincoln Laboratory Millstone Hill Observatory 84-ft', 'Millstone Hill 84-ft',type='Single Dish',wavelength=['cm'],built=1956,decommissioned=1978,diameter=26,latitude=42.6233,longitude=-71.4882,notes='Originally the Ballistic Missile Early Warning System radar antenna.  Decommissioning date is a best guess based on the installation of larger telescopes to the site at that time.')
MtWilson = Telescope('Mount Wilson 100-in','Mt. Wilson',type='Optical',wavelength=['UV','VIS'],diameter=2.54,built=1917,decommissioned=1989,notes='Now known as the Hooker Telescope.')
IRAS = Telescope('Infrared Astronomical Satellite','IRAS',type='Space',wavelength=['IR'],diameter=0.60,built=1983,decommissioned=1983)


scopes_list = [
	GBT,
	IRAM30,
	Spitzer,
	Hubble,
	ALMA,
	ISO,
	NRAO140,
	Algonquin46,
	NRAOARO12,
	NRAO36,
	Nobeyama45,
	Effelsberg100,
	Haystack37,
	PdBI,
	NOEMA,
	BIMA,
	OVRO,
	Yebes40,
	NRL85,
	ATCA,
	Parkes64,
	SMT10,
	SEST15,
	Goldstone70,
	Mitaka6,
	McMath,
	Bell7m,
	IRTF,
	KPNO4m,
	Onsala20m,
	FCRAO14m,
	APEX,
	CSO,
	MWO4m,
	HatCreek,
	SMA,
	Herschel,
	UKIRT,
	SOFIA,
	Odin,
	FUSE,
	Kuiper,
	MtHopkins,
	Aerobee,
	Millstone,
	MtWilson,
	IRAS,
	]	

#############################################################
#						Source Class  						#
#############################################################

class Source(object):

	def __init__(self,name,type=None,ra=None,dec=None,detects=0,mols=None,simbad_url=None):
	
		self.name = name
		self.type = type
		self.ra = ra
		self.dec = dec
		self.detects = detects
		self.mols = mols
		self.simbad_url = simbad_url
					
		return	
		
	def update_stats(self,list):
	
		for x in list:
		
			if self in x.sources:
			
				if self.mols is None:
				
					self.mols = [x]
					
				else:
				
					self.mols.append(x)
					
				self.detects += 1
				
		return
		

#############################################################
#							Sources  						#
#############################################################

'''
Source coordinates are generalized for simplicity, and individual detections may have 
been made (read: absolutely have been made) toward various pointing positions within 
these sources.  

Similarly, source types are highly generalized, to allow for some 
generalized, aggregate analysis.  

Finally, some detections have been made along the line of sight to these sources.  In 
these cases,  the source is catagorized as a 'LOS Cloud', regardless of the actual source
type, as it was only used as an absorbing background.

Diffuse, translucent, and dense clouds are all classified this way for simplicity.

Please check the individual papers for details.
'''

AFGL890LOS = Source("AFGL 890 LOS",type='LOS Cloud',ra='06:10:48.0',dec='-06:12:00',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=AFGL+890')
AFGL961LOS = Source("AFGL 961 LOS",type='LOS Cloud',ra='06:34:37.741',dec='+04:12:44.20',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=AFGL961')
AFGL989LOS = Source("AFGL 989 LOS",type='LOS Cloud',ra='06:41:10.06',dec='+09:29:35.8',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=AFGL989')
B1b = Source("B1-b",type='Dark Cloud',ra='03:33:20.8',dec='+31:07:34',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%5BHKM99%5D+B1-b')
CRL2688 = Source("CRL 2688",type='Carbon Star',ra='21:02:18.27',dec='+36:41:37.0',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=CRL+2688')
CRL618 = Source("CRL 618",type='Carbon Star',ra='04:42:53.62',dec='+36:06:53.40',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=CRL+618')
CasALOS = Source("Cas A LOS",type='LOS Cloud',ra='23:23:24.00',dec='+58:48:54.0',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Cas+A')
CrabNebula = Source("Crab Nebula",type='SNR',ra='05:34:31.94',dec='+22:00:52.2',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Crab+Nebula')
CygnusOB212LOS = Source("Cygnus OB2 - 12",type='LOS Cloud',ra='20:32:40.96',dec='+41:04:13.2',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%4011680932&Name=Schulte%2012')
DR21 = Source("DR 21",type='SFR',ra='20:39:01.6',dec='+42:19:38',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=DR+21')
DR21LOS = Source("DR 21 LOS", type='LOS Cloud',ra='20:39:01.6',dec='+42:19:38',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=DR+21')
DR21OH = Source("DR 21(OH)",type='SFR',ra='20:39:01.01',dec='+42:22:50.22',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=DR+21%28OH%29')
G0693 = Source("G+0.693-0.027", type='Shock',ra='17:47:21.86',dec='-28:22:43.00')
G327306LOS = Source("G327.3-0.6 LOS",type='LOS Cloud',ra='15:53:05.0',dec='-54:35:24',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=G327.3-0.6')
GL2136LOS = Source("GL2136 LOS",type='LOS Cloud',ra='18:27:18.43',dec='-25:04:02.84',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402520798&Name=GJ%20%202136%20B')
GalacticCenter = Source("Galactic Center",type='SFR')
HD124314LOS = Source("HD 124314 LOS",type='LOS Cloud',ra='14:15:01.61',dec='-61:42:24.38',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=HD+124314')
HD27778LOS = Source("HD 27778 LOS",type='LOS Cloud',ra='04:23:59.78',dec='24:18:03.53',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=HD+27778')
HorseheadPDR = Source("Horsehead PDR",type='PDR',ra='05:40:53.936',dec='-02:28:00',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%40828287&Name=NAME%20Horsehead%20Nebula')
IC443G = Source("IC 443G",type='SNR',ra='06:16:43.4',dec='+22:32:24',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=IC+443G')
IRAS16293 = Source("IRAS 16293",type='Protostar',ra='16:32:22.56',dec='-24:28:31.8',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=IRAS+16293-2422')
IRC10216 = Source("IRC+10216",type='Carbon Star',ra='09:47:57.406',dec='+13:16:43.56',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=IRC%2B10216')
K350 = Source("K3-50",type='HII',ra='20:04:45.59',dec='33:32:42.0',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402904320&Name=NAME%20K%203-50A')
L134 = Source("L134",type='Dark Cloud',ra='15:53:36.3',dec='-04:35:26.0',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402622026&Name=LDN%20%20134')
L1527 = Source("L1527",type='Dark Cloud',ra='04:39:53.0',dec='+25:45:00',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=L1527')
L1544 = Source("L1544",type='Dark Cloud',ra='05:04:16.6',dec='+25:10:48',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=L1544')
L183 = Source("L183",type='Dark Cloud',ra='15:54:12.2',dec='-02:49:42',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=L183')
L483 = Source("L483",type='Dark Cloud',ra='18:17:35.0',dec='-04:39:48',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=L483')
LOSCloud = Source("LOS Cloud",type='LOS Cloud')
Lupus1A = Source("Lupus-1A",type='Dark Cloud',ra='15:42:52.4',dec='-34:07:53.5',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%405549180&Name=NAME%20Lupus-1A')
M17LOS = Source("M17 LOS",type='LOS Cloud',ra='18:20:47',dec='-16:10:18',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=M17')
M17SW = Source("M17SW",type='PDR',ra='18:20:23.1',dec='-16:11:43',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=M17+SW')
M3LOS = Source("M3 LOS",type='LOS Cloud',ra='13:42:11.62',dec='+28:22:38.2',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=M3')
NGC2024 = Source("NGC 2024",type='PDR',ra='05:41:43',dec='-01:50:30',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+2024')
NGC2024LOS = Source("NGC 2024 LOS",type='LOS Cloud',ra='05:41:43',dec='-01:50:30',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+2024')
NGC2264 = Source("NGC 2264",type='YSO',ra='06:40:58',dec='+09:53:42',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+2264')
NGC6334 = Source("NGC 6334",type='SFR',ra='17:20:53.3',dec='-35:46:59',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402361705&Name=NAME%20NGC%206334-I')
NGC6334LOS = Source("NGC 6334 LOS",type='LOS Cloud',ra='17:20:53.3',dec='-35:46:59',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402361705&Name=NAME%20NGC%206334-I')
NGC7023 = Source("NGC 7023",type='PDR',ra='21:01:36.9',dec='+68:09:48',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7023')
NGC7027 = Source("NGC 7027",type='PN',ra='21:07:01.8',dec='+42:14:10',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7023')
NGC7538 = Source("NGC 7538",type='YSO',ra='23:13:37.2',dec='61:30:00',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7538')
NGC7538LOS = Source("NGC 7538 LOS",type='LOS Cloud',ra='23:13:37.2',dec='61:30:00',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+7538')
Orion = Source("Orion",type='SFR',ra='05:35:14.16',dec='-05:22:21.5',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Orion+KL')
OrionBar = Source("Orion Bar",type='PDR',ra='05:35:22.30',dec='-05:24:33.0',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Orion+Bar')
rhoOphA = Source("rho Ophiuchi A",type='SFR',ra='16:26:27.20',dec='-24:24:04',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402488602&Name=NAME%20rho%20Oph%20A%20SM%201')
SgrA = Source("Sgr A",type='Sgr A',ra='17:45:40.0',dec='-29:00:28.2',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+A')
SgrALOS = Source("Sgr A LOS",type='LOS Cloud',ra='17:45:40.0',dec='-29:00:28.2',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+A')
SgrB2 = Source("Sgr B2",type='SFR',ra='17:47:20.4',dec='-28:23:07',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+B2')
SgrB2LOS = Source("Sgr B2 LOS",type='LOS Cloud',ra='17:47:20.4',dec='-28:23:07',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sgr+B2')
TC1 = Source("TC 1",type='PN',ra='17:45:35.29',dec='-46:05:23.7',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=PN%20Tc%201%20')
TMC1 = Source("TMC-1",type='Dark Cloud',ra='04:41:45.9',dec='+25:41:27',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=TMC-1')
VYCaMaj = Source("VY Ca Maj",type='Oxygen Star',ra='07:22:58.3',dec='-25:46:03.2',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=VY+Canis+Majoris')
W3 = Source("W3",type='LOS Cloud',ra='02:27:04.10',dec='+61:52:27.1',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W3')
W3OH = Source("W3(OH)",type='SFR',ra='02:27:04.1',dec='+61:52:52',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W3+%28OH%29&NbIdent=1')
W31LOS = Source("W31 LOS",type='LOS Cloud',ra='18:10:28.6',dec='-19:55:51',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W31')
W33LOS = Source("W33 LOS",type='LOS Cloud',ra='18:14:14.0',dec='-17:55:50',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W33')
W43LOS = Source("W43 LOS",type='LOS Cloud',ra='18:47:32.4',dec='-01:56:31',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W43')
W44LOS = Source("W44 LOS",type='LOS Cloud',ra='18:56:10.65',dec='+01:13:21.30',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W43')
W49 = Source("W49",type='SFR',ra='19:10:19.6',dec='+09:07:42',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W49')
W49LOS = Source("W49 LOS",type='LOS Cloud',ra='19:10:19.6',dec='+09:07:42',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W49')
W51 = Source("W51",type='SFR',ra='19:23:50',dec='+14:06:0',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W51')
W51LOS = Source("W51 LOS",type='LOS Cloud',ra='19:23:50',dec='+14:06:0',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=W51')
XiPerLOS = Source("Xi Per LOS",type='LOS Cloud',ra='03:58:57.9',dec='+35:47:27.74',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Xi+Per')
rhoOphA = Source("rho Oph A",type='SFR',ra='16:25:35.14',dec='-23:26:49.9',simbad_url='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Rho+Oph+A')

source_list = [
	AFGL890LOS,
	AFGL961LOS,
	AFGL989LOS,
	B1b,
	CRL2688,
	CRL618,
	CasALOS,
	CrabNebula,
	CygnusOB212LOS,
	DR21,
	DR21LOS,
	DR21OH,
	G0693,
	G327306LOS,
	GL2136LOS,
	GalacticCenter,
	HD124314LOS,
	HD27778LOS,
	HorseheadPDR,
	IC443G,
	IRAS16293,
	IRC10216,
	K350,
	L134,
	L1527,
	L1544,
	L183,
	L483,
	LOSCloud,
	Lupus1A,
	M17LOS,
	M17SW,
	M3LOS,
	NGC2024,
	NGC2024LOS,
	NGC2264,
	NGC6334,
	NGC6334LOS,
	NGC7023,
	NGC7027,
	NGC7538,
	NGC7538LOS,
	Orion,
	OrionBar,
	rhoOphA,
	SgrA,
	SgrALOS,
	SgrB2,
	SgrB2LOS,
	TC1,
	TMC1,
	VYCaMaj,
	W3,
	W3OH,
	W31LOS,
	W33LOS,
	W43LOS,
	W44LOS,
	W49,
	W49LOS,
	W51,
	W51LOS,
	XiPerLOS,
	rhoOphA,
	]


#############################################################
#						Molecule Class 						#
#############################################################

class Molecule(object):

	def __init__(self,name,formula,year,label,sources,telescopes,wavelengths,other_names='',neutral=False,cation=False,anion=False,radical=False,cyclic=False,fullerene=False,pah=False,mass=0,du=0,natoms=0,Acon=None,Bcon=None,Ccon=None,mua=None,mub=None,muc=None,kappa=None,H=0,He=0,C=0,O=0,N=0,S=0,P=0,Si=0,Cl=0,F=0,Mg=0,Na=0,Al=0,K=0,Fe=0,Ti=0,Ar=0,V=0,Ca=0,d_ref=None,lab_ref=None,notes=None,ice=False,ice_d_ref=None,ice_l_ref=None,ppd=None,exgal=None,exo=None,isos=None,isomers=None,ppd_isos=None,ppd_d_ref=None,ppd_l_ref=None,ppd_isos_ref=None,exgal_d_ref=None,exgal_l_ref=None,exo_d_ref=None,exo_l_ref=None,exgal_sources=None,isos_d_ref=None,isos_l_ref=None):
	
		self.name = name
		self.formula = formula
		self.year = year
		self.label = label
		self.sources = sources
		self.telescopes = telescopes
		self.wavelengths = wavelengths
		self.other_names = other_names		
		self.neutral = neutral
		self.cation = cation
		self.anion = anion
		self.radical = radical
		self.cyclic = cyclic
		self.fullerene = fullerene
		self.pah = pah
		self.mass = mass
		self.du = du
		self.natoms = natoms
		self.Acon = Acon
		self.Bcon = Bcon
		self.Ccon = Ccon
		self.mua = mua
		self.mub = mub
		self.muc = muc
		self.kappa = kappa		
		self.H = H
		self.He = He
		self.C = C
		self.O = O
		self.N = N
		self.S = S
		self.P = P
		self.Si = Si
		self.Cl = Cl
		self.F = F
		self.Mg = Mg
		self.Na = Na
		self.Al = Al
		self.K = K
		self.Fe = Fe
		self.Ti = Ti
		self.Ar = Ar
		self.V = V
		self.Ca = Ca
		self.ice = ice
		self.ppd = ppd
		self.exgal = exgal
		self.exo = exo
		self.isos = isos
		self.isomers = isomers
		self.ppd_isos = ppd_isos
		self.d_ref = d_ref
		self.lab_ref = lab_ref
		self.ice_d_ref = ice_d_ref
		self.ice_l_ref = ice_l_ref
		self.ppd_d_ref = ppd_d_ref
		self.ppd_l_ref = ppd_l_ref
		self.ppd_isos_ref = ppd_isos_ref
		self.exgal_d_ref = exgal_d_ref
		self.exgal_l_ref = exgal_l_ref
		self.exgal_sources = exgal_sources
		self.exo_d_ref = exo_d_ref
		self.exo_l_ref = exo_l_ref
		self.isos_d_ref = isos_d_ref
		self.isos_l_ref = isos_l_ref
		self.notes = notes
		self.maxdu = None
		
		self.update_stats()
		
		return
		
	def update_stats(self):
	
		#calculate the number of atoms
	
		self.natoms = self.H + self.He + self.C + self.O + self.N + self.S + self.P + self.Si + self.Cl + self.F + self.Mg + self.Na + self.Al + self.K + self.Fe + self.Ti + self.Ar + self.V + self.Ca
		
		#calculate the mass
		
		self.mass = self.H + self.He*4 + self.C*12 + self.O*16 + self.N*14 + self.S*32 + self.P*31 + self.Si*28 + self.Cl*35 + self.F*19 + self.Mg*24 + self.Na*23 + self.Al*27 + self.K*39 + self.Fe*56 + self.Ti*48 + self.Ar*36 + self.V*51 + self.Ca*40

		#if this is a carbon-bearing molecule and has only H, O, N, and/or halogens, calculate the degree of unsaturation

		if self.C != 0 and all(int(i) is 0 for i in [self.Si, self.Mg, self.Na, self.Al, self.K, self.Fe, self.Ti, self.Ar, self.P, self.He, self.V, self.Ca]):
		
			self.du = 1 + 0.5*(self.H*-1 + self.C*2 + self.N*1 + self.Cl*-1 + self.F*-1)
			self.maxdu = 1 + 0.5*(self.C*2 + self.N*1)
			
		else:
		
			self.du = None
			
		#calculate kappa
		
		if True in (t != None for t in [self.Acon,self.Bcon,self.Ccon]):
	
			A = self.Acon
			B = self.Bcon
			C = self.Ccon
		
			if 	A == None and C == None:
		
				self.kappa = -1
			
			else:
		
				self.kappa = (2*B - A - C)/(A - C)		
				
		return
		
CH = Molecule('methylidyne','CH',1937,'CH',[LOSCloud],[MtWilson],['UV', 'Vis'],neutral=True,H=1,C=1,d_ref='Dunham 1937 PASP 49, 26; Swings & Rosenfeld 1937 ApJ 86, 483; McKellar 1940 PASP 52, 187',lab_ref='Jevons 1932 Phys Soc. pp 177-179; Brazier & Brown 1983 JCP 78, 1608',notes='*First radio in Rydbeck et al. 1973 Nature 246, 466',exgal=True,exgal_d_ref='Whiteoak et al. 1980 MNRAS 190, 17',exgal_sources='LMC, NGC 4945, NGC 5128',Bcon=425476,mua=1.5)
CN = Molecule('cyano radical','CN',1940,'CN',[LOSCloud],[MtWilson],['UV'],radical=True,neutral=True,C=1,N=1,d_ref='McKellar 1940 PASP 52, 187',lab_ref='Poletto and Rigutti 1965 Il Nuovo Cimento 39, 519; Dixon & Woods 1977 JCP 67, 3956; Thomas & Dalby 1968 Can. J. Phys. 46, 2815',notes='*First radio in Jefferts et al. 1970 ApJ 161, L87',ppd=True,ppd_d_ref='Kastner et al. 1997 Science 277, 67; Dutrey et al. 1997 A&A 317, L55',ppd_isos='C15N',ppd_isos_ref='[C15N] Hily-Blant et al. 2017 A&A 603, L6',exgal=True,exgal_d_ref='Henkel et al. 1988 A&A 201, L23',exgal_sources='M82, NGC 253, IC 342',Bcon=56693,mua=1.5)
CHp = Molecule('methylidyne cation','CH+',1941,'CH+',[LOSCloud],[MtWilson],['UV', 'Vis'],cation=True,H=1,C=1,d_ref='Douglas & Herzberg 1941 ApJ 94, 381; Dunham 1937 PASP 49, 26',lab_ref='Douglas & Herzberg 1941 ApJ 94, 381',notes=None,ppd=True,ppd_d_ref='Thi et al. 2011 A&A 530, L2',exgal=True,exgal_d_ref='Magain & Gillet 1987 A&A 184, L5',exgal_sources='LMC',Bcon=417617,mua=1.7)
OH = Molecule('hydroxyl radical','OH',1963,'OH',[CasALOS],[Millstone],['cm'],radical=True,neutral=True,H=1,O=1,d_ref='Weinreb et al. 1963 Nature 200, 829',lab_ref='Ehrenstein et al. 1959 PRL 3, 40',notes=None,ppd=True,ppd_d_ref='Mandell et al. 2008 ApJ 681, L25; Salyk et al. 2008 ApJ 676, L49',exgal=True,exgal_d_ref='Weliachew 1971 ApJ 167, L47',exgal_sources='M82, NGC 253',Bcon=556174,mua=1.7)
CO = Molecule('carbon monoxide','CO',1970,'CO',[Orion],[NRAO36],['mm'],neutral=True,C=1,O=1,d_ref='Wilson et al. 1970 ApJ 161, L43',lab_ref='Cord et al. 1968 Microwave Spectral Tables V5',notes=None,ice=True,ice_d_ref='Soifer et al. 1979 ApJ 232, L53',ice_l_ref='Mantz et al. 1975 JMS 57, 155',ppd=True,ppd_d_ref='Beckwith et al. 1986 ApJ 309, 755',ppd_isos='13CO, C18O, C17O',ppd_isos_ref='[13CO] Sargent & Beckwith 1987 ApJ 323, 294 [C18O] Dutrey et al. 1994 A&A 286, 149 [C17O] Smith et al. 2009 ApJ 701, 163; Guilloteau et al. 2013 A&A 549, A92',exo=True,exo_d_ref='Madhusudhan et al. 2011 Nature 469, 64; Barman et al. 2011 ApJ 733, 65; Lanotte et al. 2014 A&A 572, A73; Barman et al. 2015 ApJ 804, 61',exgal=True,exgal_d_ref='Rickard et al. 1975 ApJ 199, L75',exgal_sources='M82, NGC 253',Bcon=57636,mua=0.1)
H2 = Molecule('hydrogen','H2',1970,'H2',[XiPerLOS],[Aerobee],['UV'],neutral=True,H=2,d_ref='Carruthers 1970 ApJ 161, L81',lab_ref='Carruthers 1970 ApJ 161, L81',notes=None,ppd=True,ppd_d_ref='Thi et al. 1999 ApJ 521, L63',ppd_isos='HD',ppd_isos_ref='[HD] Bergin et al. 2013 Nature 493, 644',exgal=True,exgal_d_ref='Thompson et al. 1978 ApJ 222, L49',exgal_sources='NGC 1068',mua=0.0)
SiO = Molecule('silicon monoxide','SiO',1971,'SiO',[SgrB2],[NRAO36],['mm'],neutral=True,O=1,Si=1,d_ref='Wilson et al. 1971 ApJ 167, L97',lab_ref='Törring 1968 Z. Naturforschung 23A, 777; Raymonda et al. 1970 JCP 52, 3458',notes=None,exgal=True,exgal_d_ref='Mauersberger & Henkel 1991 A&A 245, 457',exgal_sources='NGC 253',Bcon=21712,mua=3.1)
CS = Molecule('carbon monosulfide','CS',1971,'CS',[Orion, W51, IRC10216, DR21],[NRAO36],['mm'],neutral=True,C=1,S=1,d_ref='Penzias et al. 1971 ApJ 168, L53',lab_ref='Mockler & Bird 1955 Phys Rev 98, 1837',notes=None,ppd=True,ppd_d_ref='Ohashi et al. 1991 AJ 102, 2054; Blake et al. 1992 ApJ 391, L99; Guilloteau et al. 2012 A&A 548, A70',ppd_isos='C34S',ppd_isos_ref='[C34S] Artur de la Villarmois et al. 2018 A&A 614, A26',exgal=True,exgal_d_ref='Henkel & Bally 1985 A&A 150, L25',exgal_sources='M82, IC 342',Bcon=24496,mua=2.0)
SO = Molecule('sulfur monoxide','SO',1973,'SO',[Orion],[NRAO36],['mm'],neutral=True,O=1,S=1,d_ref='Gottlieb & Ball 1973 ApJ 184, L59',lab_ref='Winnewisser et al. 1964 JCP 41, 1687',notes=None,ppd=True,ppd_d_ref='Fuente et al. 2010 A&A 524, A19',exgal=True,exgal_d_ref='Johansson 1991 Proc. IAU Symposium 146, 1; Petuchowski & Bennett 1992 ApJ 391, 137',exgal_sources='M82, NGC 253',Bcon=21524,mua=1.5)
SiS = Molecule('silicon monosulfide','SiS',1975,'SiS',[IRC10216],[NRAO36],['mm'],neutral=True,S=1,Si=1,d_ref='Morris et al. 1975 ApJ 199, L47',lab_ref='Hoeft 1965 Z. fur Naturforschung A, A20, 1327',notes=None,Bcon=9077,mua=1.7)
NS = Molecule('nitrogen monosulfide','NS',1975,'NS',[SgrB2],[NRAO36],['mm'],neutral=True,N=1,S=1,d_ref='Gottlieb et al. 1975 ApJ 200, L147; Kuiper et al. 1975 ApJ 200, L151',lab_ref='Amano et al. 1969 JMS 32, 97',notes=None,exgal=True,exgal_d_ref='Martin et al. 2003 A&A 411, L465',exgal_sources='NGC 253',Bcon=23155,mua=1.8)
C2 = Molecule('dicarbon','C2',1977,'C2',[CygnusOB212LOS],[MtHopkins],['IR'],neutral=True,C=2,d_ref='Souza and Lutz 1977 ApJ 216, L49',lab_ref='Phillips 1948 ApJ 107, 389',exgal=True,exgal_d_ref='Welty et al. 2012 MNRAS 428, 1107', exgal_sources='SMC', notes=None,mua=0.0)
NO = Molecule('nitric oxide','NO',1978,'NO',[SgrB2],[NRAO36],['mm'],neutral=True,O=1,N=1,d_ref='Liszt and Turner 1978 ApJ 224, L73',lab_ref='Gallagher & Johnson 1956 Phys Rev 103, 1727',notes=None,exgal=True,exgal_d_ref='Martin et al. 2003 A&A 411, L465',exgal_sources='NGC 253',Bcon=50849,mua=0.2)
HCl = Molecule('hydrogen chloride','HCl',1985,'HCl',[Orion],[Kuiper],['sub-mm'],neutral=True,H=1,Cl=1,d_ref='Blake et al. 1985 ApJ 295, 501',lab_ref='de Lucia et al. 1971 Phys Rev A 3, 1849',exgal=True,exgal_d_ref='Wallstrom et al. 2019 A&A 629, A128',exgal_sources='PKS 1830-211',notes=None,Bcon=312989,mua=1.1)
NaCl = Molecule('sodium chloride','NaCl',1987,'NaCl',[IRC10216],[IRAM30],['mm'],neutral=True,Cl=1,Na=1,d_ref='Cernicharo & Guélin 1987 A&A 183, L10',lab_ref='Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609',notes=None,Bcon=6513,mua=9.0)
AlCl = Molecule('aluminum chloride','AlCl',1987,'AlCl',[IRC10216],[IRAM30],['mm'],neutral=True,Cl=1,Al=1,d_ref='Cernicharo & Guélin 1987 A&A 183, L10',lab_ref='Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609',notes=None,Bcon=7289,mua='*')
KCl = Molecule('potassium chloride','KCl',1987,'KCl',[IRC10216],[IRAM30],['mm'],neutral=True,Cl=1,K=1,d_ref='Cernicharo & Guélin 1987 A&A 183, L10',lab_ref='Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609',notes=None,Bcon=3845,mua=10.3)
AlF = Molecule('aluminum fluoride','AlF',1987,'AlF',[IRC10216],[IRAM30],['mm'],neutral=True,F=1,Al=1,d_ref='Cernicharo & Guélin 1987 A&A 183, L10',lab_ref='Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609',notes='*Confirmed in 1994 ApJ 433, 729',isos='26AlF',isos_d_ref='[26AlF] Kamiński et al. 2018 Nature Astronomy 2, 778',Bcon=16488,mua=1.5)
PN = Molecule('phosphorous mononitride','PN',1987,'PN',[TMC1, Orion, W51, SgrB2],[NRAOARO12, FCRAO14m, OVRO],['mm'],neutral=True,N=1,P=1,d_ref='Sutton et al. 1985 ApJS 58, 341',lab_ref='Wyse et al. 1972 JCP 57, 1106',notes='*Confirmed in Turner & Bally 1987 ApJ 321, L75 and Ziurys 1987 ApJ 321 L81',Bcon=23495,mua=2.7)
SiC = Molecule('silicon carbide','SiC',1989,'SiC',[IRC10216],[IRAM30],['mm'],radical=True,neutral=True,C=1,Si=1,d_ref='Cernicharo et al. 1989 ApJ 341, L25',lab_ref='Cernicharo et al. 1989 ApJ 341, L25',notes=None,Bcon=20298,mua=1.7)
CP = Molecule('carbon monophosphide','CP',1990,'CP',[IRC10216],[IRAM30],['mm'],radical=True,neutral=True,C=1,P=1,d_ref='Guélin et al. 1990 A&A 230, L9',lab_ref='Saito et al. 1989 ApJ 341, 1114',notes=None,Bcon=23860,mua='*')
NH = Molecule('imidogen radical','NH',1991,'NH',[XiPerLOS, HD27778LOS],[KPNO4m, IRAM30],['UV'],radical=True,neutral=True,H=1,N=1,d_ref='Meyer & Roth 1991 ApJ 376, L49',lab_ref='Dixon 1959 Can J. Phys. 37, 1171 and Klaus et al. 1997 A&A 322, L1',notes='*First radio in Cernicharo et al. 2000 ApJ 534, L199',exgal=True,exgal_d_ref='Gonzalez-Alfonso et al. 2004 ApJ 613, 247',exgal_sources='Arp 220',Bcon=489959,mua=1.4)
SiN = Molecule('silicon nitride ','SiN',1992,'SiN',[IRC10216],[NRAOARO12],['mm'],radical=True,neutral=True,N=1,Si=1,d_ref='Turner 1992 ApJ 388, L35',lab_ref='Saito et al. 1983 JCP 78, 6447',notes=None,Bcon=21828,mua=2.6)
SOp = Molecule('sulfur monoxide cation','SO+',1992,'SO+',[IC443G],[NRAOARO12],['mm'],cation=True,radical=True,O=1,S=1,d_ref='Turner 1992 ApJ 396, L107',lab_ref='Amano et al. 1991 JMS 146, 519',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS',Bcon=23249,mua='*')
COp = Molecule('carbon monoxide cation','CO+',1993,'CO+',[M17SW, NGC7027],[NRAOARO12],['mm'],cation=True,C=1,O=1,d_ref='Latter et al. 1993 ApJ 419, L97',lab_ref='Sastry et al. 1981 ApJ 250, L91',notes=None,exgal=True,exgal_d_ref='Fuente et al. 2006 ApJ 641, L105',exgal_sources='M82',Bcon=58983,mua=2.6)
HF = Molecule('hydrogen fluoride','HF',1997,'HF',[SgrB2LOS],[ISO],['IR'],neutral=True,H=1,F=1,d_ref='Neufeld et al. 1997 ApJ 488, L141',lab_ref='Nolt et al. 1987 JMS 125, 274',notes=None,exgal=True,exgal_d_ref='van der Werf et al. 2010 A&A 518, L42; Rangwala et al. 2011 ApJ 743, 94; Monje et al. 2011 ApJL 742, L21',exgal_sources='Mrk 231, Arp 220, Cloverleaf LOS',Bcon=616365,mua=1.8)
N2 = Molecule('nitrogen','N2',2004,'N2',[HD124314LOS],[FUSE],['UV'],neutral=True,N=2,d_ref='Knauth et al. 2004 Nature 409, 636',lab_ref='Stark et al. 2000 ApJ 531, 321',notes=None,mua=0.0)
CFp = Molecule('fluoromethylidynium cation','CF+',2006,'CF+',[OrionBar],[IRAM30, APEX],['mm'],cation=True,C=1,F=1,d_ref='Neufeld et al. 2006 A&A 454, L37',lab_ref='Plummer et al. 1986 JCP 84, 2427',notes=None,exgal=True,exgal_d_ref='Muller et al. 2016 A&A 589, L5',exgal_sources='PKS 1830-211 LOS',Bcon=51294,mua=1.1)
PO = Molecule('phosphorous monoxide','PO',2007,'PO',[VYCaMaj],[SMT10],['mm'],neutral=True,O=1,P=1,d_ref='Tenenbaum et al. 2007 ApJ 666, L29',lab_ref='Bailleux et al. 2002 JMS 216, 465',notes=None,Bcon=21900,mua=1.9)
O2 = Molecule('oxygen','O2',2007,'O2',[Orion, rhoOphA],[Odin, Herschel],['mm', 'sub-mm'],neutral=True,O=2,d_ref='Larsson et al. 2007 A&A 466, 999',lab_ref='Endo & Mizushima 1982 Jpn J Appl Phys 21, L379; Drouin et al. 2010 J Quant Spec Rad Transf 111, 1167',notes='*Also Larsson et al. 2007 A&A 466, 999; Tentative in Goldsmith 2002 ApJ 576, 814',mua=0.0)
AlO = Molecule('aluminum monoxide','AlO',2009,'AlO',[VYCaMaj],[SMT10],['mm'],neutral=True,O=1,Al=1,d_ref='Tenenbaum & Ziurys 2009 ApJ 693, L59',lab_ref='Yamada et al. 1990, JCP 92, 2146',notes=None,Bcon=19142,mua=4.6)
CNm = Molecule('cyanide anion','CN-',2010,'CN-',[IRC10216],[IRAM30],['mm'],anion=True,C=1,N=1,d_ref='Agúndez et al. 2010 A&A 517, L2',lab_ref='Amano 2008 JCP 129, 244305',notes=None,Bcon=56133,mua=0.7)
OHp = Molecule('hydroxyl cation','OH+',2010,'OH+',[SgrB2LOS],[APEX],['sub-mm'],cation=True,H=1,O=1,d_ref='Wyrowski et al. 2010 A&A 518, A26; Gerin et al. 2010 A&A 518, L110; Benz et al. 2010 A&A 521, L35',lab_ref='Bekooy et al. 1985 JCP 82, 3868',notes=None,exgal=True,exgal_d_ref='van der Werf et al. 2010 A&A 518, L42; Rangwala et al. 2011 ApJ 743, 94; Gonzalez-Alfonso et al. 2013 A&A 550, A25',exgal_sources='Mrk 231, Arp 220, NGC 4418',Bcon=492346,mua=2.3)
SHp = Molecule('sulfanylium cation','SH+',2011,'SH+',[SgrB2],[Herschel],['sub-mm'],cation=True,H=1,S=1,d_ref='Benz et al. 2010 A&A 521, L35',lab_ref='Brown et al. 2009 JMS 255, 68',notes='*Also in Menten et al. 2011 A&A 525, A77',exgal=True,exgal_d_ref='Muller et al. 2017 A&A 606, A109',exgal_sources='PKS 1830-211',Bcon=273810,mua=1.3)
HClp = Molecule('hydrogen chloride cation','HCl+',2012,'HCl+',[W31LOS, W49LOS],[Herschel],['sub-mm'],cation=True,H=1,Cl=1,d_ref='de Luca et al. 2012 ApJ 751, L37',lab_ref='Gupta et al. 2012 ApJ 751, L38',notes=None,Bcon=293444,mua=1.8)
SH = Molecule('mercapto radical','SH',2012,'SH',[W49LOS],[SOFIA],['sub-mm'],radical=True,neutral=True,H=1,S=1,d_ref='Neufeld et al. 2012 A&A 542, L6',lab_ref='Morino & Kawaguchi 1995 JMS 170, 172; Klisch et al. 1996 ApJ 473, 1118',notes=None,Bcon=283588,mua=0.8)
TiO = Molecule('titanium monoxide','TiO',2013,'TiO',[VYCaMaj],[SMA],['mm'],neutral=True,O=1,Ti=1,d_ref='Kamiński et al. 2013 A&A 551, A113',lab_ref='Nakimi et al. 1998 JMS 191, 176',notes=None,exo=True,exo_d_ref='Haynes et al. 2015 ApJ 806, 146; Sedaghati et al. 2017 Nature 549, 238; Nugroho et al. 2017 ApJ 154, 221',Bcon=16004,mua=3.3)
ArHp = Molecule('argonium','ArH+',2013,'ArH+',[CrabNebula],[Herschel],['sub-mm'],cation=True,H=1,Ar=1,d_ref='Barlow et al. 2013 Science 342, 1343',lab_ref='Barlow et al. 2013 Science 342, 1343',notes=None,exgal=True,exgal_d_ref='Muller et al. 2015 A&A 582, L4',exgal_sources='PKS 1830-211 LOS',Bcon=307966,mua=2.2)
NSp = Molecule('nitrogen sulfide cation','NS+',2018,'NS+',[B1b, TMC1, L483],[IRAM30],['mm'],cation=True,N=1,S=1,d_ref='Cernicharo et al. 2018 ApJL 853, L22',lab_ref='Cernicharo et al. 2018 ApJL 853, L22',notes=None,Bcon=25050,mua=2.2)
HeHp = Molecule('helium hydride cation',
					'HeH+',
					2019,
					'HeH+',
					[NGC7027],
					[SOFIA],
					['sub-mm'],
					cation=True,
					He=1,
					H=1,
					d_ref='Gusten et al. 2019 Nature 568, 357',
					lab_ref='Perry et al. 2014 JCP 141, 101101',
					notes='Dipole moment from Engel et al. 2005 MNRAS 357, 471 using equilibrium internuclear separation value of 1.45 from Peyerimhoff 1965 JCP 43, 998',
					Bcon=1006063,
					mua=1.7)

VO = Molecule(name='vanadium oxide',
				formula='VO',
				year=2019,
				label='VO',
				sources=[VYCaMaj],
				telescopes=[Hubble],
				wavelengths=['IR'],
				neutral=True,
				O=1,
				V=1,
				d_ref='Humphreys et al. 2019 ApJL 874, L26',
				lab_ref='Adam et al. 1995 JMS 170, 94; Cheung et al. 1994 JMS 163, 443',
				notes='Some VO transitions originally observed in the source by Wallerstein & Gonzalez 2001 PASP 113, 954, Wallerstein 1971 ApJ 169, 195, and Wallerstein 1986 ApJ 164, 101, but not assigned as circumstellar until now.  B constant from Cheung et al. 1982 JMS 91, 165.  Dipole moment from Suenram et al. 1991 JMS 148, 114.',
				Bcon=16381,
				mua=3.355)

#############################################################
#					Three Atom Molecules					#
#############################################################	

H2O = Molecule('water','H2O',1969,'H2O',[SgrB2, Orion, W49],[HatCreek],['cm'],neutral=True,H=2,O=1,d_ref='Cheung et al. 1969 Nature 221, 626',lab_ref='Golden et al. 1948 Phys Rev 73, 92',notes=None,ice=True,ice_d_ref='Gillett & Forrest 1973 ApJ 179, 483',ice_l_ref='Irvine & Pollack 1968 Icarus 8, 324',isos='HDO',isos_d_ref='[HDO] Turner et al. 1975 ApJ 198, L125',isos_l_ref='[HDO] de Lucia et al. 1974 J Phys Chem Ref Data 3, 211; Erlandsson & Cox 1956 J Chem Phys 25, 778',ppd=True,ppd_d_ref='Carr et al. 2004 ApJ 603, 213; Hogerheijde et al. 2011 Science 344, 338',exo=True,exo_d_ref='Tinetti et al. 2007 Nature 448, 169; Deming et al. 2014 ApJ 774, 95; Kreidberg et al. 2014 ApJL 793, L27; Kreidberg et al. 2015 ApJ 814, 66; Lockwood et al. 2014 ApJ 783, L29',exgal=True,exgal_d_ref='Churchwell et al. 1977 A&A 54, 969',exgal_sources='M33',Acon=835840,Bcon=435352,Ccon=278139,mub=1.9)
HCOp = Molecule('formylium cation','HCO+',1970,'HCO+',[W3OH, Orion, L134, SgrA, W51],[NRAO36],['mm'],cation=True,H=1,C=1,O=1,d_ref='Buhl & Snyder 1970 Nature 228, 267',lab_ref='Woods et al. 1975 PRL 35, 1269',notes=None,ppd=True,ppd_d_ref='Kastner et al. 1997 Science 277, 67; Dutrey et al. 1997 A&A 317, L55',ppd_isos='DCO+, H13CO+',ppd_isos_ref='[DCO+] van Dishoeck et al. 2003 A&A 400, L1 [H13CO+] van Zadelhoff et al. 2001 A&A 377, 566; van Dishoeck et al. 2003 A&A 400, L1',exgal=True,exgal_d_ref='Stark et al. 1979 ApJ 229, 118',exgal_sources='M82',Bcon=44594,mua=3.9)
HCN = Molecule('hydrogen cyanide','HCN',1971,'HCN',[W3OH, Orion, SgrA, W49, W51, DR21],[NRAO36],['mm'],neutral=True,H=1,C=1,N=1,d_ref='Snyder et al. 1971 ApJ 163, L47',lab_ref='de Lucia & Gordy 1969 Phys Rev 187, 58',notes=None,ppd=True,ppd_d_ref='Kastner et al. 1997 Science 277, 67; Dutrey et al. 1997 A&A 317, L55',ppd_isos='DCN, H13CN, HC15N',ppd_isos_ref='[DCN] Qi et al. 2008 ApJ 681, 1396 [H13CN] Guzman et al. 2015 ApJ 814, 53 [HC15N] Guzman et al. 2015 ApJ 814, 53',exgal=True,exgal_d_ref='Rickard et al. 1977 ApJ 214, 390',exgal_sources='NGC 253, M82',exo=True,exo_d_ref='Hawker et al. 2018 ApJL 863, L11',Bcon=44316,mua=3.0)
OCS = Molecule('carbonyl sulfide','OCS',1971,'OCS',[SgrB2],[NRAO36],['mm'],neutral=True,C=1,O=1,S=1,d_ref='Jefferts et al. 1971 ApJ 168, L111',lab_ref='King & Gordy 1954 Phys Rev 93, 407',notes=None,ice=True,ice_d_ref='Palumbo et al. 1995 ApJ 449, 674; Palumbo et al. 1997 ApJ 479, 839',ice_l_ref='Palumbo et al. 1995 ApJ 449, 674',exgal=True,exgal_d_ref='Mauersberger et al. 1995 A&A 294, 23',exgal_sources='NGC 253',Bcon=6081,mua=0.7)
HNC = Molecule('hydrogen isocyanide','HNC',1972,'HNC',[W51, NGC2264],[NRAO36],['mm'],neutral=True,H=1,C=1,N=1,d_ref='Snyder & Buhl 1972 Annals of the New York Academy of Science 194, 17; Zuckerman et al. 1972 ApJ 173, L125',lab_ref='Blackman et al. 1976 Nature 261, 395',notes=None,ppd=True,ppd_d_ref='Dutrey et al. 1997 A&A 317, L55',exgal=True,exgal_d_ref='Henkel et al. 1988 A&A 201, L23',exgal_sources='IC 342',Bcon=45332,mua=3.1)
H2S = Molecule('hydrogen sulfide','H2S',1972,'H2S',[W3, W3OH, Orion, NGC2264, SgrB2, W51, DR21OH, NGC7538],[NRAO36],['mm'],neutral=True,H=2,S=1,d_ref='Thaddeus et al. 1972 ApJ 176, L73',lab_ref='Cupp et al. 1968 Phys Rev 171, 60',notes=None,exgal=True,exgal_d_ref='Hekkila et al. 1999 A&A 344, 817',exgal_sources='LMC',ppd=True,ppd_d_ref='Phuong et al. 2018 A&A 616, L5',Acon=310584,Bcon=270368,Ccon=141820,mub=1.0)
N2Hp = Molecule('protonated nitrogen','N2H+',1974,'N2H+',[SgrB2, DR21, NGC6334, NGC2264],[NRAO36],['mm'],cation=True,H=1,N=2,d_ref='Turner 1974 ApJ 193, L83; Green et al. 1974 ApJ 193, L89; Thaddues & Turner 1975 ApJ 201, L25',lab_ref='Saykally et al. 1976 ApJ 205, L101',notes=None,ppd=True,ppd_d_ref='Qi et al. 2003 ApJ 597, 986; Dutrey et al. 2007 A&A 464, 615',ppd_isos='N2D+',ppd_isos_ref='[N2D+] Huang et al. 2015 ApJL 809, L26',exgal=True,exgal_d_ref='Mauersberger & Henkel 1991 A&A 245, 457',exgal_sources='NGC 253, Maffei 2, IC 342, M82, NGC 6946',Bcon=46587,mua=3.4)
C2H = Molecule('ethynyl radical','C2H',1974,'C2H',[Orion],[NRAO36],['mm'],radical=True,neutral=True,H=1,C=2,d_ref='Tucker et al. 1974 ApJ 193, L115',lab_ref='Sastry et al. 1981 ApJ 251, L119',notes=None,ppd=True,ppd_d_ref='Dutrey et al. 1997 A&A 317, L55',exgal=True,exgal_d_ref='Henkel et al. 1988 A&A 201, L23',exgal_sources='M82',Bcon=43675,mua=0.8)
SO2 = Molecule('sulfur dioxide','SO2',1975,'SO2',[Orion, SgrB2],[NRAO36],['mm'],neutral=True,O=2,S=1,d_ref='Snyder et al. 1975 ApJ 198, L81',lab_ref='Steenbeckeliers 1968 Ann. Soc. Sci. Brux 82, 331',notes=None,exgal=True,exgal_d_ref='Martin et al. 2003 A&A 411, L465',exgal_sources='NGC 253',Acon=60779,Bcon=10318,Ccon=8800,mub=1.6)
HCO = Molecule('formyl radical','HCO',1976,'HCO',[W3, NGC2024, W51, K350],[NRAO36],['mm'],radical=True,neutral=True,H=1,C=1,O=1,d_ref='Snyder et al. 1976 ApJ 208, L91',lab_ref='Saito 1972 ApJ 178, L95',notes=None,exgal=True,exgal_d_ref='Sage & Ziurys 1995 ApJ 447, 625; Garcia-Burillo et al. 2002 ApJ 575, L55',exgal_sources='M82',Acon=7829365,Bcon=44788,Ccon=41930,mua=1.4,mub=0.7)
HNO = Molecule('nitroxyl radical','HNO',1977,'HNO',[SgrB2, NGC2024],[NRAO36],['mm'],neutral=True,H=1,O=1,N=1,d_ref='Ulich et al. 1977 ApJ 217, L105',lab_ref='Saito & Takagi 1973 JMS 47, 99',notes=None,Acon=553899,Bcon=42313,Ccon=39165,mua=1.0,mub=1.3)
HCSp = Molecule('protonated carbon monosulfide','HCS+',1981,'HCS+',[Orion, SgrB2],[NRAO36, Bell7m],['mm'],cation=True,H=1,C=1,S=1,d_ref='Thaddeus et al. 1981 ApJ 246, L41',lab_ref='Gudeman et al. 1981 ApJ 246, L47',notes=None,exgal=True,exgal_d_ref='Muller et al. 2013 A&A 551, A109',exgal_sources='PKS 1830-211 LOS',Bcon=10691,mua=1.9)
HOCp = Molecule('hydroxymethyliumylidene','HOC+',1983,'HOC+',[SgrB2],[FCRAO14m, Onsala20m],['mm'],cation=True,H=1,C=1,O=1,d_ref='Woods et al. 1983 ApJ 270, 583',lab_ref='Gudeman et al. 1982 PRL 48, 1344',notes='*Confirmed in 1995 ApJ 455, L73',exgal=True,exgal_d_ref='Usero et al. 2004 A&A 419, 897',exgal_sources='NGC 1068',Bcon=44744,mua=4.0)
SiC2 = Molecule('silacyclopropynylidene','SiC2',1984,'SiC2',[IRC10216],[NRAO36, Bell7m],['mm'],cyclic=True,radical=True,neutral=True,C=2,Si=1,d_ref='Thaddeus et al. 1984 ApJ 283, L45',lab_ref='Michalopoulos et al. 1984 JCP 80, 3556',notes=None,Acon=52474,Bcon=13157,Ccon=10443,mua=2.4)
C2S = Molecule('dicarbon sulfide','C2S',1987,'C2S',[TMC1, IRC10216, SgrB2],[Nobeyama45, IRAM30],['cm', 'mm'],radical=True,neutral=True,C=2,S=1,d_ref='Saito et al. 1987 ApJ 317, L115',lab_ref='Saito et al. 1987 ApJ 317, L115',notes='*Also Cernicharo et al. 1987 A&A 181, L9',exgal=True,exgal_d_ref='Martin et al. 2006 ApJS 164, 450',exgal_sources='NGC 253',Bcon=6478,mua=2.9)
C3 = Molecule('tricarbon','C3',1988,'C3',[IRC10216],[KPNO4m],['IR'],neutral=True,C=3,d_ref='Hinkle et al. 1988 Science 241, 1319',lab_ref='Gausset et al. 1965 ApJ 142, 45',exgal=True,exgal_d_ref='Welty et al. 2012 MNRAS 428, 1107', exgal_sources='SMC',isos='13CCC, C13CC',isos_d_ref='https://arxiv.org/abs/1911.09751',notes=None,mua=0.0)
CO2 = Molecule('carbon dioxide','CO2',1989,'CO2',[AFGL961LOS, AFGL989LOS, AFGL890LOS], [IRAS],['IR'],neutral=True,C=1,O=2,d_ref='d\'Hendecourt & Jourdain de Muizon 1989 A&A 223, L5; van Dishoeck et al. 1996 A&A 315, L349',lab_ref='d\'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453; Paso et al. 1980 JMS 79, 236; Reichle & Young 1972 Can J Phys 50, 2662',notes='*First detected in ices, then in gas phase',ice=True,ice_d_ref='d\'Hendecourt & Jourdain de Muizon 1989 A&A 223, L5',ice_l_ref='d\'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453',ppd=True,ppd_d_ref='Carr & Najita 2008 Science 319, 1504',exo=True,exo_d_ref='Stevenson et al. 2010 Nature 464, 1161; Madhusudhan et al. 2011 Nature 469, 64; Lanotte et al. 2014 A&A 572, A73',mua=0.0)
CH2 = Molecule('methylene','CH2',1989,'CH2',[Orion],[NRAOARO12],['mm'],radical=True,neutral=True,H=2,C=1,d_ref='Hollis et al. 1989 ApJ 346, 794',lab_ref='Lovas et al. 1983 ApJ 267, L131',notes='*Confirmed in 1995 ApJ 438, 259',Acon=2211494,Bcon=253618,Ccon=215102,mub=0.6)
C2O = Molecule('dicarbon monoxide','C2O',1991,'C2O',[TMC1],[Nobeyama45],['cm'],radical=True,neutral=True,C=2,O=1,d_ref='Ohishi et al. 1991 ApJ 380, L39',lab_ref='Yamada et al. 1985 ApJ 290, L65',notes=None,Bcon=11546,mua=1.3)
MgNC = Molecule('magnesium isocyanide','MgNC',1993,'MgNC',[IRC10216],[IRAM30],['mm'],radical=True,neutral=True,C=1,N=1,Mg=1,d_ref='Guélin et al. 1986 A&A 157, L17',lab_ref='Kawaguchi et al. 1993 ApJ 406, L39',notes='*Actually identified in Kawaguchi et al. 1993 ApJ 406, L39 and Guélin et al. 1993 A&A 280, L19',Bcon=5967,mua=5.2)
NH2 = Molecule('amidogen',
				'NH2',
				1993,
				'NH2',
				[SgrB2LOS],
				[CSO],
				['sub-mm'],
				radical=True,
				neutral=True,
				H=2,
				N=1,
				d_ref='van Dishoeck et al. 1993 ApJ 416, L83',
				lab_ref='Charo et al. 1981 ApJ 244, L111',
				notes=None,
				exgal=True,
				exgal_d_ref='Muller et al. 2014 A&A 566, A112',
				exgal_sources='PKS 1830-211 LOS',
				Acon=710302,
				Bcon=388289,
				Ccon=245014,
				mub=1.8,
				isos = 'NHD, ND2',
				isos_d_ref = 'Melosso et al. 2020 A&A 641, A153',
				isos_l_ref = 'Martin-Drumel et al. 2014 JPCA 118, 1331; Melosso et al. 2017 ApJS 233, 1; Bizzocchi et al. 2020 ApJS 247, 59',			
				)
				
NaCN = Molecule('sodium cyanide','NaCN',1994,'NaCN',[IRC10216],[NRAOARO12],['mm'],neutral=True,C=1,N=1,Na=1,d_ref='Turner et al. 1994 ApJ 426, L97',lab_ref='van Vaals et al. 1984 Chem Phys 86, 147',notes=None,Acon=57922,Bcon=8368,Ccon=7272,mua=8.9)
N2O = Molecule('nitrous oxide','N2O',1994,'N2O',[SgrB2],[NRAOARO12],['mm'],neutral=True,O=1,N=2,d_ref='Ziurys et al. 1994 ApJ 436, L181',lab_ref='Lovas 1978 J Phys Chem Ref Data 7, 1445',notes=None,Bcon=12562,mua=0.2)
MgCN = Molecule('magnesium cyanide','MgCN',1995,'MgCN',[IRC10216],[NRAOARO12, IRAM30],['mm'],radical=True,neutral=True,C=1,N=1,Mg=1,d_ref='Ziurys et al. 1995 ApJ 445, L47',lab_ref='Anderson et al. 1994 ApJ 429, L41',notes=None,Bcon=5095,mua='*')
H3p = Molecule('','H3+',1996,'H3+',[GL2136LOS, W33LOS],[UKIRT],['IR'],cation=True,H=3,d_ref='Geballe & Oka 1996 Nature 384, 334',lab_ref='Oka 1980 PRL 45, 531',notes=None,exgal=True,exgal_d_ref='Geballe et al. 2006 ApJ 644, 907',exgal_sources='IRAS 08572+3915',mua=0.0)
SiCN = Molecule('silicon monocyanide radical','SiCN',2000,'SiCN',[IRC10216],[IRAM30],['mm'],radical=True,neutral=True,C=1,N=1,Si=1,d_ref='Guélin et al. 2000 A&A 363, L9',lab_ref='Apponi et al. 2000 ApJ 536, L55',notes=None,Bcon=5543,mua=2.9)
AlNC = Molecule('aluminum isocyanide','AlNC',2002,'AlNC',[IRC10216],[IRAM30],['mm'],neutral=True,C=1,N=1,Al=1,d_ref='Ziurys et al. 2002 ApJ 564, L45',lab_ref='Robinson et al. 1997 Chem Phys Lett 278, 1',notes=None,Bcon=5985,mua=3.1)
SiNC = Molecule('silicon monoisocyanide','SiNC',2004,'SiNC',[IRC10216],[IRAM30],['mm'],radical=True,neutral=True,C=1,N=1,Si=1,d_ref='Guélin et al. 2004 A&A 426, L49',lab_ref='Apponi et al. 2000 ApJ 536, L55',notes=None,Bcon=6397,mua=2.0)
HCP = Molecule('phosphaethyne','HCP',2007,'HCP',[IRC10216],[IRAM30],['mm'],neutral=True,H=1,C=1,P=1,d_ref='Agúndez et al. 2007 ApJ 662, L91',lab_ref='Bizzocchi et al. 2001 JMS 205, 110',notes='*First attempt 1990 ApJ 365, 59. Confirmed 2008 ApJ 684, 618',Bcon=19976,mua=0.4)
CCP = Molecule('dicarbon phosphide radical','CCP',2008,'CCP',[IRC10216],[NRAOARO12],['mm'],radical=True,neutral=True,C=2,P=1,d_ref='Halfen et al. 2008 ApJ 677, L101',lab_ref='Halfen et al. 2008 ApJ 677, L101',notes=None,Bcon=6373,mua=3.4)
AlOH = Molecule('aluminum hydroxide','AlOH',2010,'AlOH',[VYCaMaj],[NRAOARO12, SMT10],['mm'],neutral=True,H=1,O=1,Al=1,d_ref='Tenenbaum & Ziurys 2010 ApJ 712, L93',lab_ref='Apponi et al. 1993 ApJ 414, L129',notes=None,Bcon=15740,mua=1.0)
H2Op = Molecule('oxidaniumyl','H2O+',2010,'H2O+',[SgrB2, SgrB2LOS, NGC6334, DR21],[Herschel],['sub-mm'],cation=True,H=2,O=1,d_ref='Ossenkopf et al. 2010 A&A 518, L111; Gerin et al. 2010 A&A 518, L110',lab_ref='Strahan et al. 1986 JCP 85, 1252; Murtz et al. 1998 JCP 109, 9744 ',notes=None,exgal=True,exgal_d_ref='Weiss et al. 2010 A&A 521, L1',exgal_sources='M82',Acon=870579,Bcon=372365,Ccon=253878,mub=2.4)
H2Clp = Molecule('chloronium','H2Cl+',2010,'H2Cl+',[SgrB2, SgrB2LOS, NGC6334, NGC6334LOS],[Herschel],['sub-mm'],cation=True,H=2,Cl=1,d_ref='Lis et al. 2010 A&A 521, L9',lab_ref='Araki et al. 2001 JMS 210, 132',notes=None,exgal=True,exgal_d_ref='Muller et al. 2014 A&A 566, L6',exgal_sources='PKS 1830-211 LOS',Acon=337352,Bcon=273588,Ccon=148101,mub=1.9)
KCN = Molecule('potassium cyanide','KCN',2010,'KCN',[IRC10216],[NRAOARO12, SMT10, IRAM30],['mm'],neutral=True,C=1,N=1,K=1,d_ref='Pulliam et al. 2010 ApJ 727, L181',lab_ref='Torring et al. 1980 JCP 73, 4875',notes=None,Acon=58266,Bcon=4940,Ccon=4536,mub=10.0)
FeCN = Molecule('iron cyanide','FeCN',2011,'FeCN',[IRC10216],[NRAOARO12],['mm'],neutral=True,C=1,N=1,Fe=1,d_ref='Zack et al. 2011 ApJ 733, L36',lab_ref='Flory & Ziurys 2011 JCP 135, 184303',notes=None,Bcon=4080,mua=4.5)
HO2 = Molecule('hydroperoxyl radical','HO2',2012,'HO2',[rhoOphA],[IRAM30, APEX],['mm'],radical=True,neutral=True,H=1,O=2,d_ref='Parise et al. 2012 A&A 541, L11',lab_ref='Beers & Howard 1975 JCP 63, 4212; Saito 1977 JMS 65, 229; Charo & de Lucia 1982 JMS 94, 426',notes=None,Acon=610273,Bcon=33514,Ccon=31672,mua=1.4,mub=1.5)
TiO2 = Molecule('titanium dioxide','TiO2',2013,'TiO2',[VYCaMaj],[SMA, PdBI],['mm'],neutral=True,O=1,Ti=1,d_ref='Kamiński et al. 2013 A&A 551, A113',lab_ref='Brunken 2008 APJ 676, 1367; Kania et al. 2011 JMS 268, 173',notes=None,Acon=30521,Bcon=8472,Ccon=6614,mua=6.3)
CCN = Molecule('cyanomethylidyne','CCN',2014,'CCN',[IRC10216],[NRAOARO12, SMT10],['mm'],radical=True,neutral=True,C=2,N=1,d_ref='Anderson & Ziurys 2014 ApJ 795, L1',lab_ref='Anderson et al. 2015 JMS 307, 1',notes=None,Bcon=11939,mua=0.4)
SiCSi = Molecule('disilicon carbide','SiCSi',2015,'SiCSi',[IRC10216],[IRAM30],['mm'],neutral=True,C=1,Si=2,d_ref='Cernicharo et al. 2015 ApJ 806, L3',lab_ref='McCarthy 2015 JPC Lett 6, 2107',notes=None,Acon=64074,Bcon=4396,Ccon=4102,mub=0.9)
S2H = Molecule('hydrogen disulfide','S2H',2017,'S2H',[HorseheadPDR],[IRAM30],['mm'],neutral=True,H=1,S=2,d_ref='Fuente et al. 2017 ApJ 851, L49',lab_ref='Tanimoto et al. 2000 JMS 199, 73',notes=None,Acon=296979,Bcon=7996,Ccon=7777,mua=1.2,mub=0.9)
HCS = Molecule('thioformyl','HCS',2018,'HCS',[L483],[IRAM30],['mm'],radical=True,neutral=True,H=1,C=1,S=1,d_ref='Agúndez et al. 2018 A&A 611, L1',lab_ref='Habara et al. 2002 JCP 116, 9232',notes=None,Acon=954000,Bcon=20359,Ccon=19970,mua=0.4,mub=0.9)
HSC = Molecule('sulfhydryl carbide','HSC',2018,'HSC',[L483],[IRAM30],['mm'],radical=True,neutral=True,H=1,C=1,S=1,d_ref='Agúndez et al. 2018 A&A 611, L1',lab_ref='Habara 2000 JCP 112, 10905',notes=None,Acon=295039,Bcon=22036,Ccon=19564,mua=2.5,mub=1.0)
NCO = Molecule('isocyanate radical','NCO',2018,'NCO',[L483],[IRAM30],['mm'],radical=True,neutral=True,C=1,O=1,N=1,d_ref='Marcelino et al. 2018 A&A 612, L10',lab_ref='Kawaguchi et al. 1985 Mol Phys 55, 341; Saito and Amano 1970 JMS34, 383',notes=None,Bcon=11677,mua=0.6)
CaNC = Molecule(name='calcium isocyanide',
					formula='CaNC',
					year=2019,
					label='CaNC',
					sources=[IRC10216],
					telescopes=[IRAM30],
					wavelengths=['mm'],
					neutral=True,
					Ca=1,
					N=1,
					C=1,
					d_ref='Cernicharo et al. 2019 A&A 627, L4',
					lab_ref='Steimle et al. 1993 ApJ 410, L49; Scurlock et al. 1994 JCP 100, 3497',
					notes='Dipole moment from Steimle et al. 1992 JCP 97, 2909',
					Bcon=4048,
					mua=6.985)

#############################################################
#					Four Atom Molecules						#
#############################################################

NH3 = Molecule('ammonia','NH3',1968,'NH3',[GalacticCenter],[HatCreek],['cm'],neutral=True,H=3,N=1,d_ref='Cheung et al. 1968 PRL 25, 1701',lab_ref='Cleeton & Williams 1934 Phys Rev 45, 234',notes=None,ice=True,ice_d_ref='Lacy et al. 1998 ApJ 501, L105',ice_l_ref='d\'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453',ppd=True,ppd_d_ref='Salinas et al. 2016 A&A 591, A122',exgal=True,exgal_d_ref='Martin & Ho 1979 A&A 74, L7',exgal_sources='IC 342, NGC 253',Acon=298193,Bcon=298193,Ccon=286696,muc=1.5)
H2CO = Molecule('formaldehyde','H2CO',1969,'H2CO',[M17LOS, M3LOS, W49LOS, NGC2024LOS, DR21LOS, W43LOS, W44LOS, W51LOS, SgrALOS, SgrB2LOS, W33LOS, NGC6334LOS, CasALOS],[NRAO140],['cm'],neutral=True,H=2,C=1,O=1,d_ref='Snyder et al. 1969 PRL 22, 679',lab_ref='Shinegari 1967 J Phys Soc Jpn 23, 404',notes=None,ice=True,ice_d_ref='Keane et al. 2001 A&A 376, 254',ice_l_ref='Schutte et al. 1993 Icarus 104, 118',ppd=True,ppd_d_ref='Dutrey et al. 1997 A&A 317, L55',exgal=True,exgal_d_ref='Gardner & Whiteoak 1974 Nature 247, 526',exgal_sources='NGC 253, NGC 4945',Acon=281971,Bcon=38834,Ccon=34004,mua=2.3)
HNCO = Molecule('isocyanic acid','HNCO',1972,'HNCO',[SgrB2],[NRAO36],['cm', 'mm'],neutral=True,H=1,C=1,O=1,N=1,d_ref='Snyder & Buhl 1972 ApJ 177, 619',lab_ref='Kewley et al. 1963 JMS 10, 418',notes=None,exgal=True,exgal_d_ref='Nguyen-Q-Rieu et al. 1991 A&A 241, L33',exgal_sources='NGC 253, Maffei 2, IC 342',Acon=912711,Bcon=11071,Ccon=10911,mua=1.6,mub=1.4)
H2CS = Molecule('thioformaldehyde','H2CS',1973,'H2CS',[SgrB2LOS],[Parkes64],['cm'],neutral=True,H=2,C=1,S=1,d_ref='Sinclair et al. 1973 Aust. J. Phys. 26, 85',lab_ref='Johnson & Powell 1970 Science 169, 679',notes=None,exgal=True,exgal_d_ref='Martin et al. 2006 ApJS 164, 450',exgal_sources='NGC 253',Acon=291292,Bcon=17700,Ccon=16652,mua=1.6)
C2H2 = Molecule('acetylene','C2H2',1976,'C2H2',[IRC10216],[KPNO4m],['IR'],neutral=True,H=2,C=2,d_ref='Ridgway et al. 1976 Nature 264, 345',lab_ref='Baldacci et al. 1973 JMS 48, 600',notes=None,ppd=True,ppd_d_ref='Lahuis et al. 2006 ApJ 66, L145',exgal=True,exgal_d_ref='Matsuura et al. 2002 ApJ 580, L133',exgal_sources='LMC',mua=0.0)
C3N = Molecule('cyanoethynyl radical','C3N',1977,'C3N',[IRC10216, TMC1],[NRAO36, Onsala20m],['cm', 'mm'],radical=True,neutral=True,C=3,N=1,d_ref='Guelin & Thaddeus 1977 ApJ 212, L81',lab_ref='Gottlieb et al. 1983 ApJ 275, 916',notes='*Confirmed in Friberg et al. 1980 ApJ 241, L99',Bcon=4968,mua=2.9)
HNCS = Molecule('isothiocyanic acid','HNCS',1979,'HNCS',[SgrB2],[Bell7m, NRAO36],['mm'],neutral=True,H=1,C=1,N=1,S=1,d_ref='Frerking et al. 1979 ApJ 234, L143',lab_ref='Kewley et al. 1963 JMS 10, 418',notes=None,Acon=1348662,Bcon=5883,Ccon=5847,mua=1.6,mub='*')
HOCOp = Molecule('protonated carbon dioxide','HOCO+',1981,'HOCO+',[SgrB2],[Bell7m],['mm'],cation=True,H=1,C=1,O=2,d_ref='Thaddeus et al. 1981 ApJ 246, L41',lab_ref='Green et al. 1976 Chem Phys 17, 479; Bogey et al. 1984 A&A 138, L11',notes=None,exgal=True,exgal_d_ref='Aladro et al. 2015 A&A 579, A101; Martin et al. 2006 ApJS 164, 450',exgal_sources='NGC 253',Acon=789951,Bcon=10774,Ccon=10609,mua=2.7,mub=1.8)
C3O = Molecule('tricarbon monoxide','C3O',1985,'C3O',[TMC1],[NRAOARO12, FCRAO14m, Nobeyama45],['cm', 'mm'],radical=True,neutral=True,C=3,O=1,d_ref='Matthews et al. 1984 Nature 310, 125',lab_ref='Brown et al. 1983 JACS 105, 6496',notes='*Confirmed in Brown et al. 1985 ApJ 297, 302 and Kaifu et al. 2004 PASJ 56, 69',Bcon=4811,mua=2.4)
lC3H = Molecule('propynylidyne radical','l-C3H',1985,'l-C3H',[TMC1, IRC10216],[NRAO36, Bell7m, FCRAO14m, Onsala20m],['cm', 'mm'],radical=True,neutral=True,H=1,C=3,d_ref='Thaddeus et al. 1985 ApJ 294, L49',lab_ref='Gottlieb et al. 1985 ApJ 294, L55',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS',Bcon=11189,mua=3.6,mub=0.5)
HCNHp = Molecule('protonated hydrogen cyanide','HCNH+',1986,'HCNH+',[SgrB2],[NRAOARO12, MWO4m],['mm'],cation=True,H=2,C=1,N=1,d_ref='Ziurys & Turner 1986 ApJ 302, L31',lab_ref='Bogey et al. 1985 JCP 83, 3703; Altman et al. 1984 JCP 80, 3911',notes=None,Bcon=37056,mua=0.3)
H3Op = Molecule('hydronium','H3O+',1986,'H3O+',[Orion, SgrB2],[NRAOARO12, MWO4m],['mm'],cation=True,H=3,O=1,d_ref='Wootten et al. 1986 A&A 166, L15; Hollis et al. 1986 Nature 322, 524',lab_ref='Plummer et al. 1985 JCP 83, 1428; Bogey et al. 1985 A&A 148, L11; Liu & Oka 1985 PRL 54, 1787',notes='*Confirmed in Wootten et al. 1991 ApJ 390, L79',exgal=True,exgal_d_ref='van der Tak et al. 2007 A&A 477, L5',exgal_sources='M82, Arp 220',Acon=334405,Bcon=334405,Ccon=184725,muc=1.4)
C3S = Molecule('tricarbon monosulfide','C3S',1987,'C3S',[TMC1, IRC10216],[Nobeyama45, IRAM30],['cm', 'mm'],radical=True,neutral=True,C=3,S=1,d_ref='Yamamoto et al. 1987 ApJ 317, L119',lab_ref='Yamamoto et al. 1987 ApJ 317, L119',notes=None,Bcon=2890,mua=3.7)
cC3H = Molecule('cyclopropenylidene radical','c-C3H',1987,'c-C3H',[TMC1],[Nobeyama45],['mm'],cyclic=True,radical=True,neutral=True,H=1,C=3,d_ref='Yamamoto et al. 1987 ApJ 322, L55',lab_ref='Yamamoto et al. 1987 ApJ 322, L55',notes=None,exgal='Tentative',exgal_d_ref='Martin et al. 2006 ApJS 164, 450',exgal_sources='NGC 253',Acon=44517,Bcon=34016,Ccon=19189,mua=2.4)
HC2N = Molecule('cyanocarbene radical','HC2N',1991,'HC2N',[IRC10216],[IRAM30],['mm'],radical=True,neutral=True,H=1,C=2,N=1,d_ref='Guélin & Cernicharo 1991 A&A 244, L21',lab_ref='Saito et al. 1984 JCP 80, 1427; Brown et al. 1990 JMS 143, 203',notes=None,Bcon=10986,mua=3.0)
H2CN = Molecule('methylene amidogen radical','H2CN',1994,'H2CN',[TMC1],[NRAOARO12],['cm'],radical=True,neutral=True,H=2,C=1,N=1,d_ref='Ohishi et al. 1994 ApJ 427, L51',lab_ref='Yamamoto & Saito 1992 JCP 96, 4157',notes=None,Acon=284343,Bcon=39158,Ccon=34246,mua=2.5)
SiC3 = Molecule('silicon tricarbide','SiC3',1999,'SiC3',[IRC10216],[NRAOARO12],['mm'],neutral=True,cyclic=True,C=3,Si=1,d_ref='Apponi et al. 1999 ApJ 516, L103',lab_ref='Apponi et al. 1999 JCP 111, 3911; McCarthy et al. JCP 110, 1064',notes=None,Acon=37944,Bcon=6283,Ccon=5387,mua=4.0)
CH3 = Molecule('methyl radical','CH3',2000,'CH3',[SgrALOS],[ISO],['IR'],radical=True,neutral=True,H=3,C=1,d_ref='Feuchtgruber et al. 2000 ApJ 535, L111',lab_ref='Yamada et al. 1981 JCP 75, 5256',notes=None,mua=0.0)
C3Nm = Molecule('cyanoethynyl anion','C3N-',2008,'C3N-',[IRC10216],[IRAM30],['mm'],anion=True,C=3,N=1,d_ref='Thaddeus et al. 2008 ApJ 677, 1132',lab_ref='Thaddeus et al. 2008 ApJ 677, 1132',notes=None,Bcon=4852,mua=3.1)
PH3 = Molecule('phosphine','PH3',2008,'PH3',[IRC10216, CRL2688],[IRAM30, Herschel, SMT10, CSO],['mm', 'sub-mm'],neutral=True,H=3,P=1,d_ref='Agúndez et al. 2008 A&A 485, L33',lab_ref='Cazzoli & Puzzarini 2006 JMS 239, 64; Sousa-Silva et al. 2013 JMS 288, 28; Muller 2013 JQSRT 130, 335',notes='*Confirmed in Agúndez et al. 2014 ApJL 790, L27',Acon=133480,Bcon=133480,Ccon=117488)
HCNO = Molecule('fulminic acid','HCNO',2009,'HCNO',[B1b, L1544, L183, L1527],[IRAM30],['mm'],neutral=True,H=1,C=1,O=1,N=1,d_ref='Marcelino et al. 2009 ApJ 690, L27',lab_ref='Winnewisser & Winnewisser 1971 Z Naturforsch 26, 128',notes=None,Bcon=11469,mua=3.1)
HOCN = Molecule('cyanic acid','HOCN',2009,'HOCN',[SgrB2],[Bell7m, NRAO36, NRAOARO12],['mm'],neutral=True,H=1,C=1,O=1,N=1,d_ref='Brünken et al. 2009 ApJ 697, 880',lab_ref='Brünken et al. 2009 ApJ 697, 880',notes='*Confirmed in Brünken et al. 2010 A&A 516, A109',Acon=681000,Bcon=10577,Ccon=10398,mua=3.7)
HSCN = Molecule('thiocyanic acid','HSCN',2009,'HSCN',[SgrB2],[NRAOARO12],['mm'],neutral=True,H=1,C=1,N=1,S=1,d_ref='Halfen et al. 2009 ApJ 702, L124',lab_ref='Brunken et al. 2009 ApJ 706, 1588',notes=None,Acon=289830,Bcon=5795,Ccon=5675,mua=3.3)
HOOH = Molecule('hydrogen peroxide','HOOH',2011,'HOOH',[rhoOphA],[APEX],['mm'],neutral=True,H=2,O=2,d_ref='Bergman et al. 2011 A&A 531, L8',lab_ref='Petkie et al. 1995 JMS 171, 145; Helminger et al. 1981 JMS 85, 120',notes=None,Acon=301878,Bcon=26212,Ccon=25099,muc=1.6)
lC3Hp = Molecule('cyclopropynylidynium cation','l-C3H+',2012,'l-C3H+',[HorseheadPDR],[IRAM30],['mm'],cation=True,H=1,C=3,d_ref='Pety et al. 2012 A&A 549, A68',lab_ref='Brunken et al. 2014 ApJ 783, L4',notes=None,Bcon=11245,mua=3.0)
HMgNC = Molecule('hydromagnesium isocyanide','HMgNC',2013,'HMgNC',[IRC10216],[IRAM30],['mm'],neutral=True,H=1,C=1,N=1,Mg=1,d_ref='Cabezas et al. 2013 ApJ 75, 133',lab_ref='Cabezas et al. 2013 ApJ 75, 133',notes=None,Bcon=5481,mua=3.5)
HCCO = Molecule('ketenyl radical','HCCO',2015,'HCCO',[Lupus1A, L483],[IRAM30],['mm'],radical=True,neutral=True,H=1,C=2,O=1,d_ref='Agúndez et al. 2015 A&A 577, L5',lab_ref='Endo & Hirota 1987 JCP 86, 4319; Oshima & Endo 1993 JMS 159, 458',notes=None,Bcon=10831,mua=1.6)
CNCN = Molecule('isocyanogen','CNCN',2018,'CNCN',[L483],[IRAM30],['mm'],neutral=True,C=2,N=2,d_ref='Agundez et al. 2018 ApJL 861, L22',lab_ref='Gerry et al. 1990 JMS 140, 147; Winnewisser et al. 1992 JMS 153, 635',notes=None,Bcon=5174,mua=0.7)
HONO = Molecule('nitrous acid','HONO',2019,'HONO',[IRAS16293],[ALMA],['sub-mm'],neutral=True,H=1,O=2,N=1,d_ref='Coutens et al. 2019 A&A 623, L13',lab_ref='Guilmot et al. 1993 JMS 160, 387; Guilmot et al. 1993 JMS 160, 401; Dehayem-Kamadjeu et al. 2005 JMS 234, 182',notes='Only lines of trans-HONO are claimed as detected. As such, constants for this entry are for trans-HONO.',Acon=92892,Bcon=12525,Ccon=11017,mua=1.378,mub=1.242)
MgCCH = Molecule('magnesium ethynyl radical','MgCCH',2019,'MgCCH',[IRC10216],[IRAM30],['mm'],neutral=True,radical=True,Mg=1,C=2,H=1,d_ref='Agundez et al. 2014 A&A 570, A45 (tentative); Cernicharo et al. 2019 A&A 630, L2 (confirmation)',lab_ref='Brewster et al. 1999 Chem. Phys. Lett. 310, 411',notes='Dipole from Woon 1996 ApJ 456, 602',Bcon=4965,mua=1.68)

#############################################################
#					Five Atom Molecules						#
#############################################################

HC3N = Molecule('cyanoacetylene','HC3N',1971,'HC3N',[SgrB2],[NRAO140],['cm'],neutral=True,H=1,C=3,N=1,d_ref='Turner 1971 ApJ 163, L35',lab_ref='Tyler & Sheridan 1963 Trans Faraday Soc 59, 2661',notes='*Confirmed in Dickinson 1972 AL 12, 235',ppd=True,ppd_d_ref='Chapillon et al. 2012 ApJ 756, 58',exgal=True,exgal_d_ref='Mauersberger et al. 1990 A&A 236, 63; Henkel et al. 1988 A&A 201, L23',exgal_sources='NGC 253',Bcon=4549,mua=3.7)
HCOOH = Molecule('formic acid','HCOOH',1971,'HCOOH',[SgrB2],[NRAO140],['cm'],neutral=True,H=2,C=1,O=2,d_ref='Zukerman et al. 1971 ApJ 163, L41',lab_ref='Zukerman et al. 1971 ApJ 163, L41; Bellet et al. 1971 J Mol Struct 9, 49; Bellet et al. 1971 J Mol Struct 9, 65',notes='*Confirmed in Winnewisser & Churchwell 1975 ApJ 200, L33',ice=True,ice_d_ref='Schutte et al. 1999 A&A 343, 966',ice_l_ref='Schutte et al. 1999 A&A 343, 966',ppd=True,ppd_d_ref='Favre et al. 2018 ApJL 862, L2',Acon=77512,Bcon=12055,Ccon=10416,mua=1.4,mub=0.2)
CH2NH = Molecule('methanimine','CH2NH',1973,'CH2NH',[SgrB2],[Parkes64],['cm'],neutral=True,H=3,C=1,N=1,d_ref='Godfrey et al. 1973 ApL 13, 119',lab_ref='Godfrey et al. 1973 ApL 13, 119; Johnson & Lovas 1972 CPL 15, 65',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS',Acon=196211,Bcon=34532,Ccon=29352,mua=1.3,mub=1.5)
NH2CN = Molecule('cyanamide','NH2CN',1975,'NH2CN',[SgrB2],[NRAO36],['mm'],neutral=True,H=2,C=1,N=2,d_ref='Turner et al. 1975 ApJ 201, L149',lab_ref='Tyler et al. 1972 JMS 43, 248; Miller et al. 1962 JMS 8, 153; Lide 1962 JMS 8, 142; Johnson & Suenram 1976 ApJ 208, 245',notes=None,exgal=True,exgal_d_ref='Martin et al. 2006 ApJS 164, 450',exgal_sources='NGC 253',Acon=312142,Bcon=10130,Ccon=9866,mua=4.3,muc=1.0)
H2CCO = Molecule('ketene','H2CCO',1977,'H2CCO',[SgrB2],[NRAO36],['mm'],neutral=True,H=2,C=2,O=1,d_ref='Turner 1977 ApJ 213, L75',lab_ref='Johnson & Strandberg 1952 JCP 20, 687; Johns et al. 1972 JMS 42, 523',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS',Acon=282473,Bcon=10294,Ccon=9916,mua=1.4)
C4H = Molecule('butadiynyl radical','C4H',1978,'C4H',[IRC10216],[NRAO36],['mm'],radical=True,neutral=True,H=1,C=4,d_ref='Guélin et al. 1978 ApJ 224, L27',lab_ref='Gottlieb et al. 1983 ApJ 275, 916',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS',Bcon=4759,mua=0.9)
SiH4 = Molecule('silane','SiH4',1984,'SiH4',[IRC10216],[IRTF],['IR'],neutral=True,H=4,Si=1,d_ref='Goldhaber and Betz 1977 ApJ 279, L55',lab_ref='Goldhaber and Betz 1977 ApJ 279, L55',notes=None,mua=0.0)
cC3H2 = Molecule('cyclopropenylidene','c-C3H2',1985,'c-C3H2',[SgrB2, Orion, TMC1],[Bell7m],['cm', 'mm'],neutral=True,cyclic=True,H=2,C=3,d_ref='Thaddeus et al. 1985 ApJ 299, L63',lab_ref='Thaddeus et al. 1985 ApJ 299, L63',notes='*See also Vrtilek et al. 1987 ApJ 314, 716',ppd=True,ppd_d_ref='Qi et al. 2013 ApJL 765, L14',exgal=True,exgal_d_ref='Seaquist & Bell 1986 ApJ 303, L67',exgal_sources='NGC 5128',Acon=35093,Bcon=32213,Ccon=16749,mub=3.4)
CH2CN = Molecule('cyanomethyl radical','CH2CN',1988,'CH2CN',[TMC1, SgrB2],[FCRAO14m, NRAO140, Onsala20m, Nobeyama45],['cm'],radical=True,neutral=True,H=2,C=2,N=1,d_ref='Irvine et al. 1988 ApJ 334, L107',lab_ref='Saito et al. 1988 ApJ 334, L113',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS',Acon=285130,Bcon=10246,Ccon=9877,mua=1.6)
C5 = Molecule('pentacarbon','C5',1989,'C5',[IRC10216],[KPNO4m],['IR'],neutral=True,C=5,d_ref='Bernath et al. 1989 Science 244, 562',lab_ref='Vala et al. 1989 JCP 90, 595',notes=None,mua=0.0)
SiC4 = Molecule('silicon tetracarbide','SiC4',1989,'SiC4',[IRC10216],[Nobeyama45],['cm', 'mm'],neutral=True,C=4,Si=1,d_ref='Ohishi et al. 1989 ApJ 345, L83',lab_ref='Ohishi et al. 1989 ApJ 345, L83',notes=None,Bcon=1534,mua=6.4)
H2CCC = Molecule('propadienylidene','H2CCC',1991,'H2CCC',[TMC1],[IRAM30,Effelsberg100],['cm', 'mm'],neutral=True,H=2,C=3,d_ref='Cernicharo et al. 1991 ApJ 368, L39',lab_ref='Vrtilek et al. 1990 ApJ 364, L53',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS',Acon=288775,Bcon=10589,Ccon=10204,mua=4.1)
CH4 = Molecule('methane','CH4',1991,'CH4',[NGC7538LOS],[IRTF],['IR'],neutral=True,H=4,C=1,d_ref='Lacy et al. 1991 ApJ 376, 556',lab_ref='Champion et al. 1989 JMS 133, 256; d\'Hendecourt & Allamandola 1986 A&A Supp Ser. 64, 453 ',notes=None,ice=True,ice_d_ref='Lacy et al. 1991 ApJ 376, 556',ice_l_ref='d\'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453',ppd=True,ppd_d_ref='Gibb et al. 2013 ApJL 776, L28',exo=True,exo_d_ref='Swain et al. 2008 Nature 452, 329; Barman et al. 2011 ApJ 733, 65; Stevenson et al. 2014 ApJ 791, 36; Barman et al. 2015 ApJ 804, 61',mua=0.0)
HCCNC = Molecule('isocyanoacetylene','HCCNC',1992,'HCCNC',[TMC1],[Nobeyama45],['cm', 'mm'],neutral=True,H=1,C=3,N=1,d_ref='Kawaguchi et al. 1992 ApJ 386, L51',lab_ref='Kruger et al. 2010 Ang. Chem. 23, 1644',notes=None,Bcon=4968,mua=2.9)
HNCCC = Molecule('','HNCCC',1992,'HNCCC',[TMC1],[Nobeyama45],['cm'],neutral=True,H=1,C=3,N=1,d_ref='Kawaguchi et al. 1992 ApJ 396, L49',lab_ref='Kawaguchi et al. 1992 ApJ 396, L49',notes=None,Bcon=4668,mua=5.7)
H2COHp = Molecule('protonated formaldehyde','H2COH+',1996,'H2COH+',[SgrB2, Orion, W51],[Nobeyama45, NRAOARO12],['cm', 'mm'],cation=True,H=3,C=1,O=1,d_ref='Ohishi et al. 1996 ApJ 471, L61',lab_ref='Chomiak et al. 1994 Can J Phys 72, 1078',notes=None,Acon=197582,Bcon=34351,Ccon=29173,mua=1.4,mub=1.8)
C4Hm = Molecule('butadiynyl anion','C4H-',2007,'C4H-',[IRC10216],[IRAM30],['mm'],anion=True,H=1,C=4,d_ref='Cernicharo et al. 2007 A&A 467, L37',lab_ref='Gupta et al. 2007 ApJ 655, L57',notes=None,Bcon=4655,mua=6.2)
CNCHO = Molecule('cyanoformaldehyde','CNCHO',2007,'CNCHO',[SgrB2],[GBT],['cm'],neutral=True,H=1,C=2,O=1,N=1,d_ref='Remijan et al. 2008 ApJ 675, L85',lab_ref='Bogey et al. 1988 CPL 146, 227; Bogey et al. 1995 JMS 172, 344',notes=None,Acon=67470,Bcon=5010,Ccon=4657,mua=0.8,mub=1.9)
HNCNH = Molecule('carbodiimide','HNCNH',2012,'HNCNH',[SgrB2],[GBT],['cm'],neutral=True,H=2,C=1,N=2,d_ref='McGuire et al. 2012 ApJ 758, L33',lab_ref='Birk et al. 1989 JMS 135, 402; Wagener et al. 1995 JMS 170, 323; Jabs et al. 1997 Chem Phys 225, 77',notes=None,Acon=379244,Bcon=10367,Ccon=10366,mub=1.9)
CH3O = Molecule('methoxy radical','CH3O',2012,'CH3O',[B1b],[IRAM30],['mm'],radical=True,neutral=True,H=3,C=1,O=1,d_ref='Cernicharo et al. 2012 ApJ 759, L43',lab_ref='Momose et al. 1988 JCP 88, 5338; Endo et al. 1984 JCP 81, 122',notes=None,Acon=513887,Bcon=27930,Ccon=27930,mua=2.1)
NH3Dp = Molecule('ammonium ion','NH3D+',2013,'NH3D+',[Orion, B1b],[IRAM30],['mm'],cation=True,H=4,N=1,d_ref='Gupta et al. 2013 ApJ 778, L1',lab_ref='Gupta et al. 2013 ApJ 778, L1',notes='*Confirmed in Marcelino et al. 2018 A&A 612, L10',Acon=175439,Bcon=131412,Ccon=131412,mua=0.3)
H2NCOp = Molecule('protonated isocyanic acid','H2NCO+',2013,'H2NCO+',[SgrB2, L483],[GBT],['cm'],cation=True,H=2,C=1,O=1,N=1,d_ref='Cernicharo et al. 2013 ApJ 771, L10',lab_ref='Cernicharo et al. 2013 ApJ 771, L10',notes='*See also Doménech et al. 2013 ApJ 77, L11',Acon=319800,Bcon=10279,Ccon=9949,mua=4.1)
NCCNHp = Molecule('protonated cyanogen','NCCNH+',2015,'NCCNH+',[TMC1, L483],[IRAM30, Yebes40],['cm', 'mm'],cation=True,H=1,C=2,N=2,d_ref='Agúndez et al. 2015 A&A 579, L10',lab_ref='Amano & Scappini 1991 JCP 95, 2280; Gottlieb et al. 200 JCP 113, 1910',notes=None,Bcon=4438,mua=6.5)
CH3Cl = Molecule('chloromethane','CH3Cl',2017,'CH3Cl',[IRAS16293],[ALMA],['mm'],neutral=True,H=3,C=1,Cl=1,d_ref='Fayolle et al. 2017 Nature Astron. 1, 702',lab_ref='Wlodarczak et al. 1986 JMS 116, 251',notes=None,Acon=156051,Bcon=13293,Ccon=13293,mua=1.9)
MgC3N = Molecule('magnesium cyanoethynyl radical','MgC3N',2019,'MgC3N',[IRC10216],[IRAM30,Yebes40],['cm', 'mm'],neutral=True,radical=True,Mg=1,C=3,N=1,d_ref='Cernicharo et al. 2019 A&A 630, L2',lab_ref='Cernicharo et al. 2019 A&A 630, L2',notes='Assigned based entirely on quantum chemistry; no lab work.',Bcon=1381,mua=6.3)
HC3Op = Molecule('protonated tricarbon monoxide',
					'HC3O+',
					2020,
					'HC3O+',
					[TMC1],
					[IRAM30, Yebes40],
					['cm', 'mm'],
					cation=True,
					H=1,
					C=3,
					O=1,
					d_ref='Cernicharo et al. 2020 A&A 642, L17',
					lab_ref='Cernicharo et al. 2020 A&A 642, L17',
					notes=None,
					Bcon=4461,
					mua=3.4)


#############################################################
#					Six Atom Molecules						#
#############################################################


CH3OH = Molecule('methanol','CH3OH',1970,'CH3OH',[SgrA, SgrB2],[NRAO140],['cm'],neutral=True,H=4,C=1,O=1,d_ref='Ball et al. 1970 ApJ 162, L203',lab_ref='Ball et al. 1970 ApJ 162, L203',notes=None,ice=True,ice_d_ref='Grim et al. 1991 A&A 243, 473',ice_l_ref='d\'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453',ppd=True,ppd_d_ref='Walsh et al. 2016 ApJL 823, L10',exgal=True,exgal_d_ref='Henkel et al. 1987 A&A 188, L1',exgal_sources='NGC 253, IC 342',Acon=127523,Bcon=24690,Ccon=23760,mua=0.9,mub=1.4)
CH3CN = Molecule('methyl cyanide','CH3CN',1971,'CH3CN',[SgrA, SgrB2],[NRAO36],['mm'],neutral=True,H=3,C=2,N=1,d_ref='Solomon et al. 1971 ApJ 168, L107',lab_ref='Cord et al. 1968 Microwave Spectral Tables V5; Kessler et al. Phys Rev 79, 54',notes=None,ppd=True,ppd_d_ref='Oberg et al. 2015 Nature 520, 198',exgal=True,exgal_d_ref='Mauersberger et al. 1991 A&A 247, 307',exgal_sources='NGC 253',Acon=158099,Bcon=9199,Ccon=9199,mua=3.9)
NH2CHO = Molecule('formamide','NH2CHO',1971,'NH2CHO',[SgrB2],[NRAO140],['cm'],neutral=True,H=3,C=1,O=1,N=1,d_ref='Rubin et al. 1971 ApJ 169, L39',lab_ref='Rubin et al. 1971 ApJ 169, L39',notes=None,exgal=True,exgal_d_ref='Muller et al. 2013 A&A 551, A109',exgal_sources='PKS 1830-211 LOS',Acon=72717,Bcon=11373,Ccon=9834,mua=3.6,mub=0.9)
CH3SH = Molecule('methyl mercaptan','CH3SH',1979,'CH3SH',[SgrB2],[Bell7m],['mm'],neutral=True,H=4,C=1,S=1,d_ref='Linke et al. 1979 ApJ 234, L139',lab_ref='Kilb 1955 JCP 23, 1736',notes=None,Acon=102771,Bcon=12952,Ccon=12400,mua=1.3,mub=0.8)
C2H4 = Molecule('ethylene','C2H4',1981,'C2H4',[IRC10216],[McMath],['IR'],neutral=True,H=4,C=2,d_ref='Betz 1981 ApJ 244, L103',lab_ref='Lambeau et al. 1980 JMS 81, 227',notes=None,mua=0.0)
C5H = Molecule('pentynylidyne radical','C5H',1986,'C5H',[IRC10216],[IRAM30],['cm'],radical=True,neutral=True,H=1,C=5,d_ref='Cernicharo et al. 1986 A&A 164, L1',lab_ref='Gottlieb et al. 1986 A&A 164, L5',notes='*See also Cernicharo et al. 1986 A&A 167, L5 and Cernicharo et al. 1987 A&A 172, L5',Bcon=2395,mua=4.9)
CH3NC = Molecule('methyl isocyanide','CH3NC',1988,'CH3NC',[SgrB2],[IRAM30],['mm'],neutral=True,H=3,C=2,N=1,d_ref='Cernicharo et al. 1988 A&A 189, L1',lab_ref='Kukolich 1972 JCP 57, 869; Ring et al. 1947 Phys Rev 72, 1262',notes='*Confirmed in Remijan et al. 2005 ApJ 632, 333 and Gratier et al. 2013 557, A101',Acon=157151,Bcon=10053,Ccon=10053,mua=3.9)
HC2CHO = Molecule('propynal','HC2CHO',1988,'HC2CHO',[TMC1],[NRAO140, Nobeyama45],['cm'],neutral=True,H=2,C=3,O=1,d_ref='Irvine et al. 1988 ApJ 335, L89',lab_ref='Winnewisser 1973 JMS 46, 16',notes=None,Acon=68035,Bcon=4826,Ccon=4500,mua=2.4,mub=0.6)
H2C4 = Molecule('butatrienylidene','H2C4',1991,'H2C4',[IRC10216],[IRAM30],['mm'],neutral=True,H=2,C=4,d_ref='Cernicharo et al. 1991 ApJ 368, L43',lab_ref='Killian et al. 1990 ApJ 365, L89',notes=None,Acon=286234,Bcon=4503,Ccon=4429,mua=4.1)
C5S = Molecule('pentacarbon monosulfide radical','C5S',1993,'C5S',[IRC10216],[NRAO140],['cm'],radical=True,neutral=True,C=5,S=1,d_ref='Bell et al. 1993 ApJ 417, L37',lab_ref='Kasai et al. 1993 ApJ 410, L45; Gordon et al. 2001 ApJS 134, 311',notes='*Confirmed in Agúndez et al. 2014 A&A 570, A45',Bcon=923,mua=5.1)
HC3NHp = Molecule('protonated cyanoacetylene','HC3NH+',1994,'HC3NH+',[TMC1],[Nobeyama45],['cm'],cation=True,H=2,C=3,N=1,d_ref='Kawaguchi et al. 1994 ApJ 420, L95',lab_ref='Lee & Amano 1987 ApJ 323',notes=None,Bcon=4329,mua=1.6)
C5N = Molecule('cyanobutadiynyl radical','C5N',1998,'C5N',[TMC1],[IRAM30,Effelsberg100],['cm'],radical=True,neutral=True,C=5,N=1,d_ref='Guélin et al. 1998 A&A 355, L1',lab_ref='Kasai et al. 1997 ApJ 477, L65',notes=None,Bcon=1403,mua=3.4)
HC4H = Molecule('diacetylene','HC4H',2001,'HC4H',[CRL618],[ISO],['IR'],neutral=True,H=2,C=4,d_ref='Cernicharo et al. 2001 ApJ 546, L123',lab_ref='Arie & Johns 1992 JMS 155, 195',notes='*Confirmed in 2018 ApJ 852, 80',exgal=True,exgal_d_ref='Bernard-Salas et al. 2006 ApJ 652, L29',exgal_sources='SMP LMC 11',mua=0.0)
HC4N = Molecule('','HC4N',2004,'HC4N',[IRC10216],[IRAM30],['mm'],radical=True,neutral=True,H=1,C=4,N=1,d_ref='Cernicharo et al. 2004 ApJ 615, L145',lab_ref='Tang et al. 1999 CPL 315, 69',notes=None,Bcon=2302,mua=4.3)
cH2C3O = Molecule('cyclopropenone','c-H2C3O',2006,'c-H2C3O',[SgrB2],[GBT],['cm'],neutral=True,cyclic=True,H=2,C=3,O=1,d_ref='Hollis et al. 2006 ApJ 642, 933',lab_ref='Benson et al. 1973 JACS 95, 2772; Guillemin et al. 1990 JMS 140, 190',notes=None,Acon=32041,Bcon=7825,Ccon=6281,mua=4.4)
CH2CNH = Molecule('ketenimine','CH2CNH',2006,'CH2CNH',[SgrB2],[GBT],['cm'],neutral=True,H=3,C=2,N=1,d_ref='Lovas et al. 2006 ApJ 645, L137',lab_ref='Rodler et al. 1984 CPL 110, 447; Rodler et al. 1986 JMS 118, 267',notes=None,Acon=201444,Bcon=9663,Ccon=9470,mua=0.4,mub=1.4)
C5Nm = Molecule('cyanobutadiynyl anion','C5N-',2008,'C5N-',[IRC10216],[IRAM30],['mm'],anion=True,C=5,N=1,d_ref='Cernicharo et al. 2008 ApJ 688, L83',lab_ref='Botschwina & Oswald 2008 JCP 129, 044305',notes=None,Bcon=1389,mua=5.2)
HNCHCN = Molecule('E-cyanomethanimine','HNCHCN',2013,'HNCHCN',[SgrB2],[GBT],['cm'],neutral=True,H=2,C=2,N=2,d_ref='Zaleski et al. 2013 ApJ 765, L9',lab_ref='Zaleski et al. 2013 ApJ 765, L9',notes=None,Bcon=1389,mua=3.3,mub=2.5)
SiH3CN = Molecule('silyl cyanide','SiH3CN',2014,'SiH3CN',[IRC10216],[IRAM30],['mm'],neutral=True,H=3,C=1,N=1,Si=1,d_ref='Agúndez et al. 2014 A&A 570, A45',lab_ref='Priem et al. 1998 JMS 191, 183',notes='*Confirmed in Cernicharo et al. 2017 A&A 606, L5',Acon=62695,Bcon=4972,Ccon=4600,mua=3.4)
MgC4H = Molecule('magnesium butadiynyl raidcal','MgC4H',2019,'MgC4H',[IRC10216],[IRAM30,Yebes40],['cm', 'mm'],neutral=True,radical=True,Mg=1,C=4,H=1,d_ref='Cernicharo et al. 2019 A&A 630, L2',lab_ref='Forthomme et al. 2010 Chem. Phys. Lett. 488, 116',notes='Lab spectroscopy is electronic - no pure rotational spectra are available for this species.  Assignment was made based on quantum chemical calculations performed in Cernicharo et al. 2019.',Bcon=1381,mua=2.1)

#############################################################
#					Seven Atom Molecules					#
#############################################################

CH3CHO = Molecule('acetaldehyde','CH3CHO',1973,'CH3CHO',[SgrB2],[NRAO140],['cm'],neutral=True,H=4,C=2,O=1,d_ref='Gottlieb 1973 Molecules in the Galactic Environment 181; Fourikis et al. 1974 Aust J Phys 27, 425; Gilmore et al. 1976 ApJ 204, 43',lab_ref='Kilb et al. 1957 JCP 26, 1695; Souter & Wood 1970 JCP 52, 674',notes=None,ice='Tentative',ice_d_ref='Schutte et al. 1999 A&A 343, 966',ice_l_ref='Schutte et al. 1999 A&A 343, 966',exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS',Acon=56449,Bcon=10160,Ccon=9101,mua=2.4,mub=1.3)
CH3CCH = Molecule('methylacetylene','CH3CCH',1973,'CH3CCH',[SgrB2],[NRAO36],['mm'],neutral=True,H=4,C=3,d_ref='Buhl & Snyder 1973 Molecules in the Galactic Environment 187',lab_ref='Trambarulo et al. 1950 JCP 18, 1613',notes=None,exgal=True,exgal_d_ref='Mauersberger et al. 1991 A&A 247, 307',exgal_sources='NGC 253, M82',Acon=158590,Bcon=8546,Ccon=8546,mua=0.8)
CH3NH2 = Molecule('methylamine','CH3NH2',1974,'CH3NH2',[SgrB2, Orion],[Mitaka6, NRAO36, Parkes64],['cm', 'mm'],neutral=True,H=5,C=1,N=1,d_ref='Fourikis et al. 1974 ApJ 191, L139; Kaifu et al. 1974 ApJ 191, L135',lab_ref='Takagi & Kojima 1973 ApJ 181, L91',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS',Acon=103156,Bcon=22608,Ccon=21730,mua=0.3,mub=1.3)
CH2CHCN = Molecule('vinylcyanide','CH2CHCN',1975,'CH2CHCN',[SgrB2],[Parkes64],['cm'],neutral=True,H=3,C=3,N=1,d_ref='Gardner & Winnewisser 1975 ApJ 195, L127',lab_ref='Gerry & Winnewisser 1973 JMS 48, 1',notes=None,Acon=49851,Bcon=4971,Ccon=4514,mua=3.8,mub=0.9)
HC5N = Molecule('cyanodiacetylene','HC5N',1976,'HC5N',[SgrB2],[Algonquin46],['cm'],neutral=True,H=1,C=5,N=1,d_ref='Broten et al. 1976 ApJ 209, L143; Avery et al. 1976 ApJ 205 L173',lab_ref='Alexander et al. 1976 JMS 62, 175',notes=None,exgal='Tentative',exgal_d_ref='Aladro et al. 2015 A&A 579, A101',exgal_sources='NGC 253',Bcon=1331,mua=4.3)
C6H = Molecule('hexatriynyl radical','C6H',1986,'C6H',[TMC1],[Nobeyama45],['cm'],radical=True,neutral=True,H=1,C=6,d_ref='Suzuki et al. 1986 PASJ 38, 911',lab_ref='Pearson et al. 1988 A&A 189, L13',notes=None,Bcon=1391,mua=5.5)
cC2H4O = Molecule('ethylene oxide','c-C2H4O',1997,'c-C2H4O',[SgrB2],[Haystack37, Nobeyama45, SEST15],['cm', 'mm'],neutral=True,cyclic=True,H=4,C=2,O=1,d_ref='Dickens et al. 1997 ApJ 489, 753',lab_ref='Hirose 1974 ApJ 189, L145',notes=None,Acon=25484,Bcon=22121,Ccon=14098,mub=1.9)
CH2CHOH = Molecule('vinyl alcohol','CH2CHOH',2001,'CH2CHOH',[SgrB2],[NRAOARO12],['mm'],neutral=True,H=4,C=2,O=1,d_ref='Turner & Apponi 2001 ApJ 561, L207',lab_ref='Rodler 1985 JMS 114, 23; Kaushik 1977 CPL 49, 90',notes=None,Acon=62868,Bcon=10456,Ccon=8963,mua=0.5,mub=1.7)
C6Hm = Molecule('hexatriynyl anion','C6H-',2006,'C6H-',[TMC1, IRC10216],[GBT],['cm'],anion=True,H=1,C=6,d_ref='McCarthy et al. 2006 ApJ 652, L141',lab_ref='McCarthy et al. 2006 ApJ 652, L141',notes='*First gas-phase molecular anion',Bcon=1377,mua=8.2)
CH3NCO = Molecule('methyl isocyanate','CH3NCO',2015,'CH3NCO',[SgrB2, Orion],[NRAOARO12, SMT10],['mm'],neutral=True,H=3,C=2,O=1,N=1,d_ref='Halfen et al. 2015 ApJ 812, L5',lab_ref='Halfen et al. 2015 ApJ 812, L5',notes='*see also Cernicharo et al. 2016 A&A 587, L4',Acon=128400,Bcon=4415,Ccon=4257,mua=2.9)
HC5O = Molecule('butadiynylformyl radical','HC5O',2017,'HC5O',[TMC1],[GBT],['cm'],radical=True,neutral=True,H=1,C=5,O=1,d_ref='McGuire et al. 2017 ApJ 843, L28',lab_ref='Mohamed et al. 2005 JCP 123, 234301',notes=None,Bcon=1294,mua=2.2)
HOCH2CN = Molecule('glycolonitrile','HOCH2CN',2019,'HOCH2CN',[IRAS16293],[ALMA],['mm'],neutral=True,H=3,C=2,N=1,O=1,d_ref='Zeng et al. 2019 MNRAS 484, L43',lab_ref='Margules et al. 2017 A&A 601, A50',notes=None,Acon=33610,Bcon=4838,Ccon=4377,mua=2.32,mub=1.31,muc=1.23)
HC4NC = Molecule('isocyanoacetylene',
					'HC4NC',
					2020,
					'HC4NC',
					[TMC1],
					[GBT],
					['mm'],
					neutral=True,
					H=1,
					C=5,
					N=1,
					d_ref='Xue et al. 2020 ApJL 900, L9',
					lab_ref='Botschwina et al. 1998 JCP 109, 3108',
					notes='Also known as isocyanobutadiyne',
					Bcon=1402,
					mua=3.24)

HC3HNH = Molecule(name = 'propargylamine',
					formula = 'HC3HNH',
					year = 2020,
					label = 'HC3HNH',
					sources = [G0693],
					telescopes = [IRAM30],
					wavelengths = ['mm'],
					neutral=True,
					H=3,
					C=3,
					N=1,
					d_ref='Bizzocchi et al. 2020 A&A 640, A98',
					lab_ref='Bizzocchi et al. 2020 A&A 640, A98; Kroto et al. 1984 J. Chem. Soc. Chem. Comm. 993; Sugie et al. 1985 JMS 111, 83; McNaughton et al. 1988 J. Mol. Struct. 190, 195.',
					Acon=54640,
					Bcon=4862,
					Ccon=4458,
					mua=2.14,
					mub=0.17,)



#############################################################
#					Eight Atom Molecules					#
#############################################################

HCOOCH3 = Molecule('methyl formate','HCOOCH3',1975,'HCOOCH3',[SgrB2],[Parkes64,Effelsberg100],['cm'],neutral=True,H=4,C=2,O=2,d_ref='Churchwell & Winnewisser 1975 A&A 45, 229; Brown et al. 1975 ApJ 197, L29',lab_ref='Brown et al. 1975 ApJ 197, L29',notes='*t-mf detected 2012 ApJ 755, 143',exgal=True,exgal_d_ref='Sewiło et al. 2018 ApJL 853, L19',exgal_sources='LMC',Acon=17630,Bcon=9243,Ccon=5318,mua=1.6,mub=0.7)
CH3C3N = Molecule('methylcyanoacetylene','CH3C3N',1984,'CH3C3N',[TMC1],[NRAO140],['cm'],neutral=True,H=3,C=4,N=1,d_ref='Broten et al. 1984 ApJ 276, L25',lab_ref='Moises et al. 1982 JMS 92, 497',notes=None,Acon=158099,Bcon=2066,Ccon=2066,mua=4.8)
C7H = Molecule('heptatriynylidyne radical','C7H',1997,'C7H',[IRC10216],[IRAM30],['mm'],radical=True,neutral=True,H=1,C=7,d_ref='Guélin et al. 1997 A&A 317, L1',lab_ref='Travers et al. 1996 ApJ 465, L77',notes=None,Bcon=875,mua=5.9)
CH3COOH = Molecule('acetic acid','CH3COOH',1997,'CH3COOH',[SgrB2],[BIMA, OVRO],['mm'],neutral=True,H=4,C=2,O=2,d_ref='Mehringer et al. 1997 ApJ 480, L71',lab_ref='Tabor 1957 JCP 27, 974',notes=None,Acon=11335,Bcon=9479,Ccon=5325,mua=2.9,mub=4.9)
H2C6 = Molecule('hexapentaenylidene','H2C6',1997,'H2C6',[TMC1],[Goldstone70],['cm'],neutral=True,H=2,C=6,d_ref='Langer et al. 1997 ApJ 480, L63',lab_ref='McCarthy et al. 1997 Science 275, 518',notes=None,Acon=268400,Bcon=1348,Ccon=1341,mua=6.2)
CH2OHCHO = Molecule('glycolaldehyde','CH2OHCHO',2000,'CH2OHCHO',[SgrB2],[NRAOARO12],['mm'],neutral=True,H=4,C=2,O=2,d_ref='Hollis et al. 2000 ApJ 540, L107',lab_ref='Marstokk & Mollendal 1973 J Mol Struct 16, 259',notes=None,Acon=18446,Bcon=6526,Ccon=4969,mua=0.3,mub=2.3)
HC6H = Molecule('triacetylene','HC6H',2001,'HC6H',[CRL618],[ISO],['IR'],neutral=True,H=2,C=6,d_ref='Cernicharo et al. 2001 ApJ 546, L123',lab_ref='Haas etal. 1994 JMS 167, 176',notes=None,exgal=True,exgal_d_ref='Bernard-Salas et al. 2006 ApJ 652, L29',exgal_sources='SMP LMC 11',mua=0.0)
CH2CHCHO = Molecule('propenal','CH2CHCHO',2004,'CH2CHCHO',[SgrB2, G327306LOS],[NRAOARO12, SEST15],['cm'],neutral=True,H=4,C=3,O=1,d_ref='Hollis et al. 2004 ApJ 610, L21',lab_ref='Winnewisser et al. 1975 Z Naturforsch 30, 1001',notes=None,Acon=47354,Bcon=4660,Ccon=4243,mua=3.1,mub=0.6)
CH2CCHCN = Molecule('cyanoallene','CH2CCHCN',2006,'CH2CCHCN',[TMC1],[GBT],['cm'],neutral=True,H=3,C=4,N=1,d_ref='Lovas et al. 2006 ApJ 637, L37',lab_ref='Bouche et al. 1973 J Mol Struct 18, 211',notes='*Also Chin et al. 2006 AIP Conf. Proc. 855, 149',Acon=25981,Bcon=2689,Ccon=2475,mua=4.1,mub=1.3)
NH2CH2CN = Molecule('aminoacetonitrile','NH2CH2CN',2008,'NH2CH2CN',[SgrB2],[IRAM30, PdBI, ATCA],['mm'],neutral=True,H=4,C=2,N=2,d_ref='Belloche et al. 2008 A&A 482, 179',lab_ref='Bogey et al. 1990 JMS 143, 180',notes=None,Acon=30246,Bcon=4761,Ccon=4311,mua=2.6,mub=0.6)
CH3CHNH = Molecule('ethanimine','CH3CHNH',2013,'CH3CHNH',[SgrB2],[GBT],['cm'],neutral=True,H=5,C=2,N=1,d_ref='Loomis et al. 2013 ApJL 765, L10',lab_ref='Loomis et al. 2013 ApJL 765, L10',notes=None,Acon=49961,Bcon=9828,Ccon=8650,mua=0.8,mub=1.9)
CH3SiH3 = Molecule('methyl silane','CH3SiH3',2017,'CH3SiH3',[IRC10216],[IRAM30],['mm'],neutral=True,H=6,C=1,Si=1,d_ref='Cernicharo et al. 2017 A&A 606, L5',lab_ref='Wong et al. 1983 JMS 102, 89',notes=None,Acon=56189,Bcon=10986,Ccon=10986,mua=0.7)
NH2CONH2 = Molecule('urea','NH2CONH2',2019,'NH2CONH2',[SgrB2],[ALMA],['mm'],neutral=True,H=4,N=2,C=1,O=1,d_ref='Belloche et al. 2019 A&A 628, A10',lab_ref='Brown et al. 1975 JMS 58, 445; Kasten & Dreizler 1986 Z. Naturforsch A. 41, 1173; Kretschmer et al. 1996 Mol. Phys. 87, 1159; Godfrey et al. 1997 J. Mol. Struct. 413-414, 405; Remijan et al. 2014 ApJ 783, 77; Additional work used in Belloche et al. 2019 A&A 628, A10 to be reported in Medvedev et al. in prep as of 9/16/2019.',notes='Evidence for the detection, but no claim, made in Remijan et al. 2014 ApJ 783, 77.  Dipole is from Brown et al. 1975 JMS 58, 445',Acon=11233,Bcon=10369,Ccon=5417,mub=3.83)
HCCCH2CN = Molecule(name ='propargyl cyanide',
						formula='HCCCH2CN',
						year=2020,
						label='HCCCH2CN',
						sources=[TMC1],
						telescopes=[GBT],
						wavelengths=['cm'],
						neutral=True,
						H=3,
						C=4,
						N=1,
						d_ref='McGuire et al. 2020 ApJL 900, L10',
						lab_ref='Jones & Sheridan 1982 J Mol Struct 78, 303; Demaison et al. 1985 JMS 114, 210; McNaughton et al. 1988 JMS 132, 407; Jager et al. 1990 JMS 143, 50; McGuire et al. 2020 ApJL 900, L10',
						notes='Also known as 3-butynenitrile and 1-cyanoprop-2-yne',
						Acon=19820,
						Bcon=2910,
						Ccon=2573,
						mua=3.23,
						mub=2.34)


#############################################################
#					Nine Atom Molecules						#
#############################################################

CH3OCH3 = Molecule('dimethyl ether','CH3OCH3',1974,'CH3OCH3',[Orion],[NRAO36, NRL85],['cm', 'mm'],neutral=True,H=6,C=2,O=1,d_ref='Snyder et al. 1974 ApJ 191, L79',lab_ref='Kasai & Myers JCP 30, 1096; Blukis et al. 1963 JCP 38, 2753',notes=None,exgal=True,exgal_d_ref='Qiu et al. 2018 A&A 613, A3; Sewiło et al. 2018 ApJL 853, L19',exgal_sources='NGC 1068, LMC',Acon=38788,Bcon=10057,Ccon=8887,mub=1.3)
CH3CH2OH = Molecule('ethanol','CH3CH2OH',1975,'CH3CH2OH',[SgrB2],[NRAO36],['mm'],neutral=True,H=6,C=2,O=1,d_ref='Zukerman et al. 1975 ApJ 196, L99',lab_ref='Takano et al. 1986 JMS 26, 157',notes='*g-ethanol detected 1997 ApJ 480, 420',Acon=34892,Bcon=9351,Ccon=8135,mua=0.1,mub=1.4)
CH3CH2CN = Molecule('ethyl cyanide','CH3CH2CN',1977,'CH3CH2CN',[SgrB2, Orion],[NRAO36],['mm'],neutral=True,H=5,C=3,N=1,d_ref='Johnson et al. 1977 ApJ 218, 370',lab_ref='Johnson et al. 1977 ApJ 218, 370',notes=None,Acon=27664,Bcon=4714,Ccon=4235,mua=3.9,mub=1.2)
HC7N = Molecule('cyanotriacetylene','HC7N',1977,'HC7N',[TMC1],[Algonquin46, Haystack37],['cm'],neutral=True,H=1,C=7,N=1,d_ref='Kroto et al. 1977 Bull. Am. As. Soc. 9, 303',lab_ref='Kirby et al. 1980 JMS 83, 261',notes=None,Bcon=564,mua=4.8)
CH3C4H = Molecule('methyldiacetylene','CH3C4H',1984,'CH3C4H',[TMC1],[Haystack37, NRAO140, Effelsberg100],['cm'],neutral=True,H=4,C=5,d_ref='Walmsley et al. 1984 A&A 134, L11',lab_ref='Heath et al. 1955 Faraday Discuss. 19, 38',notes=None,Acon=159140,Bcon=2036,Ccon=2036,mua=1.2)
C8H = Molecule('octatriynyl radical','C8H',1996,'C8H',[IRC10216],[IRAM30, Nobeyama45],['mm'],radical=True,neutral=True,H=1,C=8,d_ref='Cernicharo & Guélin 1996 A&A 309, L27',lab_ref='Pauzat et al. 1991 ApJ 369, L13',notes=None,Bcon=587,mua=6.5)
CH3CONH2 = Molecule('acetamide','CH3CONH2',2006,'CH3CONH2',[SgrB2],[GBT],['cm'],neutral=True,H=5,C=2,O=1,N=1,d_ref='Hollis et al. 2006 ApJ 643, L25',lab_ref='Suenram et al. 2001 JMS 208, 188',notes=None,Acon=10788,Bcon=9331,Ccon=5157,mua=1.1,mub=3.5)
C8Hm = Molecule('octatriynyl anion','C8H-',2007,'C8H-',[TMC1, IRC10216],[GBT],['cm'],anion=True,H=1,C=8,d_ref='Brünken et al. 2007 ApJ 664, L43; Remijan et al. 2007 ApJ 664, L47',lab_ref='Gupta et al. 2007 ApJ 655, L57',notes=None,Bcon=583,mua=10.4)
CH2CHCH3 = Molecule('propylene','CH2CHCH3',2007,'CH2CHCH3',[TMC1],[IRAM30],['mm'],neutral=True,H=6,C=3,d_ref='Marcelino et al. 2007 ApJ 665, L127',lab_ref='Pearson et al. 1994 JMS 166, 120; Wlodarczak et al. 1994 JMS 167, 239',notes=None,Acon=46281,Bcon=9308,Ccon=8130,mua=0.4,mub=0.1)
CH3CH2SH = Molecule('ethyl mercaptan','CH3CH2SH',2014,'CH3CH2SH',[Orion],[IRAM30],['mm'],neutral=True,H=6,C=2,S=1,d_ref='Kolesniková et al. 2014 ApJ 784, L7',lab_ref='Kolesniková et al. 2014 ApJ 784, L7',notes=None,Acon=28747,Bcon=5295,Ccon=4846,mua=1.5,mub=0.2,muc=0.6)
HC7O = Molecule('hexadiynylformyl radical','HC7O',2017,'HC7O',[TMC1],[GBT],['cm'],radical=True,neutral=True,H=1,C=7,O=1,d_ref='McGuire et al. 2017 ApJ 843, L28',lab_ref='Mohamed et al. 2005 JCP 123, 234301',notes='*Confirmed in Cordiner et al. 2017 ApJ 850, 194',Bcon=549,mua=2.2)

#############################################################
#					Ten Atom Molecules						#
#############################################################

acetone = Molecule('acetone','(CH3)2CO',1987,'acetone',[SgrB2],[IRAM30, NRAO140, NRAOARO12],['cm', 'mm'],neutral=True,H=6,C=3,O=1,d_ref='Combes et al. 1987 A&A 180, L13',lab_ref='Vacherand et al. 1986 JMS 118, 355',notes='*Confirmed in 2002 ApJ 578, 245',Acon=10165,Bcon=8515,Ccon=4910,mub=2.9)
HOCH2CH2OH = Molecule('ethylene glycol','HOCH2CH2OH',2002,'HOCH2CH2OH',[SgrB2],[NRAOARO12],['mm'],neutral=True,H=6,C=2,O=2,d_ref='Hollis et al. 2002 ApJ 571, L59',lab_ref='Christen et al. 1995 JMS 172, 57',notes='*aGg\' conformer in 2017 A&A 598, A59',Acon=15361,Bcon=5588,Ccon=4614,mua=2.1,mub=0.9)
CH3CH2CHO = Molecule('propanal','CH3CH2CHO',2004,'CH3CH2CHO',[SgrB2],[GBT],['cm'],neutral=True,H=6,C=3,O=1,d_ref='Hollis et al. 2004 ApJ 610, L21',lab_ref='Butcher & Wilson 1964 JCP 40, 1671',notes=None,Acon=16712,Bcon=5969,Ccon=4648,mua=1.7,mub=1.9)
CH3C5N = Molecule('methylcyanodiacetylene','CH3C5N',2006,'CH3C5N',[TMC1],[GBT],['cm'],neutral=True,H=4,C=6,N=1,d_ref='Snyder et al. 2006 ApJ 647, 412',lab_ref='Chen et al. 1998 JMS 192, 1',notes=None,Acon=158099,Bcon=778,Ccon=778,mua=5.4)
CH3CHCH2O = Molecule('propylene oxide','CH3CHCH2O',2016,'CH3CHCH2O',[SgrB2],[GBT],['cm'],neutral=True,cyclic=True,H=6,C=3,O=1,d_ref='McGuire & Carroll et al. 2016 Science 352, 1449',lab_ref='McGuire & Carroll et al. 2016 Science 352, 1449',notes='*First chiral molecule',Acon=18024,Bcon=6682,Ccon=5951,mua=1.0,mub=1.7,muc=0.6)
CH3OCH2OH = Molecule('methoxymethanol','CH3OCH2OH',2017,'CH3OCH2OH',[NGC6334],[ALMA],['mm'],neutral=True,H=6,C=2,O=2,d_ref='McGuire et al. 2017 ApJ 851, L46',lab_ref='Motiyenko et al. 2018 PCCP 20, 5509',notes=None,Acon=17238,Bcon=5568,Ccon=4813,mua=0.2,mub=0.1,muc=0.1)

#############################################################
#					Eleven Atom Molecules					#
#############################################################

HC9N = Molecule('cyanotetraacetylene','HC9N',1978,'HC9N',[TMC1],[Algonquin46, NRAO140],['cm'],neutral=True,H=1,C=9,N=1,d_ref='Broten et al. 1978 ApJ 223, L105',lab_ref='Iida et al. 1991 ApJ 371, L45',notes=None,Bcon=291,mua=5.2)
CH3C6H = Molecule('methyltriacetylene','CH3C6H',2006,'CH3C6H',[TMC1],[GBT],['cm'],neutral=True,H=4,C=7,d_ref='Remijan et al. 2006 ApJ 643, L37',lab_ref='Alexander et al. 1978 JMS 70, 84',notes=None,Acon=159140,Bcon=778,Ccon=778,mua=1.5)
C2H5OCHO = Molecule('ethyl formate','C2H5OCHO',2009,'C2H5OCHO',[SgrB2],[IRAM30],['mm'],neutral=True,H=6,C=3,O=2,d_ref='Belloche et al. 2009 A&A 499, 215',lab_ref='Medvedev et al. 2009 ApJS 181, 433',notes=None,Acon=17747,Bcon=2905,Ccon=2579,mua=1.9,mub=0.7,muc=0.0)
CH3COOCH3 = Molecule('methyl acetate','CH3COOCH3',2013,'CH3COOCH3',[Orion],[IRAM30],['mm'],neutral=True,H=6,C=3,O=2,d_ref='Tercero et al. 2013 ApJ 770, L13',lab_ref='Tudorie et al. 2011 JMS 269, 211',notes=None,Acon=10247,Bcon=4170,Ccon=3077,mua=0.0,mub=1.6)
CH3COCH2OH = Molecule('hydroxyacetone','CH3COCH2OH',2020,'CH3COCH2OH',[IRAS16293],[ALMA],['mm'],neutral=True,H=6,C=3,O=2,
	d_ref='Zhou et al. 2020 Res. Astron. & Astrophys. 20, 125',
	lab_ref='Kattija-Ari & Harmony et al. 1980 Int. J. Quant. Chem.: Quant. Chem Symp. 14, 18, 443; Apponi et al. 2006 ApJ 652, 1787; Braakman et al. 2010 JMS 264, 43',
	notes=None,Acon=10074,Bcon=3817,Ccon=2867,mua=2.22,mub=2.17)

#############################################################
#					Twelve Atom Molecules					#
#############################################################

C6H6 = Molecule('benzene','C6H6',2001,'C6H6',[CRL618],[ISO],['IR'],neutral=True,cyclic=True,H=6,C=6,d_ref='Cernicharo et al. 2001 ApJ 546, L123',lab_ref='Lindenmayer et al. 1988 JMS 128 172',notes=None,exgal=True,exgal_d_ref='Bernard-Salas et al. 2006 ApJ 652, L29',exgal_sources='SMP LMC 11',mua=0.0)
nC3H7CN = Molecule('n-propyl cyanide','n-C3H7CN',2009,'n-C3H7CN',[SgrB2],[IRAM30],['mm'],neutral=True,H=7,C=4,N=1,d_ref='Belloche et al. 2009 A&A 499, 215',lab_ref='Belloche et al. 2009 A&A 499, 215',notes=None,Acon=23668,Bcon=2268,Ccon=2153,mua=4.0,mub=1.0,muc=0.0)
iC3H7CN = Molecule('isopropyl cyanide','i-C3H7CN',2014,'i-C3H7CN',[SgrB2],[ALMA],['mm'],neutral=True,H=7,C=4,N=1,d_ref='Belloche et al. 2014 Science 345, 1584',lab_ref='Muller et al. 2011 JMS 267, 100',notes=None,Acon=7941,Bcon=3968,Ccon=2901,mua=4.0,mub=0.6)
C5H5CN1 = Molecule('1-cyano-1,3-cyclopentadiene','C5H5CN',2020,'C5H5CN',[TMC1],[GBT],['cm'],neutral=True,cyclic=True,C=6,H=5,N=1,d_ref='McCarthy et al. 2020 Nature Astronomy, doi:10.1038/s41550-020-01213-y.',lab_ref='McCarthy et al. 2020 Nature Astronomy, doi:10.1038/s41550-020-01213-y.',notes=None,Acon=8353,Bcon=1904,Ccon=1565,mua=4.15)
#C5H5CN2 = Molecule('2-cyano-1,3-cyclopentadiene','C5H5CN',2020,'C5H5CN',[TMC1],[GBT],['cm'],neutral=True,cyclic=True,C=6,H=5,N=1,d_ref='Lee et al. 2021 ApJL, in press',lab_ref='McCarthy et al. 2020 Nature Astronomy, doi:10.1038/s41550-020-01213-y, Lee et al. 2021 ApJL, in press',notes=None,Acon=8236,Bcon=1902,Ccon=1560,mua=4.36)


#############################################################
#					Thirteen Atom Molecules					#
#############################################################

cC6H5CN = Molecule('benzonitrile','C6H5CN',2018,'C6H5CN',[TMC1],[GBT],['cm'],neutral=True,cyclic=True,H=5,C=7,N=1,d_ref='McGuire et al. 2018 Science 359, 202',lab_ref='Wohlfart et al. 2008 JMS 247, 119',notes=None,Acon=5655,Bcon=1547,Ccon=1214,mua=4.5)
HC11N = Molecule('cyanopentaacetylene','HC11N',2020,'HC11N',[TMC1],[GBT],['cm'],neutral=True,H=1,C=11,N=1,d_ref='',lab_ref='Travers et al. 1996 ApJL 469, L65',notes=None,Bcon=169,mua=5.47)

#############################################################
#						PAH Molecules						#
#############################################################

CNN1 = Molecule(name='1-cyanonaphthalene',
				formula='C10H7CN',
				year=2020,
				label='CNN1',
				sources=[TMC1],
				telescopes=[GBT],
				wavelengths=['cm'],
				neutral=True,
				cyclic=True,
				pah=True,
				C=11,
				H=7,
				N=1,
				d_ref='', 
				lab_ref='McNaughton et al. 2018 MNRAS 476, 5268', 
				notes=None,
				Acon=1479,
				Bcon=957,
				Ccon=581,
				mua=3.56,
				mub=2.96)
CNN2 = Molecule('2-cyanonaphthalene','C10H7CN',2020,'CNN2',[TMC1],[GBT],['cm'],neutral=True,cyclic=True,pah=True,C=11,H=7,N=1,d_ref='', lab_ref='McNaughton et al. 2018 MNRAS 476, 5268', notes=None, Acon=2707,Bcon=606,Ccon=495,mua=5.09,mub=0.98)

#############################################################
#					Fullerene Molecules						#
#############################################################

C60 = Molecule('buckminsterfullerene','C60',2010,'C60',[TC1, NGC7023],[Spitzer],['IR'],neutral=True,fullerene=True,cyclic=True,C=60,d_ref='Cami et al. 2010 Science 329, 1180',lab_ref='Nemes et al. 1994 CPL 218, 295',notes='*See also Sellgren et al. 2010 ApJ 722, L54 and Werner 2004b, Sellgren 2007 therein',mua=0.0)
C60p = Molecule('buckminsterfullerene cation','C60+',2013,'C60+',[NGC7023],[Spitzer],['IR'],cation=True,fullerene=True,cyclic=True,C=60,d_ref='Berné et al. 2013 A&A 550, L4',lab_ref='Kern et al. 2013 JPCA 117, 8251',notes='*See also Campbell et al. 2015 Nature 523, 322',mua=0.0)
C70 = Molecule('rugbyballene','C70',2010,'C70',[TC1],[Spitzer],['IR'],neutral=True,fullerene=True,cyclic=True,C=70,d_ref='Cami et al. 2010 Science 329, 1180',lab_ref='Nemes et al. 1994 CPL 218, 295',notes=None,mua=0.0)

#############################################################
#							Full List						#
#############################################################

full_list = [
	CH,
	CN,
	CHp,
	OH,
	CO,
	H2,
	SiO,
	CS,
	SO,
	SiS,
	NS,
	C2,
	NO,
	HCl,
	NaCl,
	AlCl,
	KCl,
	AlF,
	PN,
	SiC,
	CP,
	NH,
	SiN,
	SOp,
	COp,
	HF,
	N2,
	CFp,
	PO,
	O2,
	AlO,
	CNm,
	OHp,
	SHp,
	HClp,
	SH,
	TiO,
	ArHp,
	NSp,
	HeHp,
	VO,
	H2O,
	HCOp,
	HCN,
	OCS,
	HNC,
	H2S,
	N2Hp,
	C2H,
	SO2,
	HCO,
	HNO,
	HCSp,
	HOCp,
	SiC2,
	C2S,
	C3,
	CO2,
	CH2,
	C2O,
	MgNC,
	NH2,
	NaCN,
	N2O,
	MgCN,
	H3p,
	SiCN,
	AlNC,
	SiNC,
	HCP,
	CCP,
	AlOH,
	H2Op,
	H2Clp,
	KCN,
	FeCN,
	HO2,
	TiO2,
	CCN,
	SiCSi,
	S2H,
	HCS,
	HSC,
	NCO,
	CaNC,
	NH3,
	H2CO,
	HNCO,
	H2CS,
	C2H2,
	C3N,
	HNCS,
	HOCOp,
	C3O,
	lC3H,
	HCNHp,
	H3Op,
	C3S,
	cC3H,
	HC2N,
	H2CN,
	SiC3,
	CH3,
	C3Nm,
	PH3,
	HCNO,
	HOCN,
	HSCN,
	HOOH,
	lC3Hp,
	HMgNC,
	HCCO,
	CNCN,
	HONO,
	HC3N,
	HCOOH,
	CH2NH,
	NH2CN,
	H2CCO,
	C4H,
	SiH4,
	cC3H2,
	CH2CN,
	C5,
	SiC4,
	H2CCC,
	CH4,
	HCCNC,
	HNCCC,
	H2COHp,
	C4Hm,
	CNCHO,
	HNCNH,
	CH3O,
	NH3Dp,
	H2NCOp,
	NCCNHp,
	CH3Cl,
	CH3OH,
	CH3CN,
	NH2CHO,
	CH3SH,
	C2H4,
	C5H,
	CH3NC,
	HC2CHO,
	H2C4,
	C5S,
	HC3NHp,
	C5N,
	HC4H,
	HC4N,
	cH2C3O,
	CH2CNH,
	C5Nm,
	HNCHCN,
	SiH3CN,
	CH3CHO,
	CH3CCH,
	CH3NH2,
	CH2CHCN,
	HC5N,
	C6H,
	cC2H4O,
	CH2CHOH,
	C6Hm,
	CH3NCO,
	HC5O,
	HOCH2CN,
	HCOOCH3,
	CH3C3N,
	C7H,
	CH3COOH,
	H2C6,
	CH2OHCHO,
	HC6H,
	CH2CHCHO,
	CH2CCHCN,
	NH2CH2CN,
	CH3CHNH,
	CH3SiH3,
	NH2CONH2,
	CH3OCH3,
	CH3CH2OH,
	CH3CH2CN,
	HC7N,
	CH3C4H,
	C8H,
	CH3CONH2,
	C8Hm,
	CH2CHCH3,
	CH3CH2SH,
	HC7O,
	acetone,
	HOCH2CH2OH,
	CH3CH2CHO,
	CH3C5N,
	CH3CHCH2O,
	CH3OCH2OH,
	HC9N,
	CH3C6H,
	C2H5OCHO,
	CH3COOCH3,
	C6H6,
	nC3H7CN,
	iC3H7CN,
	cC6H5CN,
	C60,
	C60p,
	C70,
	CNN1,
	CNN2,
	HCCCH2CN,
	HC4NC,
	HC11N,
	C5H5CN1,
	#C5H5CN2,
	MgC3N,
	MgC4H,
	HC3HNH,
	HC3Op
	]


#############################################################
#				     Do Some Loop Updates	 				#
#############################################################

for x in scopes_list:

	x.update_stats(full_list)
	
for x in source_list:

	x.update_stats(full_list)	

#############################################################
#						Functions	 						#
#############################################################	

def update_plots(list):

	'''
	A meta function that, when run, will call every plot command and generate new plots based on the input list using default parameters.  Useful for rapidly re-generating all figures.
	'''
	cumu_det_plot(full_list)
	cumu_det_natoms_plot(full_list)
	det_per_year_per_atoms(full_list)
	facility_shares(scopes_list, full_list)
	cumu_det_facility(full_list)
	periodic_heatmap(full_list)
	mass_by_wavelengths(full_list)
	mols_waves_by_atoms(full_list)
	du_histogram(full_list)
	type_pie_chart(full_list)
	source_pie_chart(full_list)
	mol_type_by_source_type(full_list)
	du_by_source_type(full_list)
	rel_du_by_source_type(full_list)
	mass_by_source_type(full_list)
	waves_by_source_type(full_list)
	
	return

def summary(y):

	'''
	Prints a summary of the information in the database for either a molecule or a source to the terminal.  Requires the variable name (usually intuitive, but if not, do a plain text search...).
	
	'''
	
	if isinstance(y,Molecule):
	
		n_dash = len(y.name) + len(y.formula) +3 
	
		dashes = '-' * n_dash

		print('\n' + dashes)
		print('{} ({})' .format(y.name,y.formula))
		print(dashes + '\n')
		print('Atoms:\t{}' .format(y.natoms))
		print('Mass:\t{} amu' .format(y.mass))
		print('Year Detected:\t{}' .format(y.year))
		sources = [x.name for x in y.sources]
		sources_str = ', '.join(sources)		
		print('Source(s):\t{}' .format(sources_str))
		scopes = [x.shortname for x in y.telescopes]
		scopes_str = ', '.join(scopes)
		print('Telescope(s) Used:\t{}' .format(scopes_str))
	
		attr_str = ''
	
		if y.neutral == True:
		
			attr_str += 'Neutral, '
		
		if y.cation == True:
	
			attr_str += 'Cation, '
		
		if y.anion == True:
	
			attr_str += 'Anion, '
		
		if y.cyclic == True:
	
			attr_str += 'Cyclic'
		
		if y.radical == True:
	
			attr_str += 'Radical, '
	
		attr_str = attr_str.strip(' ').strip(',')
	
		print('Attributes:\t{}\n' .format(attr_str))
		
		if y.isos != None:
		
			print('Known Isotopologues:\t{}\n' .format(y.isos))
		
		other_envs = [y.ice,y.ppd,y.exgal,y.exo]
		
		if any(other_envs) == True:
		
			other_str = ''
			
			if y.ice == True:
			
				other_str += 'Ices, '
				
			if y.ice == 'Tentative':
			
				other_str += 'Ices (Tentative), '	
				
			if y.ppd == True:
			
				other_str += 'Protoplanetary Disks, '
				
			if y.ppd == 'Tentative':
			
				other_str += 'Protoplanetary Disks (Tentative), '				
				
			if y.exgal == True:
			
				other_str += 'External Galaxies, '
				
			if y.exgal == 'Tentative':
			
				other_str += 'External Galaxies (Tentative), '				
				
			if y.exo == True:
			
				other_str += 'Exoplanetary Atmospheres, '
				
			if y.exo == 'Tentative':
			
				other_str += 'Exoplanetary Atmospheres (Tentative), '				
				
			other_str = other_str.strip().strip(',') + '\n'
		
			print('Also Detected In:\t{}' .format(other_str))
			
		if y.exgal == True or y.exgal == 'Tentative':
		
			print('Sources of External Galaxy Detections:\t{}\n' .format(y.exgal_sources))	

		if y.ppd_isos != None:
		
			print('Isotopologues Also Detected in Protoplanetary Disks:\t{}\n' .format(y.ppd_isos))	
			

	
		refs(y)
		
	elif isinstance(y,Source):
	
		n_dash = len(y.name)
		
		dashes = '-' * n_dash
		
		print('\n' + dashes)
		print('{}' .format(y.name))
		print(dashes + '\n')

		print('RA (J2000):\t{}' .format(y.ra))
		print('DEC (J2000):\t{}\n' .format(y.dec))
		
		print('Generalized Type:\t{}\n' .format(y.type))
		
		print('Number of Detections:\t{}\n' .format(y.detects))
		
		print('Simbad URL:\t{}' .format(y.simbad_url))
			
		mol_str = 'Molecules Detected in {}' .format(y.name)
		
		dashes = '-' * len(mol_str)
			
		print('\n' + dashes)	
		print(mol_str)
		print(dashes)
		
		print('{}' .format(mols_in_source(y.name)).replace("'",'').strip(']').strip('['))
	
def refs(y):

	'''
	
	Prints a nicely formatted list of references and notes for a molecule to the terminal.
	
	'''

	lab_refs = y.lab_ref.split(';')
	d_refs = y.d_ref.split(';')
	
	if y.notes != None:
		notes = y.notes.strip('*')
	
	print('Detection Reference(s)')
	
	for x in range(len(d_refs)):
	
		print('[{}] {}' .format(x+1,d_refs[x].strip()))
		
	print('\nLaboratory Reference(s)')	
	
	for x in range(len(lab_refs)):
	
		print('[{}] {}' .format(x+1,lab_refs[x].strip()))	
		
	if y.notes!= None:
		print('\nNotes')
		print(notes)
		
	if y.isos != None:
	
		iso_d_refs = y.isos_d_ref.split('[')
		
		#Not implemented yet
		#iso_l_refs = y.isos_l_ref.split('[')
		
		del iso_d_refs[0]
		#del iso_l_refs[0]
	
		print('\nIsotopologue Detection Reference(s)')
		
		for x in iso_d_refs:
		
			print('[' + x.strip())
		
		#print('\nIsotopologue Laboratory Reference(s)')	
			
		#for x in iso_l_refs:
		
			#print('[' + x.strip())	

	if y.ice == True or y.ice == 'Tentative':
	
		print('\nIce Reference(s)')
		
		print('[Det] {}' .format(y.ice_d_ref))
		print('[Lab] {}' .format(y.ice_l_ref))

	if y.ppd == True or y.ppd == 'Tentative':
	
		print('\nProtoplanetary Disks Reference(s)')
		
		print('[{}] {}' .format(y.formula,y.ppd_d_ref))
		
		
		#to be enabled later if lab references are cataloged.
		#print('[Lab] {}' .format(y.ppd_l_ref))	
		
		if y.ppd_isos != None:
		
			ppd_isos_refs = y.ppd_isos_ref.split('[')
			del ppd_isos_refs[0]
			
			for x in ppd_isos_refs:
			
				print('[' + x.strip())
			
	if y.exgal == True or y.exgal == 'Tentative':
	
		print('\nExternal Galaxies Reference(s)')
		
		print('[{}] {}' .format(y.formula,y.exgal_d_ref))
		
		#to be enabled later if lab references are cataloged.
		#print('[Lab] {}' .format(y.exgal_l_ref))		
		
	if y.exo == True or y.exo == 'Tentative':
	
		print('\nExoplanetary Atmospheres Reference(s)')
		
		print('[{}] {}' .format(y.formula,y.exo_d_ref))
		
def output_summary(y,filename=None):

	'''
	Writes out an ascii file containing the output of summary(y) for molecule y or a list of molecules y.  A filename can optionally be specified.
	
	'''

	#if we're using a default filename and only one molecule has been specified, name the output file that molecules formula

	if filename == None and type(y) != list:
	
		filename = y.formula + '.txt'
		
	elif filename == None and type(y) == list:
	
		filename = 'summary.txt'
	
	if type(y) != list:
	
		y = [y]
		
	with open(filename,'w') as output:
	
		for x in y:
		
			out_list = summary_list(x)
			
			for z in out_list:
			
				output.write(z + '\n')	

def cumu_det_plot(list,syear=None,eyear=None):

	'''
	Makes a plot of the cumulative detections by year using 'list', which is usually 'full_list'.  The start year and end year are defined by default to be the earliest year in the list and the current year, but these are overridable. 
	
	'''
	
	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Cumulative Detections')
	
	fig = plt.figure(num='Cumulative Detections',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')	
	
	#get the starting and ending years, if they aren't set by the user
	
	if syear is None:
	
		#grab the earliest year in the list of molecules
	
		syear = min([x.year for x in list])
		
	if eyear is None:
	
		#grab the current year
	
		eyear = date.today().year
		
	#make the x-axis array of years
	
	years = np.arange(syear,eyear+1)
		
	#make an array to hold the detections
	
	dets = np.copy(years)*0
	
	#loop through the years and the list and add everything up.  There's gotta be a better way to do this, but who cares, its fast enough.
	
	for x in range(len(dets)):
	
		i = 0
		
		for mol in list:
		
			if mol.year < years[x]+1:
			
				i += 1
				
		dets[x] = i
		
	#get some year indicies for years we care about
	
	def iyear(x):
	
		return np.argwhere(years==x)[0][0]

	
	#do linear fits to the data for the two ranges we care about (1968-Present, 2005-Present)
	#get the slope of the fit to detections since 1968 and since 2005
	
	trend1968 = np.polynomial.polynomial.Polynomial.fit(years[iyear(1968):],dets[iyear(1968):],1).convert().coef[1]
	trend2005 = np.polynomial.polynomial.Polynomial.fit(years[iyear(2005):],dets[iyear(2005):],1).convert().coef[1]

	#load up an axis
	
	ax = fig.add_subplot(111)
	
	#label the axes
	
	plt.xlabel('Year')
	plt.ylabel('Cumulative Number of Detected Molecules')
	
	#customize tick marks
	
	ax.tick_params(axis='x', which='both', direction='in',length=15,width=1)
	ax.tick_params(axis='y', which='both', direction='in',length=15,width=1)
	
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')
	
	#plot it
	
	ax.plot(years,dets,color='dodgerblue',linewidth=4)
	
	#add annotations
	
	#first, the detections per year
	
	args = {'ha' : 'left','size' : '24'}
	
	det_str = '\\noindent Since 1968: {:.1f} detections/year \\\\ Since 2005: {:.1f} detections/year' .format(trend1968,trend2005)
	
	ax.annotate(det_str, xy=(0.05,0.85),xycoords='axes fraction',**args)
	
	#Now the total number of detections
		
	ax.annotate('Total: {}' .format(dets[iyear(eyear)]),xy=(eyear-2,dets[iyear(eyear)]),xycoords='data',va='center',ha='right',size=16)
	
	#Now the facilities
	
	args = {'size' : '16'}
	
	arrowprops = {'arrowstyle' : '-|>', 'connectionstyle' : 'arc3', 'facecolor' : 'black'}
	
	ax.annotate('NRAO 36-foot (1968)',xy=(1968,dets[iyear(1968)]+4),xycoords='data',xytext=(0,35),textcoords='offset points',rotation=90,arrowprops=arrowprops,va='bottom',ha='center',**args)
	ax.annotate('',xy=(1982,dets[iyear(1982)]-4),xycoords='data',xytext=(0,-35),textcoords='offset points',rotation=0,arrowprops=arrowprops,va='bottom',ha='left',**args)
	ax.annotate('Nobeyama (1982)',xy=(1982,dets[iyear(1982)]-31),xycoords='data',**args)
	ax.annotate('IRAM (1984)',xy=(1984,dets[iyear(1984)]+4),xycoords='data',xytext=(0,35),textcoords='offset points',rotation=90,arrowprops=arrowprops,va='bottom',ha='center',**args)
	ax.annotate('GBT (2001)',xy=(2001,dets[iyear(2001)]+4),xycoords='data',xytext=(0,35),textcoords='offset points',rotation=90,arrowprops=arrowprops,va='bottom',ha='center',**args)
	ax.annotate('ALMA (2011)',xy=(2011,dets[iyear(2011)]-4),xycoords='data',xytext=(0,-35),textcoords='offset points',rotation=90,arrowprops=arrowprops,va='top',ha='center',**args)
		
	#show the plot
	
	plt.show()
	
	#write out the figure
	
	plt.savefig('cumulative_detections.pdf',format='pdf',transparent=True,bbox_inches='tight')
	
	return
	
def cumu_det_natoms_plot(list,syear=None,eyear=None):

	'''
	Makes a plot of the cumulative detections (sorted by atoms) by year using 'list', which is usually 'full_list'.  The start year and end year are defined by default to be the earliest year in the list and the current year + 20 (to give room for labels), but these are overridable. 
	
	'''

	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Cumulative Detections By Atoms')
	
	fig = plt.figure(num='Cumulative Detections By Atoms',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')	
	
	#get the starting and ending years, if they aren't set by the user
	
	if syear is None:
	
		#grab the earliest year in the list of molecules
	
		syear = min([x.year for x in list])
		
	if eyear is None:
	
		#grab the current year
	
		eyear = date.today().year
		
	#make the x-axis array of years
	
	years = np.arange(syear,eyear+1)
	
	#now we make a dictionary of the different traces we're gonna want, loop through the list, and populate an array of detections for that number of atoms or special case
	
	dets_dict = {}
	
	for natoms in range(2,14):
	
		#fill it with an empty array
	
		dets_dict[natoms] = np.copy(years)*0
		
		#loop through the years and the list and add everything up.  There's gotta be a better way to do this, but who cares, its fast enough.
	
		for x in range(len(dets_dict[natoms])):
	
			i = 0
		
			for mol in list:
		
				if mol.year < years[x]+1 and mol.natoms == natoms:
			
					i += 1
				
			dets_dict[natoms][x] = i
		
	#do the fullerenes and pahs
	
	dets_dict['fullerenes'] = np.copy(years)*0
	
	for x in range(len(dets_dict['fullerenes'])):

		i = 0
	
		for mol in list:
	
			if mol.year < years[x]+1 and mol.fullerene is True:
		
				i += 1
			
		dets_dict['fullerenes'][x] = i	
		
	dets_dict['pahs'] = np.copy(years)*0
	
	for x in range(len(dets_dict['pahs'])):

		i = 0
	
		for mol in list:
	
			if mol.year < years[x]+1 and mol.pah is True:
		
				i += 1
			
		dets_dict['pahs'][x] = i
		
	#load up an axis
	
	ax = fig.add_subplot(111)
	
	#label the axes
	
	plt.xlabel('Year')
	plt.ylabel('Cumulative Number of Detected Molecules')
	
	#customize tick marks
	
	ax.tick_params(axis='x', which='both', direction='in',length=15,width=1)
	ax.tick_params(axis='y', which='both', direction='in',length=15,width=1)
	
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')		
	ax.set_xticks([1940,1960, 1980,2000,2020])
	
	ax.set_xlim([syear,eyear+35])
	
	#plot the traces
	
	ax.plot(years,dets_dict[2],color='#000000')
	ax.plot(years,dets_dict[3],color='#800000')
	ax.plot(years,dets_dict[4],color='#f032e6')
	ax.plot(years,dets_dict[5],color='#9A6324')
	ax.plot(years,dets_dict[6],color='dodgerblue')
	ax.plot(years,dets_dict[7],color='#e6194B')
	ax.plot(years,dets_dict[8],color='#469990')
	ax.plot(years,dets_dict[9],color='#f58231')
	ax.plot(years,dets_dict[10],color='#42d4f4')
	ax.plot(years,dets_dict[11],color='#ffe119')
	ax.plot(years,dets_dict[12],color='#3cb44b')
	ax.plot(years,dets_dict[13],color='#e6beff')
	ax.plot(years,dets_dict['fullerenes'],color='#000075')
	ax.plot(years,dets_dict['pahs'],color='#aaffc3')
	
	#add the annotations	
	
	ax.annotate(r'\textbf{2 atoms}',xy=(eyear+2,dets_dict[2][-1]),xycoords='data',size=16,color='#000000',va='top')
	ax.annotate(r'\textbf{3 atoms}',xy=(eyear+2,dets_dict[3][-1]),xycoords='data',size=16,color='#800000',va='bottom')
	ax.annotate(r'\textbf{4 atoms}',xy=(eyear+2,dets_dict[4][-1]),xycoords='data',size=16,color='#f032e6',va='center')
	ax.annotate(r'\textbf{5 atoms}',xy=(eyear+2,dets_dict[5][-1]),xycoords='data',size=16,color='#9A6324',va='center')
	ax.annotate(r'\textbf{6 atoms}',xy=(eyear+2,dets_dict[6][-1]),xycoords='data',size=16,color='dodgerblue',va='center')
	ax.annotate(r'\textbf{7 atoms}',xy=(eyear+2,dets_dict[7][-1]),xycoords='data',size=16,color='#e6194B',va='center')
	ax.annotate(r'\textbf{8 atoms}',xy=(eyear+2,dets_dict[8][-1]),xycoords='data',size=16,color='#469990',va='bottom')
	ax.annotate(r'\textbf{9 atoms}',xy=(eyear+2,dets_dict[9][-1]),xycoords='data',size=16,color='#f58231',va='center')
	ax.annotate(r'\textbf{10 \&}',xy=(eyear+2,dets_dict[10][-1]),xycoords='data',size=16,color='#42d4f4',va='bottom')
	ax.annotate(r'\textbf{11 atoms}',xy=(eyear+11,dets_dict[11][-1]),xycoords='data',size=16,color='#ffe119',va='bottom')
	ax.annotate(r'\textbf{12 atoms}',xy=(eyear+2,dets_dict[12][-1]),xycoords='data',size=16,color='#3cb44b',va='center')
	ax.annotate(r'\textbf{13 atoms \& }',xy=(eyear+2,dets_dict[13][-1]),xytext=(eyear+2,dets_dict[13][-1]-0.5),xycoords='data',size=16,color='#e6beff',va='center')	
	arrowprops = {'arrowstyle' : '-|>', 'connectionstyle' : 'arc3', 'facecolor' : '#000075'}
	ax.annotate(r'\textbf{Fullerenes}',xy=(eyear+2,dets_dict['fullerenes'][-1]),xycoords='data',xytext=(eyear+17,dets_dict['fullerenes'][-1]),textcoords='data',arrowprops=arrowprops,size=16,color='#000075',va='center')
	ax.annotate(r'\textbf{PAHs}',xy=(eyear+22,dets_dict['pahs'][-1]-0.5),xycoords='data',size=16,color='#aaffc3',va='center')

	#show the plot
	
	plt.show()		
	
	#write out the figure
	
	plt.savefig('cumulative_by_atoms.pdf',format='pdf',transparent=True,bbox_inches='tight')
	
	return

def change_color(color, amount=1.0):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])		
	
def det_per_year_per_atom(list):

	'''
	Makes a plot of the average number of detections per year (y) for a molecule with (x) atoms, starting in the year they were first detected.  Has the ability to plot PAHs and fullerenes, but doesn't.
	'''

	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Detects Per Year Per Atom')
	
	fig = plt.figure(num='Detects Per Year Per Atom',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')	
	
	#We're going to cheat later on and lump PAHs and Fullerenes together, but these need fake natoms.  They'll always be +2 and +4, respectively, beyond the maximum of 'normal' molecules, and we'll leave +1 blank for visual separation.  For now, though, time to make a list and loop through our molecules.
	
	natoms = np.arange(2,18)
	
	#and an array of the average number of detections
	
	avg_dets = np.copy(natoms)*0.0
	n_dets = np.copy(natoms)*0.0
	
	#get the years that are being spanned
	
	eyear = max([x.year for x in list])
	
	for x in range(len(natoms)):
	
		i = 0
		
		years = []
		
		if natoms[x] < 14:
		
			for mol in list:
			
				if mol.natoms == natoms[x]:
			
					i += 1
					
					years.append(mol.year)
			
		elif natoms[x] == 15:
		
			for mol in list:
			
				if mol.pah is True:		
				
					i += 1
					
					years.append(mol.year)
					
		elif natoms[x] == 17:
		
			for mol in list:
			
				if mol.fullerene is True:
				
					i += 1
					
					years.append(mol.year)
							
		if i == 0:
			
			avg_dets[x] = np.NaN
			
		else:
			
			avg_dets[x] = i/(eyear-min(years)+1)
			n_dets[x] = i
			
	n_dets[n_dets == 0] = np.nan		
		
	#load up an axis
	
	ax = fig.add_subplot(111)
	
	#label the axes
	
	plt.xlabel('Number of Atoms')
	plt.ylabel('Detections/Year*')
	
	#customize tick marks
	
	ax.tick_params(axis='x', which='both', direction='in',length=15,width=1)
	ax.tick_params(axis='y', which='both', direction='in',length=15,width=1)
	
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')		
	ax.set_xticks([2,4,6,8,10,12,15,17])
	ax.set_xticklabels(['2','4','6','8','10','12','PAHs','Fullerenes'])
	
	ax.set_xlim([1,12.8])
	ax.set_ylim([0,1.])
	
	#plot the data
	
	marker = 'o'		
	mfc = 'dodgerblue'
	mec = change_color(mfc,1.5)	
	
	sizes = np.copy(n_dets) * 200/np.nanmin(n_dets)
	
# 	sizes = [x*200/np.amin(n_dets) for x in sizes if x > 0]
# 	
# 	ddd
	
	ax.scatter(natoms,avg_dets,marker='o',c=mfc,edgecolors=mec,s=sizes)
	
	#add some annotations
	
	ax.annotate('*Since year of first detection',xy=(0.35,0.9),xycoords='axes fraction',ha='left')
	ax.annotate('\\noindent Marker size proportional to\\\\ total \# of detections',xy=(0.05,0.1),xycoords='axes fraction',ha='left')
	
	for label in ax.get_xmajorticklabels():
	
		if label._text == 'Fullerenes' or label._text == 'PAHs':
		
			label.set_rotation(-45)
			
	fig.tight_layout()		

	#show the plot
	
	plt.show()
	
	#write out the figure
	
	plt.savefig('rate_by_atoms.pdf',format='pdf',transparent=True,bbox_inches='tight')	

	return	
	
def facility_shares(scopes_list,mols_list):

	'''
	Generates a plot of the percentage share of yearly detections that a facility contributed over its operational lifetime for the top 9 facilities
	'''

	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Facility Shares')
	
	fig,axs = plt.subplots(3,3,num='Facility Shares')
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':14, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')	

	#we need to generate the data now, which we'll store in a dictionary for each telescope.  Each entry will be [syear,eyear,ndetects,ntotal,shortname] for the start year, end year, number of detections, and total number of detections over those years
	
	detects = [x.ndetects for x in scopes_list]
	
	detects.sort(reverse=True)
	
	detects = detects[:9]
	
	min_allowed = min(detects)
	
	my_dict = {}
	
	for scope in scopes_list:
	
		ndetects = scope.ndetects #number of detections
		
		if ndetects >= min_allowed:
			#print(ndetects,scope.shortname)
			pass
			
		else:
		
			
			continue
		
		syear = scope.built #year it was built
		eyear = scope.decommissioned if scope.decommissioned is not None else date.today().year #year it was decommissioned or this year if it's still in operation

		#now we go get the total number of detections in that time
		
		ntotal = 0
		
		for mol in mols_list:
		
			if mol.year <= eyear and mol.year >= syear:
			
				ntotal += 1
		
		my_dict[scope.shortname] = [syear,eyear,ndetects,ntotal,scope.shortname]
		
	my_list = [my_dict[x] for x in my_dict]	
	
	my_list.sort(key = lambda x: x[2]/x[3],reverse=True)
	
# 	for x in my_list:
# 		print(x, x[2]/x[3])

	#load some axes
	
	idx = 0
	
	for i in range(3):
		
		for j in range(3)[::-1]:
		
			fracs = [my_list[idx][2]/my_list[idx][3],1.-my_list[idx][2]/my_list[idx][3]]
			
			label = '\\textbf' + '{' +'{}' .format(int(my_list[idx][2]/my_list[idx][3] * 100)) +'\%'
			
			if my_list[idx][1] == date.today().year:
			
				color = 'dodgerblue'
				
			else:
			
				color = '#F87070'
		


			if idx > 5:
			
				slices, labels = axs[i,j].pie(fracs, colors=[color,'#F5F6FF'], wedgeprops = {'linewidth' : 1.0, 'edgecolor' : 'black'})
			
				axs[i,j].annotate(my_list[idx][-1],xy=(0.5,1.08),xycoords='axes fraction',ha='center',size=12)
				axs[i,j].annotate('{} - {}' .format(my_list[idx][0],my_list[idx][1]),xy=(0.5,.95),xycoords='axes fraction',ha='center',size=10)
			
				kw = dict(arrowprops=dict(arrowstyle="-"), zorder=0, va="center")
				ang = (slices[0].theta2 - slices[0].theta1)/2. + slices[0].theta1
				y = np.sin(np.deg2rad(ang))
				x = np.cos(np.deg2rad(ang))
				horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
				connectionstyle = "angle,angleA=0,angleB={}".format(ang)
				kw["arrowprops"].update({"connectionstyle": connectionstyle})
				axs[i,j].annotate(label, xy=(x, y), xytext=(1.1*np.sign(x), y),
					horizontalalignment=horizontalalignment, size=10, **kw)
					
			else:
			
				my_labels = ['\\textbf' + '{' +'{}' .format(int(my_list[idx][2]/my_list[idx][3] * 100)) +'\%','']
			
				slices, labels = axs[i,j].pie(fracs, labels=my_labels, colors=[color,'#F5F6FF'], wedgeprops = {'linewidth' : 1.0, 'edgecolor' : 'black'}, labeldistance=0.7)
			
				axs[i,j].annotate(my_list[idx][-1],xy=(0.5,1.08),xycoords='axes fraction',ha='center',size=12)
				axs[i,j].annotate('{} - {}' .format(my_list[idx][0],my_list[idx][1]),xy=(0.5,.95),xycoords='axes fraction',ha='center',size=10)
				
				for label in labels:
				
					label.set_horizontalalignment('center')
					label.set_color('white')
					label.set_fontsize(10)
				
			
			#draw a line around the outside of the pie
						
			# get the center and radius of the pie wedges			
			center = slices[0].center
			r = slices[0].r
			
			# create a new circle with the desired properties
			circle = matplotlib.patches.Circle(center, r, fill=False, edgecolor="black", linewidth=1)
			# add the circle to the axes
			axs[i,j].add_patch(circle)
			
			idx += 1

	fig.tight_layout()
	
	fig.subplots_adjust(wspace=-.5, hspace=.15)
	
	plt.show()
	
	plt.savefig('facility_shares.pdf',format='pdf',transparent=True,bbox_inches='tight')	

	return		
	
def cumu_det_facility(list):

	'''
	Makes a plot of the cumulative number of detections of a facility with time.
	'''

	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Detections Per Facility Over Time')
	
	fig = plt.figure(num='Detections Per Facility Over Time',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')	
	
	#We're only going to do facilities with 10 or more total detections.  Right now that's the GBT, IRAM, NRAO 140-ft, NRAO/ARO 12-m, NRAO 36ft, and Nobeyama.
	
	#set this thing up to kick off a few years before 1968 and run until today
	
	years = np.arange(1965,date.today().year+1)
	
	scopes = [GBT,IRAM30,NRAO140,NRAOARO12,NRAO36,Nobeyama45]
	
	my_dict = {}
	
	for scope in scopes:
	
		tmp_years = np.copy(years)*0
		
		i = 0
		
		for x in range(len(years)):
	
			for mol in list:
	
				if mol.year == years[x] and scope in mol.telescopes:
		
					i += 1
			
			tmp_years[x] = i
		
		my_dict[scope.shortname] = tmp_years
		
	
	ax = fig.add_subplot(111)
	
	#label the axes
	
	plt.xlabel('Year')
	plt.ylabel('Cumulative Number of Detected Molecules')
	
	#customize tick marks
	
	ax.tick_params(axis='x', which='both', direction='in',length=15,width=1)
	ax.tick_params(axis='y', which='both', direction='in',length=15,width=1)
	
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')		
	ax.set_xticks([1970,1980,1990,2000,2010,2020])
	
	ax.plot(years,my_dict['GBT'],color='#000000')
	ax.plot(years,my_dict['IRAM'],color='#800000')
	ax.plot(years,my_dict['NRAO 140-ft'],color='#f032e6')
	ax.plot(years,my_dict['NRAO/ARO 12-m'],color='dodgerblue')
	ax.plot(years,my_dict['NRAO 36-ft'],color='#e6194B')
	ax.plot(years,my_dict['Nobeyama'],color='#469990')
	
	ax.annotate(r'\textbf{GBT}',xy=(date.today().year+1,my_dict['GBT'][-1]),xycoords='data',size=16,color='#000000',va='center',zorder=100)
	ax.annotate(r'\textbf{IRAM}',xy=(date.today().year+1,my_dict['IRAM'][-1]),xycoords='data',size=16,color='#800000',va='center',zorder=100)
	ax.annotate(r'\textbf{NRAO 140-ft}',xy=(date.today().year+1,my_dict['NRAO 140-ft'][-1]),xycoords='data',size=16,color='#f032e6',va='top',zorder=100)
	ax.annotate(r'\textbf{NRAO/ARO 12-m}',xy=(date.today().year+1,my_dict['NRAO/ARO 12-m'][-1]),xycoords='data',size=16,color='dodgerblue',va='center',zorder=100)
	ax.annotate(r'\textbf{NRAO 36-ft}',xy=(date.today().year+1,my_dict['NRAO 36-ft'][-1]),xycoords='data',size=16,color='#e6194B',va='center',zorder=100)
	ax.annotate(r'\textbf{Nobeyama}',xy=(date.today().year+1,my_dict['Nobeyama'][-1]),xycoords='data',size=16,color='#469990',va='bottom',zorder=100)
	
	ax.set_xlim([1965,date.today().year+19])
	
	ax.set_ylim(0,max([my_dict[x][-1] for x in my_dict])+5)
	
	#do linear fits to the data for the ranges we care about for each facility:
	
	#get some year indicies for years we care about
	
	def iyear(x):
	
		return np.argwhere(years==x)[0][0]	
	
	trendGBT = np.polynomial.polynomial.Polynomial.fit(years[iyear(GBT.built):],my_dict['GBT'][iyear(GBT.built):],1).convert().coef[1]
	ax.annotate('{:.1f}/yr' .format(trendGBT),xy=(2014,7.5),xycoords='data',size=16,color='#000000',ha='center')
	ax.annotate('{}{: <7}' .format(GBT.built,' - '),xy=(2014,4.5),xycoords='data',size=16,color='#000000',ha='right')
	
	trendIRAMold = np.polynomial.polynomial.Polynomial.fit(years[iyear(IRAM30.built):iyear(2006)],my_dict['IRAM'][iyear(IRAM30.built):iyear(2006)],1).convert().coef[1]
	ax.annotate('{:.1f}/yr' .format(trendIRAMold),xy=(1990,21),xycoords='data',size=16,color='#800000',ha='center')
	ax.annotate('{} - 2006' .format(IRAM30.built),xy=(1990,18),xycoords='data',size=16,color='#800000',ha='center')

	trendIRAMnew = np.polynomial.polynomial.Polynomial.fit(years[iyear(2006):],my_dict['IRAM'][iyear(2006):],1).convert().coef[1]
	ax.annotate('{:.1f}/yr' .format(trendIRAMnew),xy=(2020,45),xycoords='data',size=16,color='#800000',ha='center')
	ax.annotate('2006{: <7}' .format(' - '),xy=(2020,42),xycoords='data',size=16,color='#800000',ha='right')	
	
	trend140 = np.polynomial.polynomial.Polynomial.fit(years[iyear(NRAO140.built):1993],my_dict['NRAO 140-ft'][iyear(NRAO140.built):1993],1).convert().coef[1]
	ax.annotate('{:.1f}/yr' .format(trend140),xy=(1980,12.5),xycoords='data',size=16,color='#f032e6',ha='center')
	ax.annotate('{} - 1993' .format(NRAO140.built),xy=(1980,9.5),xycoords='data',size=16,color='#f032e6',ha='center')	
	
	trend12 = np.polynomial.polynomial.Polynomial.fit(years[iyear(NRAOARO12.built):],my_dict['NRAO/ARO 12-m'][iyear(NRAOARO12.built):],1).convert().coef[1]
	ax.annotate('{:.1f}/yr' .format(trend12),xy=(2014.8,30.2),xycoords='data',size=16,color='dodgerblue',ha='center')
	ax.annotate('{}{: <7}' .format(NRAOARO12.built,' - '),xy=(2014.8,27.2),xycoords='data',size=16,color='dodgerblue',ha='right')	
	
	trend36 = np.polynomial.polynomial.Polynomial.fit(years[iyear(NRAO36.built):1985],my_dict['NRAO 36-ft'][iyear(NRAO36.built):1985],1).convert().coef[1]
	ax.annotate('{:.1f}/yr' .format(trend36),xy=(1975,34),xycoords='data',size=16,color='#e6194B',ha='center')
	ax.annotate('{} - 1985' .format(NRAO36.built),xy=(1975,31),xycoords='data',size=16,color='#e6194B',ha='center')		
	
	trendNobeyama = np.polynomial.polynomial.Polynomial.fit(years[iyear(Nobeyama45.built):],my_dict['Nobeyama'][iyear(Nobeyama45.built):],1).convert().coef[1]
	ax.annotate('{:.1f}/yr' .format(trendNobeyama),xy=(2012,19),xycoords='data',size=16,color='#469990',ha='center')
	ax.annotate('{}{: <7}' .format(Nobeyama45.built,' - '),xy=(2012,16),xycoords='data',size=16,color='#469990',ha='right')		
	
	plt.show()
	
	plt.savefig('scopes_by_year.pdf',format='pdf',transparent=True,bbox_inches='tight')					
	
	return
	
def periodic_heatmap(mol_list):


	'''
	Makes a periodic table heat map
	'''

	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Periodic Heatmap')
	
	fig = plt.figure(num='Periodic Heatmap',figsize=(20,9.5))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')	
	
	#make a dictionary of detections from the list
	
	census = {
		'H': 0,
		'He': 0,
		'C': 0,
		'O': 0,
		'N': 0,
		'S': 0,
		'P': 0,
		'Si': 0,
		'Cl': 0,
		'F': 0,
		'Mg': 0,
		'Na': 0,
		'Al': 0,
		'K': 0,
		'Fe': 0,
		'Ti': 0,
		'Ar': 0,
		'V': 0,
		'Ca': 0,
		}
		
	els = ['H', 'He', 'C', 'O', 'N', 'S', 'P', 'Si', 'Cl', 'F', 'Mg', 'Na', 'Al', 'K', 'Fe', 'Ti', 'Ar', 'V', 'Ca']
	
	for mol in mol_list:
	
		for el in els:
			
			if getattr(mol,el) > 0:
		
				census[el] += 1
				
	maxdets = max([census[x] for x in census])
	
	map_colors = list(Color("#f1fb53").range_to(Color("#f00707"),maxdets))

	map_colors = [str(x) for x in map_colors]		
	
	#Dictionary for the periodic table
	
	elements = {
		'H'		:	[pt.H,1,6.9],
		'He'	:	[pt.He,18,6.9],
		'Li'	:	[pt.Li,1,5.75],
		'Be'	:	[pt.Be,2,5.75],
		'B'		:	[pt.B,13,5.75],
		'C'		:	[pt.C,14,5.75],
		'N'		:	[pt.N,15,5.75],
		'O'		:	[pt.O,16,5.75],
		'F'		:	[pt.F,17,5.75],
		'Ne'	:	[pt.Ne,18,5.75],
		'Na'	:	[pt.Na,1,4.6],
		'Mg'	:	[pt.Mg,2,4.6],
		'Al'	:	[pt.Al,13,4.6],
		'Si'	:	[pt.Si,14,4.6],
		'P'		:	[pt.P,15,4.6],
		'S'		:	[pt.S,16,4.6],
		'Cl'	:	[pt.Cl,17,4.6],
		'Ar'	:	[pt.Ar,18,4.6],
		'K'		:	[pt.K,1,3.45],
		'Ca'	:	[pt.Ca,2,3.45],
		'Sc'	:	[pt.Sc,3,3.45],
		'Ti'	:	[pt.Ti,4,3.45],
		'V'		:	[pt.V,5,3.45],
		'Cr'	:	[pt.Cr,6,3.45],
		'Mn'	:	[pt.Mn,7,3.45],
		'Fe'	:	[pt.Fe,8,3.45],
		'Co'	:	[pt.Co,9,3.45],
		'Ni'	:	[pt.Ni,10,3.45],
		'Cu'	:	[pt.Cu,11,3.45],
		'Zn'	:	[pt.Zn,12,3.45],
		'Ga'	:	[pt.Ga,13,3.45],
		'Ge'	:	[pt.Ge,14,3.45],
		'As'	:	[pt.As,15,3.45],
		'Se'	:	[pt.Se,16,3.45],
		'Br'	:	[pt.Br,17,3.45],
		'Kr'	:	[pt.Kr,18,3.45],
		'Rb'	:	[pt.Rb,1,2.3],
		'Sr'	:	[pt.Sr,2,2.3],
		'Y'		:	[pt.Y,3,2.3],
		'Zr'	:	[pt.Zr,4,2.3],
		'Nb'	:	[pt.Nb,5,2.3],
		'Mo'	:	[pt.Mo,6,2.3],
		'Tc'	:	[pt.Tc,7,2.3],
		'Ru'	:	[pt.Ru,8,2.3],
		'Rh'	:	[pt.Rh,9,2.3],
		'Pd'	:	[pt.Pd,10,2.3],
		'Ag'	:	[pt.Ag,11,2.3],
		'Cd'	:	[pt.Cd,12,2.3],
		'In'	:	[pt.In,13,2.3],
		'Sn'	:	[pt.Sn,14,2.3],
		'Sb'	:	[pt.Sb,15,2.3],
		'Te'	:	[pt.Te,16,2.3],
		'I'		:	[pt.I,17,2.3],
		'Xe'	:	[pt.Xe,18,2.3],
		'Cs'	:	[pt.Cs,1,1.15],
		'Ba'	:	[pt.Ba,2,1.15],
		'Hf'	:	[pt.Hf,4,1.15],
		'Ta'	:	[pt.Ta,5,1.15],
		'W'		:	[pt.W,6,1.15],
		'Re'	:	[pt.Re,7,1.15],
		'Os'	:	[pt.Os,8,1.15],
		'Ir'	:	[pt.Ir,9,1.15],
		'Pt'	:	[pt.Pt,10,1.15],
		'Au'	:	[pt.Au,11,1.15],
		'Hg'	:	[pt.Hg,12,1.15],
		'Tl'	:	[pt.Tl,13,1.15],
		'Pb'	:	[pt.Pb,14,1.15],
		'Bi'	:	[pt.Bi,15,1.15],
		'Po'	:	[pt.Po,16,1.15],
		'At'	:	[pt.At,17,1.15],
		'Rn'	:	[pt.Rn,18,1.15],
		'Fr'	:	[pt.Fr,1,0.],
		'Ra'	:	[pt.Ra,2,0.],
		'Rf'	:	[pt.Rf,4,0.],
		'Db'	:	[pt.Db,5,0.],
		'Sg'	:	[pt.Sg,6,0.],
		'Bh'	:	[pt.Bh,7,0.],
		'Hs'	:	[pt.Hs,8,0.],
		'Mt'	:	[pt.Mt,9,0.],
		'Ds'	:	[pt.Ds,10,0.],
		'Rg'	:	[pt.Rg,11,0.],
		'Cn'	:	[pt.Cn,12,0.],
		'Nh'	:	[pt.Nh,13,0.],
		'Fl'	:	[pt.Fl,14,0.],
		'Mc'	:	[pt.Mc,15,0.],
		'Lv'	:	[pt.Lv,16,0.],
		'Ts'	:	[pt.Ts,17,0.],
		'Og'	:	[pt.Og,18,0.],
		}
	
	#load up an axis
	
	ax = plt.axes([0,0,1,1])
	
	ax.set_xlim([0,18])
	ax.set_ylim([0,8])
	
	for el in elements:
		
		x = elements[el][1]-1
		y = elements[el][2]
		
		sym = '\\textbf{' + elements[el][0].symbol + '}'
		num = elements[el][0].number
		mass = elements[el][0].mass
		name = elements[el][0].name.capitalize()
		
		this_color = 'white'
			
		if el in census:
		
			this_color = map_colors[census[el]-1]
		
			ndets = '\\textbf{' + str(census[el]) + '}'
		
			ax.annotate(ndets,xy=(x+0.8,y+.95),xycoords='data',size=14,color='black',ha='right',va='top')

		rect = patches.Rectangle((x,y),0.9,1.05,linewidth=1,edgecolor='black',facecolor=this_color,alpha=0.5)
		ax.add_patch(rect)
		ax.annotate(num,xy=(x+0.1,y+.95),xycoords='data',size=14,color='black',ha='left',va='top')
		ax.annotate(sym,xy=(x+0.1,y+.70),xycoords='data',size=20,color='black',ha='left',va='top')
		ax.annotate(mass,xy=(x+0.1,y+.42),xycoords='data',size=8,color='black',ha='left',va='top')
		ax.annotate(name,xy=(x+0.1,y+.29),xycoords='data',size=8,color='black',ha='left',va='top')
			
	
	plt.axis('equal')
	plt.axis('off')	
	plt.show()
	
	#write out the figure
	
	plt.savefig('periodic_heatmap.pdf',format='pdf',transparent=True,bbox_inches='tight',pad_inches=0)
	
	#the bit below crops off extra white space.  This only works on Macs with the TexLive pdfcrop utility installed.  Comment out if not desired.
	
	os.system('pdfcrop --margins -0 periodic_heatmap.pdf periodic_heatmap.pdf')
	
	return	
	
def mass_by_wavelength(list):

	'''
	Makes a KDE plot of detections at each wavelength vs mass
	'''
	
	#gather the data
	
	my_dict = {
	
		'UV'	: [],
		'Vis'	: [],
		'IR'	: [],
		'sub-mm': [],
		'mm'	: [],
		'cm'	: [],
		'UV-Vis': [],
		
	}
	
	for x in list:
	
		for y in my_dict:
		
			if y in x.wavelengths:
			
				my_dict[y].append(x.mass)
				
	for x in my_dict['UV']:
	
		my_dict['UV-Vis'].append(x)
		
	for x in my_dict['Vis']:
	
		my_dict['UV-Vis'].append(x)
	
	plt.close('Detections at Wavelengths by Mass')
	
	fig = plt.figure(num='Detections at Wavelengths by Mass',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')		

	#load up an axis
	
	ax = fig.add_subplot(111)

	ax.tick_params(axis='x', which='both', direction='in',length=5,width=1)
	ax.tick_params(axis='y', which='both', direction='in',length=5,width=1)	

	#label the axes
	
	plt.xlabel('Atomic Mass (amu)')
	plt.ylabel('Kernel Density Estimate')
	
	#Do the estimates and plot them
	
	xvals = np.arange(0,160)
	
	density_cm = gkde(my_dict['cm'])
	density_cm.covariance_factor = lambda : .5
	density_cm._compute_covariance()
	
	density_mm = gkde(my_dict['mm'])
	density_mm.covariance_factor = lambda : .5
	density_mm._compute_covariance()
	
	density_submm = gkde(my_dict['sub-mm'])
	density_submm.covariance_factor = lambda : .5
	density_submm._compute_covariance()	
	
	density_IR = gkde(my_dict['IR'])
	density_IR.covariance_factor = lambda : .5
	density_IR._compute_covariance()	
	
	density_UV = gkde(my_dict['UV-Vis'])
	density_UV.covariance_factor = lambda : .5
	density_UV._compute_covariance()	
	
	ax.plot(xvals,density_cm(xvals),color='dodgerblue')
	ax.fill_between(xvals,density_cm(xvals),0,facecolor='dodgerblue',alpha=0.25,zorder=4)
	ax.annotate('{}' .format(len(my_dict['cm'])), xy=(75,0.01),xycoords='data',ha='left',va='bottom',color='dodgerblue')	
	
	max_mm = max([x for x in my_dict['mm'] if x < 160])
	ax.plot(xvals[:max_mm+1],density_mm(xvals[:max_mm+1]),color='darkorange')
	ax.fill_between(xvals[:max_mm+1],density_mm(xvals[:max_mm+1]),0,facecolor='darkorange',alpha=0.25)
	ax.annotate('{}' .format(len(my_dict['mm'])), xy=(53.5,0.022),xycoords='data',ha='left',va='bottom',color='darkorange')
	
	ax.plot(xvals,density_submm(xvals),color='forestgreen')
	ax.fill_between(xvals,density_submm(xvals),0,facecolor='forestgreen',alpha=0.25)
	ax.annotate('{}' .format(len(my_dict['sub-mm'])), xy=(40,0.038),xycoords='data',ha='left',va='bottom',color='forestgreen')
	
	max_IR = max([x for x in my_dict['IR'] if x < 160])
	ax.plot(xvals[:max_IR+1],density_IR(xvals[:max_IR+1]),color='black')
	ax.fill_between(xvals[:max_IR+1],density_IR(xvals[:max_IR+1]),0,facecolor='black',alpha=0.25,zorder=5)
	ax.annotate('{}' .format(len(my_dict['IR'])), xy=(68,0.0025),xycoords='data',ha='left',va='bottom',color='black')		
	
	max_UV = max(my_dict['UV-Vis'])
	ax.plot(xvals[:max_UV+1],density_UV(xvals[:max_UV+1]),color='violet')
	ax.fill_between(xvals[:max_UV+1],density_UV(xvals[:max_UV+1]),0,facecolor='violet',alpha=0.25)
	ax.annotate('{}' .format(len(my_dict['UV-Vis'])), xy=(16.5,0.057),xycoords='data',ha='left',va='bottom',color='violet')
	
	ax.annotate(r'\underline{Detection Wavelengths}',xy=(158,0.06),xycoords='data',color='black',ha='right',va='top')
	ax.annotate('centimeter',xy=(158,0.056),xycoords='data',color='dodgerblue',ha='right',va='top')
	ax.annotate('millimeter',xy=(158,0.052),xycoords='data',color='darkorange',ha='right',va='top')
	ax.annotate('sub-millimeter',xy=(158,0.048),xycoords='data',color='forestgreen',ha='right',va='top')
	ax.annotate('infrared',xy=(158,0.044),xycoords='data',color='black',ha='right',va='top')
	ax.annotate('visible/ultraviolet',xy=(158,0.04),xycoords='data',color='violet',ha='right',va='top')		
	
	plt.show()
	
	plt.savefig('mass_by_wavelengths_kde.pdf',format='pdf',transparent=True,bbox_inches='tight',pad_inches=0)
	
	return			

def mols_waves_by_atoms(list):

	'''
	Makes six histogram plots of molecules detected in each wavelength range by number of atoms, excepting fullerenes
	'''

	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Molecules Detected in Each Wavelength by Number of Atoms')
	
	fig = plt.figure(num='Molecules Detected in Each Wavelength by Number of Atoms',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':18, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')	
	
	#gather the data
	
	my_dict = {
	
		'UV'	: [],
		'Vis'	: [],
		'IR'	: [],
		'sub-mm': [],
		'mm'	: [],
		'cm'	: [],	
		
	}
	
	max_n = []
	
	for x in list:
	
		for y in my_dict:
		
			if y in x.wavelengths:
			
				if x.fullerene is True:
				
					continue
					
				else:
			
					my_dict[y].append(x.natoms)
					
					max_n.append(x.natoms)
		
					
	max_n = max(max_n)				
				
	ax1 = plt.subplot(231)
	ax2 = plt.subplot(232)
	ax3 = plt.subplot(233)
	ax4 = plt.subplot(234)
	ax5 = plt.subplot(235)
	ax6 = plt.subplot(236)
	
	n_bins = max_n
	
	xvals = np.arange(0,max_n+1,0.5)
	
	ax1.tick_params(axis='x', which='both', direction='in',length=5,width=1)
	ax1.tick_params(axis='y', which='both', direction='in',length=5,width=1)
	ax2.tick_params(axis='x', which='both', direction='in',length=5,width=1)
	ax2.tick_params(axis='y', which='both', direction='in',length=5,width=1)
	ax3.tick_params(axis='x', which='both', direction='in',length=5,width=1)
	ax3.tick_params(axis='y', which='both', direction='in',length=5,width=1)
	ax4.tick_params(axis='x', which='both', direction='in',length=5,width=1)
	ax4.tick_params(axis='y', which='both', direction='in',length=5,width=1)
	ax5.tick_params(axis='x', which='both', direction='in',length=5,width=1)
	ax5.tick_params(axis='y', which='both', direction='in',length=5,width=1)
	ax6.tick_params(axis='x', which='both', direction='in',length=5,width=1)
	ax6.tick_params(axis='y', which='both', direction='in',length=5,width=1)					
		
	density_cm = gkde(my_dict['cm'])
	density_cm.covariance_factor = lambda : .5
	density_cm._compute_covariance()
	ax1.plot(xvals,density_cm(xvals))	
	ax1.fill_between(xvals,density_cm(xvals),0,facecolor='dodgerblue',alpha=0.25)
	ax1.annotate('cm',xy=(0.95,0.96),xycoords='axes fraction',ha='right',va='top',size=24)
	ax1.set_ylim([0,0.65])
	ax1.set_xticks([0,5,10,15,20])	
	ax1.set_xticklabels([])
	ax1.set_ylabel('Kernel Density Estimate')

	density_mm = gkde(my_dict['mm'])
	density_mm.covariance_factor = lambda : .5
	density_mm._compute_covariance()
	ax2.plot(xvals,density_mm(xvals))
	ax2.fill_between(xvals,density_mm(xvals),0,facecolor='dodgerblue',alpha=0.25)
	ax2.annotate('mm',xy=(0.95,0.96),xycoords='axes fraction',ha='right',va='top',size=24)
	ax2.set_ylim([0,0.65])
	ax2.set_xticks([0,5,10,15,20])
	ax2.set_xticklabels([])

	density_submm = gkde(my_dict['sub-mm'])
	density_submm.covariance_factor = lambda : .5
	density_submm._compute_covariance()
	ax3.plot(xvals,density_submm(xvals))
	ax3.fill_between(xvals,density_submm(xvals),0,facecolor='dodgerblue',alpha=0.25)
	ax3.annotate('sub-mm',xy=(0.95,0.96),xycoords='axes fraction',ha='right',va='top',size=24)
	ax3.set_ylim([0,0.65])
	ax3.set_xticks([0,5,10,15,20])
	ax3.set_xticklabels([])

	density_IR = gkde(my_dict['IR'])
	density_IR.covariance_factor = lambda : .5
	density_IR._compute_covariance()
	ax4.plot(xvals,density_IR(xvals))
	ax4.fill_between(xvals,density_IR(xvals),0,facecolor='dodgerblue',alpha=0.25)
	ax4.annotate('IR',xy=(0.95,0.96),xycoords='axes fraction',ha='right',va='top',size=24)
	ax4.set_ylim([0,0.65])
	ax4.set_xticks([0,5,10,15,20])
	ax4.set_xlabel('\# of Atoms')
	ax4.set_ylabel('Kernel Density Estimate')

	
	ax5.hist(my_dict['Vis'],bins=[.5,1.5,2.5,3.5,4.5])
	ax5.annotate('Vis',xy=(0.95,0.96),xycoords='axes fraction',ha='right',va='top',size=24)
	ax5.set_xlim([0,20])
	ax5.set_ylim([0,8])	
	ax5.set_xticks([0,5,10,15,20])
	ax5.set_xticklabels([])
	ax5.set_ylabel('\# of Detected Molecules')
	
	ax6.hist(my_dict['UV'],bins=[.5,1.5,2.5,3.5,4.5])
	ax6.annotate('UV',xy=(0.95,0.96),xycoords='axes fraction',ha='right',va='top',size=24)
	ax6.set_xlim([0,20])
	ax6.set_ylim([0,8])	
	ax6.set_xticks([0,5,10,15,20])
	ax6.set_xticklabels([])
	
	plt.tight_layout()
	
	plt.show()
	
	plt.savefig('mols_waves_by_atoms.pdf',format='pdf',transparent=True,bbox_inches='tight',pad_inches=0)

	return
	
def du_histogram(list):

	'''
	Makes a histogram of the degree of unsaturation of molecules containing only H, O, N, C, Cl, or F.
	'''

	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Degree of Unsaturation Histogram')
	
	fig = plt.figure(num='Degree of Unsaturation Histogram',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')
	
	#gather the data
	
	dus = []
	
	for x in list:
	
		if (x.H + x.O + x.N + x.C + x.Cl + x.F + x.S) == x.natoms and x.du is not None and x.fullerene is not True:
		
			dus.append(x.du)
			
	#set up a plot
	
	ax = plt.subplot(111)
	
	ax.tick_params(axis='x', which='both', direction='in',length=5,width=1)
	ax.tick_params(axis='y', which='both', direction='in',length=5,width=1)
	
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')
	
	plt.xlabel('Degree of Unsaturation')
	plt.ylabel('\# of Detected Molecules')
	
	bins = np.arange(-0.25,12.5,0.5)	
	(n,bins,patches) = ax.hist(dus,bins=bins,facecolor='dodgerblue',alpha=0.25)
	ax.hist(dus,bins=bins,edgecolor='royalblue',linewidth=1.5,fill=False)
	
	ax.annotate(r'\ce{CH4}, \ce{CH3OH}, \ce{CH3Cl}, \ce{CH3NH2}, ...',xy=(0,n[0]+1),xycoords='data',rotation=90,size=16,ha='center',va='bottom')		
	ax.annotate(r'\ce{HC11N}',xy=(12,n[12*2]+1),xycoords='data',rotation=90,size=16,ha='center',va='bottom')		
	
	plt.tight_layout()
	plt.show()
	
	plt.savefig('du_histogram.pdf',format='pdf',transparent=True,bbox_inches='tight')

	return	
	
def type_pie_chart(my_list):

	'''
	Makes a pie chart of the fraction of interstellar molecules that are neutral, radical, cation, cyclic, pahs, fullerenes, or anions
	'''

	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Type Pie Chart')
	
	fig = plt.figure(num='Type Pie Chart',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')
	
	#gather the data

	my_dict = {
	
		'Neutral' 	: 	[0],
		'Radical'	:	[0],
		'Cation'	:	[0],
		'Cyclic'	:	[0],
		'Anion'		:	[0],
		'Fullerene'	:	[0],
		'PAH'		:	[0],
		
		}
		
	for mol in my_list:
	
		if mol.neutral is True:
		
			my_dict['Neutral'][0] += 1 
			
		if mol.radical is True:
		
			my_dict['Radical'][0] += 1 
			
		if mol.cation is True:
		
			my_dict['Cation'][0] += 1 
			
		if mol.cyclic is True:
		
			my_dict['Cyclic'][0] += 1
			
		if mol.anion is True:
		
			my_dict['Anion'][0] += 1
			
		if mol.fullerene is True:
		
			my_dict['Fullerene'][0] += 1
			
		if mol.pah is True:
		
			my_dict['PAH'][0] += 1		
	
	nmols = len(my_list)
	
	for type in my_dict:
	
		my_dict[type].append(my_dict[type][0]/nmols)		
		
	labels = ['Neutral','Radical','Cation','Cyclic','Anion','Fullerene','PAH']
	
	fracs = [my_dict[x][1] for x in labels]			

	#set up a plot
	
	ax = plt.subplot(111)
	
	size = 0.1
	
	def getshift(x):
	
		return -(90-(360 - 360*x)/2)
		
	def getper(x):
	
		return '{}' .format(int(100*x)) + r'\%'
	
	ax.pie([fracs[0],1.-fracs[0]],  colors=['dodgerblue','#EEEEEE'], radius=1, startangle=getshift(fracs[0]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[1],1.-fracs[1]],  colors=['darkorange','#EEEEEE'], radius=1-size-.02, startangle=getshift(fracs[1]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[2],1.-fracs[2]],  colors=['forestgreen','#EEEEEE'], radius=1-2*size-.04, startangle=getshift(fracs[2]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[3],1.-fracs[3]],  colors=['violet','#EEEEEE'], radius=1-3*size-.06, startangle=getshift(fracs[3]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[4],1.-fracs[4]],  colors=['red','#EEEEEE'], radius=1-4*size-.08, startangle=getshift(fracs[4]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[5],1.-fracs[5]],  colors=['goldenrod','#EEEEEE'], radius=1-5*size-.1, startangle=getshift(fracs[5]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[6],1.-fracs[6]],  colors=['royalblue','#EEEEEE'], radius=1-6*size-.12, startangle=getshift(fracs[6]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))

	ax.annotate(r'\textbf{Neutrals}',xy=(0.5,0.11),xycoords='axes fraction',color='dodgerblue', ha='center',size=14)
	ax.annotate(r'\textbf{Radicals}',xy=(0.5,0.16),xycoords='axes fraction',color='darkorange', ha='center',size=14)
	ax.annotate(r'\textbf{Cations}',xy=(0.5,0.205),xycoords='axes fraction',color='forestgreen', ha='center',size=14)
	ax.annotate(r'\textbf{Cyclics}',xy=(0.5,0.255),xycoords='axes fraction',color='violet', ha='center',size=14)
	ax.annotate(r'\textbf{Anions}',xy=(0.5,0.305),xycoords='axes fraction',color='red', ha='center',size=14)
	ax.annotate(r'\textbf{Fullerenes}',xy=(0.5,0.3575),xycoords='axes fraction',color='goldenrod', ha='center',size=14)
	ax.annotate(r'\textbf{PAHs}',xy=(0.5,0.40),xycoords='axes fraction',color='royalblue', ha='center',size=14)
	
	percents = ['\\textbf{' + '{:.1f}' .format((x*100)) + '}\%'  for x in fracs]
	
	ax.annotate(percents[0],xy=(0.775,0.775),xycoords='axes fraction',color='white', ha='center',size=12,rotation=-45)
	ax.annotate(percents[1],xy=(0.74,0.74),xycoords='axes fraction',color='darkorange', ha='center',size=12,rotation=-45)
	ax.annotate(percents[2],xy=(0.705,0.705),xycoords='axes fraction',color='forestgreen', ha='center',size=12,rotation=-45)
	ax.annotate(percents[3],xy=(0.67,0.67),xycoords='axes fraction',color='violet', ha='center',size=12,rotation=-45)
	ax.annotate(percents[4],xy=(0.635,0.635),xycoords='axes fraction',color='red', ha='center',size=12,rotation=-45)
	ax.annotate(percents[5],xy=(0.6,0.6),xycoords='axes fraction',color='goldenrod', ha='center',size=12,rotation=-45)
	ax.annotate(percents[6],xy=(0.5675,0.5675),xycoords='axes fraction',color='royalblue', ha='center',size=12,rotation=-45)
	
	
	plt.tight_layout()
	plt.show()
	
	plt.savefig('type_pie_chart.pdf',format='pdf',transparent=True,bbox_inches='tight',pad_inches=-.65)

	return	

def source_pie_chart(my_list):

	'''
	Makes a pie chart of the fraction of interstellar molecules detected in carbon stars, dark clouds, los clouds, star forming regions, and other types of sources
	'''

	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Source Pie Chart')
	
	fig = plt.figure(num='Source Pie Chart',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')
	
	#gather the data

	my_dict = {
	
		'Carbon Star' 	: 	[0],
		'Dark Cloud'	:	[0],
		'LOS Cloud'	:	[0],
		'SFR'	:	[0],
		'Other'		:	[0],		
		}
		
	#we have to be a little careful here, because for a given source, there can be two SFRs listed, and we only want to credit it once
		
	for mol in my_list:
	
		#we'll make a dictionary here to flag if we've credited things already
		
		credit_dict = {
	
			'Carbon Star' 	: 	False,
			'Dark Cloud'	:	False,
			'LOS Cloud'	:	False,
			'SFR'	:	False,
			'Other'		:	False,		
			}		
	
		for source in mol.sources:
		
			if source.type in my_dict:
			
				if credit_dict[source.type] is False:
			
					my_dict[source.type][0] += 1
					
					credit_dict[source.type] = True
				
			else:
			
				if credit_dict['Other'] is False:
			
					my_dict['Other'][0] += 1	
					
					credit_dict['Other'] = True
	
	nmols = len(my_list)
	
	for type in my_dict:
	
		my_dict[type].append(my_dict[type][0]/nmols)		
		
	labels = ['SFR','Carbon Star','Dark Cloud','Other','LOS Cloud']
	
	fracs = [my_dict[x][1] for x in labels]			

	#set up a plot
	
	ax = plt.subplot(111)
	
	size = 0.1
	
	def getshift(x):
	
		return -(90-(360 - 360*x)/2)
		
	def getper(x):
	
		return '{}' .format(int(100*x)) + r'\%'
	
	ax.pie([fracs[0],1.-fracs[0]],  colors=['dodgerblue','#EEEEEE'], radius=1, startangle=getshift(fracs[0]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[1],1.-fracs[1]],  colors=['darkorange','#EEEEEE'], radius=1-size-.02, startangle=getshift(fracs[1]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[2],1.-fracs[2]],  colors=['forestgreen','#EEEEEE'], radius=1-2*size-.04, startangle=getshift(fracs[2]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[3],1.-fracs[3]],  colors=['violet','#EEEEEE'], radius=1-3*size-.06, startangle=getshift(fracs[3]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[4],1.-fracs[4]],  colors=['red','#EEEEEE'], radius=1-4*size-.08, startangle=getshift(fracs[4]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))

	ax.annotate(r'\textbf{SFR}',xy=(0.5,0.11),xycoords='axes fraction',color='dodgerblue', ha='center',size=14)
	ax.annotate(r'\textbf{Carbon Star}',xy=(0.5,0.16),xycoords='axes fraction',color='darkorange', ha='center',size=14)
	ax.annotate(r'\textbf{Dark Cloud}',xy=(0.5,0.205),xycoords='axes fraction',color='forestgreen', ha='center',size=14)
	ax.annotate(r'\textbf{Other}',xy=(0.5,0.255),xycoords='axes fraction',color='violet', ha='center',size=14)
	ax.annotate(r'\textbf{LOS Cloud}',xy=(0.5,0.305),xycoords='axes fraction',color='red', ha='center',size=14)
	
	percents = ['\\textbf{' + '{:.1f}' .format((x*100)) + '}\%'  for x in fracs]
	
	ax.annotate(percents[4],xy=(0.51,0.68),xycoords='axes fraction',color='white', ha='center',size=12,)
	ax.annotate(percents[3],xy=(0.51,0.725),xycoords='axes fraction',color='white', ha='center',size=12,)
	ax.annotate(percents[2],xy=(0.51,0.775),xycoords='axes fraction',color='white', ha='center',size=12,)
	ax.annotate(percents[1],xy=(0.51,0.825),xycoords='axes fraction',color='white', ha='center',size=12,)
	ax.annotate(percents[0],xy=(0.51,0.87),xycoords='axes fraction',color='white', ha='center',size=12,)
	
	
	plt.tight_layout()
	plt.show()
	
	plt.savefig('source_pie_chart.pdf',format='pdf',transparent=True,bbox_inches='tight',pad_inches=-.65)

	return	
	
def indiv_source_pie_chart(my_list):

	'''
	Makes a pie chart of the fraction of interstellar molecules detected in IRC+10216, TMC-1, Orion, and Sgr
	'''

	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Individual Source Pie Chart')
	
	fig = plt.figure(num='Individual Source Pie Chart',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')
	
	#gather the data

	my_dict = {
	
		'IRC+10216' 	: 	[0],
		'TMC-1'	:	[0],
		'Orion'	:	[0],
		'Sgr B2'	:	[0],
		'Other'		:	[0],		
		}
		
	for mol in my_list:
	
		for source in mol.sources:
		
			other = False
		
			if source == IRC10216:
				my_dict['IRC+10216'][0] += 1
				other=True
			if source == TMC1:
				my_dict['TMC-1'][0] += 1
				other=True
			if source == Orion:
				my_dict['Orion'][0] += 1
				other=True
			if source == SgrB2:
				my_dict['Sgr B2'][0] += 1
				other=True
			if other is False:
				my_dict['Other'][0] += 1		
								
	nmols = len(my_list)
	
	for type in my_dict:
	
		my_dict[type].append(my_dict[type][0]/nmols)		
		
	labels = ['Other','Sgr B2','IRC+10216','TMC-1','Orion']
	
	fracs = [my_dict[x][1] for x in labels]			
	
	print(my_dict)

	#set up a plot
	
	ax = plt.subplot(111)
	
	size = 0.1
	
	def getshift(x):
	
		return -(90-(360 - 360*x)/2)
		
	def getper(x):
	
		return '{}' .format(int(100*x)) + r'\%'
	
	ax.pie([fracs[0],1.-fracs[0]],  colors=['dodgerblue','#EEEEEE'], radius=1, startangle=getshift(fracs[0]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[1],1.-fracs[1]],  colors=['darkorange','#EEEEEE'], radius=1-size-.02, startangle=getshift(fracs[1]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[2],1.-fracs[2]],  colors=['forestgreen','#EEEEEE'], radius=1-2*size-.04, startangle=getshift(fracs[2]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[3],1.-fracs[3]],  colors=['violet','#EEEEEE'], radius=1-3*size-.06, startangle=getshift(fracs[3]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))
	ax.pie([fracs[4],1.-fracs[4]],  colors=['red','#EEEEEE'], radius=1-4*size-.08, startangle=getshift(fracs[4]), wedgeprops=dict(width=size,edgecolor='w',linewidth=1))

	ax.annotate(r'\textbf{Other}',xy=(0.5,0.11),xycoords='axes fraction',color='dodgerblue', ha='center',size=14)
	ax.annotate(r'\textbf{Sgr B2}',xy=(0.5,0.16),xycoords='axes fraction',color='darkorange', ha='center',size=14)
	ax.annotate(r'\textbf{IRC+10216}',xy=(0.5,0.205),xycoords='axes fraction',color='forestgreen', ha='center',size=14)
	ax.annotate(r'\textbf{TMC-1}',xy=(0.5,0.255),xycoords='axes fraction',color='violet', ha='center',size=14)
	ax.annotate(r'\textbf{Orion}',xy=(0.5,0.305),xycoords='axes fraction',color='red', ha='center',size=14)
	
	percents = ['\\textbf{' + '{:.1f}' .format((x*100)) + '}\%'  for x in fracs]
	
	ax.annotate(percents[4],xy=(0.51,0.68),xycoords='axes fraction',color='white', ha='center',size=12,)
	ax.annotate(percents[3],xy=(0.51,0.725),xycoords='axes fraction',color='white', ha='center',size=12,)
	ax.annotate(percents[2],xy=(0.51,0.775),xycoords='axes fraction',color='white', ha='center',size=12,)
	ax.annotate(percents[1],xy=(0.51,0.825),xycoords='axes fraction',color='white', ha='center',size=12,)
	ax.annotate(percents[0],xy=(0.51,0.87),xycoords='axes fraction',color='white', ha='center',size=12,)
	
	
	plt.tight_layout()
	plt.show()
	
	plt.savefig('indiv_source_pie_chart.pdf',format='pdf',transparent=True,bbox_inches='tight',pad_inches=-.65)

	return		

def mol_type_by_source_type(my_list):

	'''
	Generates four pie charts, one for each generalized source type, with the wedges for the types of molecules detected first in each type
	'''
	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Molecule Type by Source Type')
	
	fig,axs = plt.subplots(2,2,num='Molecule Type by Source Type',figsize=(15,12))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':26, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')	

	#collect the data
	
	type_dict = {
	
		'Carbon Star' 	: 	{'Anion': 0, 'Cation': 0, 'Cyclic': 0, 'Neutral': 0, 'Radical': 0},
		'Dark Cloud'	:	{'Anion': 0, 'Cation': 0, 'Cyclic': 0, 'Neutral': 0, 'Radical': 0},
		'LOS Cloud'	:	{'Anion': 0, 'Cation': 0, 'Cyclic': 0, 'Neutral': 0, 'Radical': 0},
		'SFR'	:	{'Anion': 0, 'Cation': 0, 'Cyclic': 0, 'Neutral': 0, 'Radical': 0},
		
		}
	
	for mol in my_list:
	
		#we'll make a dictionary here to flag if we've credited things already
		
		credit_dict = {
	
			'Carbon Star' 	: 	False,
			'Dark Cloud'	:	False,
			'LOS Cloud'	:	False,
			'SFR'	:	False,
			
			}		
	
		for source in mol.sources:
		
			if source.type in type_dict:
			
				if credit_dict[source.type] is False:
				
					if mol.anion is True:
					
						type_dict[source.type]['Anion'] += 1
						
					if mol.cation is True:
					
						type_dict[source.type]['Cation'] += 1

					if mol.cyclic is True:
					
						type_dict[source.type]['Cyclic'] += 1	
						
					if mol.neutral is True:
					
						type_dict[source.type]['Neutral'] += 1
						
					if mol.radical is True:
					
						type_dict[source.type]['Radical'] += 1																

					credit_dict[source.type] = True
					
	#make the pie charts
	
	#Carbon Stars
	
	carbon_data = [
	
		type_dict['Carbon Star']['Anion'], 
		#type_dict['Carbon Star']['Cation'],
		type_dict['Carbon Star']['Cyclic'], 
		type_dict['Carbon Star']['Neutral'], 
		type_dict['Carbon Star']['Radical'],
		
		]
		
	carbon_colors = [
	
		'darkorange',
		#'forestgreen',
		'violet',
		'dodgerblue',
		'red',
	
		]		
		
	#Dark Clouds
	
	dark_data = [
	
		type_dict['Dark Cloud']['Anion'], 
		type_dict['Dark Cloud']['Cation'], 
		type_dict['Dark Cloud']['Cyclic'], 
		type_dict['Dark Cloud']['Neutral'], 
		type_dict['Dark Cloud']['Radical'],
		
		]
		
	dark_colors = [
	
		'darkorange',
		'forestgreen',
		'violet',
		'dodgerblue',
		'red',
	
		]		
		
	#LOS Clouds
	
	los_data = [
	
		#type_dict['LOS Cloud']['Anion'], 
		type_dict['LOS Cloud']['Cation'], 
		#type_dict['LOS Cloud']['Cyclic'], 
		type_dict['LOS Cloud']['Neutral'], 
		type_dict['LOS Cloud']['Radical'],
		
		]
		
	los_colors = [
	
		#'darkorange',
		'forestgreen',
		#'violet',
		'dodgerblue',
		'red',
	
		]		
		
	#SFRs
	
	sfr_data = [
	
		#type_dict['SFR']['Anion'], 
		type_dict['SFR']['Cation'], 
		type_dict['SFR']['Cyclic'], 
		type_dict['SFR']['Neutral'], 
		type_dict['SFR']['Radical'],
		
		]		
						
	sfr_colors = [
	
		#'darkorange',
		'forestgreen',
		'violet',
		'dodgerblue',
		'red',
	
		]
	
	types = [
	
		'Anion',
		'Cation',
		'Cyclic',
		'Neutral',
		'Radical',
	
		]		
	
	axs[0,0].pie(carbon_data,labels=carbon_data,colors=carbon_colors, labeldistance=0.8, wedgeprops = {'linewidth' : 1.0, 'edgecolor' : 'black', 'alpha' : 0.5})	
	axs[0,0].annotate(r'\textbf{Carbon Stars}',xy=(0.5,1.0),xycoords='axes fraction',color='black', ha='center',va='top',size=30)
	
	wedges, texts = axs[0,1].pie(dark_data,labels=dark_data,colors=dark_colors, labeldistance=0.8, wedgeprops = {'linewidth' : 1.0, 'edgecolor' : 'black', 'alpha' : 0.5})	
	axs[0,1].annotate(r'\textbf{Dark Clouds}',xy=(0.5,1.0),xycoords='axes fraction',color='black', ha='center',va='top',size=30)

	axs[1,0].pie(los_data,labels=los_data,colors=los_colors, labeldistance=0.8, wedgeprops = {'linewidth' : 1.0, 'edgecolor' : 'black', 'alpha' : 0.5})	
	axs[1,0].annotate(r'\textbf{LOS Clouds}',xy=(0.5,1.0),xycoords='axes fraction',color='black', ha='center',va='top',size=30)
	
	axs[1,1].pie(sfr_data,labels=sfr_data,colors=sfr_colors, labeldistance=0.8, wedgeprops = {'linewidth' : 1.0, 'edgecolor' : 'black', 'alpha' : 0.5})		
	axs[1,1].annotate(r'\textbf{SFRs}',xy=(0.5,1.0),xycoords='axes fraction',color='black', ha='center',va='top',size=30)
	
	
	axs[0,1].legend(wedges, types,
          title="Molecule Types",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
	
	fig.tight_layout()
	
	fig.subplots_adjust(wspace=-0.3, hspace=0)
	
	plt.show()
	
	plt.savefig('mol_type_by_source_type.pdf',format='pdf',transparent=True,bbox_inches='tight')			

	return 	
	
def du_by_source_type(my_list):

	'''
	Makes a KDE plot of the dus in each source type
	'''
	
	#gather the data

	my_dict = {
	
		'Carbon Star' 	: 	[],
		'Dark Cloud'	:	[],
		'LOS Cloud'	:	[],
		'SFR'	:	[],	
		}
		
	#we have to be a little careful here, because for a given source, there can be two SFRs listed, and we only want to credit it once
		
	for mol in my_list:
	
		#we'll make a dictionary here to flag if we've credited things already
		
		credit_dict = {
	
			'Carbon Star' 	: 	False,
			'Dark Cloud'	:	False,
			'LOS Cloud'	:	False,
			'SFR'	:	False,	
			}		
	
		if mol.du is None or mol.fullerene is True:
		
			continue
		
		for source in mol.sources:
		
			if source.type in my_dict:
			
				if credit_dict[source.type] is False:
			
					my_dict[source.type].append(mol.du)
					
					credit_dict[source.type] = True
	
	plt.close('DU by Source Type')
	
	fig = plt.figure(num='DU by Source Type',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')		

	#load up an axis
	
	ax = fig.add_subplot(111)

	ax.tick_params(axis='x', which='both', direction='in',length=5,width=1)
	ax.tick_params(axis='y', which='both', direction='in',length=5,width=1)	

	#label the axes
	
	plt.xlabel('Degree of Unsaturation')
	plt.ylabel('Kernel Density Estimate')
	
	#Do the estimates and plot them
	
	xvals = np.arange(0,15,0.1)
	
	density_carbon = gkde(my_dict['Carbon Star'])
	density_carbon.covariance_factor = lambda : 0.5
	density_carbon._compute_covariance()
	
	density_dark = gkde(my_dict['Dark Cloud'])
	density_dark.covariance_factor = lambda : 0.5
	density_dark._compute_covariance()
	
	density_los = gkde(my_dict['LOS Cloud'])
	density_los.covariance_factor = lambda : 0.5
	density_los._compute_covariance()	
	
	density_sfr = gkde(my_dict['SFR'])
	density_sfr.covariance_factor = lambda : 0.5
	density_sfr._compute_covariance()	
	
	x_ann = 0.97
	y_ann = 0.96
	y_sep = 0.06	
		
	ax.plot(xvals,density_carbon(xvals),color='darkorange')
	ax.fill_between(xvals,density_carbon(xvals),0,facecolor='darkorange',alpha=0.25,zorder=4)

	ax.plot(xvals,density_dark(xvals),color='forestgreen')
	ax.fill_between(xvals,density_dark(xvals),0,facecolor='forestgreen',alpha=0.25,zorder=4)

	ax.plot(xvals,density_los(xvals),color='red')
	ax.fill_between(xvals,density_los(xvals),0,facecolor='red',alpha=0.25,zorder=4)

	ax.plot(xvals,density_sfr(xvals),color='dodgerblue')
	ax.fill_between(xvals,density_sfr(xvals),0,facecolor='dodgerblue',alpha=0.25,zorder=4)
		
	ax.annotate(r'\underline{Source Types}',xy=(x_ann,y_ann-0*y_sep),xycoords='axes fraction',color='black',ha='right',va='top')

	ax.annotate('SFR',xy=(x_ann,y_ann-y_sep),xycoords='axes fraction',color='dodgerblue',ha='right',va='top')
	ax.annotate('{}' .format(len(my_dict['SFR'])), xy=(2.7,0.3),xycoords='data',ha='left',va='bottom',color='dodgerblue')	

	ax.annotate('Carbon Star',xy=(x_ann,y_ann-2*y_sep),xycoords='axes fraction',color='darkorange',ha='right',va='top')
	ax.annotate('{}' .format(len(my_dict['Carbon Star'])), xy=(6,0.14),xycoords='data',ha='left',va='bottom',color='darkorange')	

	ax.annotate('Dark Cloud',xy=(x_ann,y_ann-3*y_sep),xycoords='axes fraction',color='forestgreen',ha='right',va='top')
	ax.annotate('{}' .format(len(my_dict['Dark Cloud'])), xy=(10.9,0.02),xycoords='data',ha='left',va='bottom',color='forestgreen')	

	ax.annotate('LOS Cloud',xy=(x_ann,y_ann-4*y_sep),xycoords='axes fraction',color='red',ha='right',va='top')	
	ax.annotate('{}' .format(len(my_dict['LOS Cloud'])), xy=(2.3,0.4),xycoords='data',ha='left',va='bottom',color='red')	

	
	plt.show()
	
	plt.savefig('du_by_source_type_kde.pdf',format='pdf',transparent=True,bbox_inches='tight',pad_inches=0)
	
	return			

def rel_du_by_source_type(my_list):

	'''
	Makes a KDE plot of the relative dus in each source type
	'''
	
	#gather the data

	my_dict = {
	
		'Carbon Star' 	: 	[],
		'Dark Cloud'	:	[],
		'LOS Cloud'	:	[],
		'SFR'	:	[],	
		}
		
	#we have to be a little careful here, because for a given source, there can be two SFRs listed, and we only want to credit it once
		
	for mol in my_list:
	
		#we'll make a dictionary here to flag if we've credited things already
		
		credit_dict = {
	
			'Carbon Star' 	: 	False,
			'Dark Cloud'	:	False,
			'LOS Cloud'	:	False,
			'SFR'	:	False,	
			}		
	
		if mol.du is None or mol.fullerene is True:
		
			continue
		
		for source in mol.sources:
		
			if source.type in my_dict:
			
				if credit_dict[source.type] is False:
			
					my_dict[source.type].append(mol.du/mol.maxdu)
					
					credit_dict[source.type] = True
	
	plt.close('Relative DU by Source Type')
	
	fig,axs = plt.subplots(2,2,num='Relative DU by Source Type',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':18, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')		

	#axes ticks and limits
	
	for ax in axs.flat:
		ax.set_ylim([0,3.75])
		ax.tick_params(axis='y', which='both', direction='in',length=5,width=1)
		ax.tick_params(axis='x', which='both', direction='in',length=5,width=1)
		
	axs[0,0].set_yticklabels(axs[0,0].get_yticklabels(), visible=False)			
	axs[0,0].set_xticklabels(axs[0,0].get_xticklabels(), visible=False)	

	axs[1,1].set_yticklabels(axs[1,0].get_yticklabels(), visible=False)	
	axs[1,1].set_xticklabels(axs[1,0].get_xticklabels(), visible=False)	
	
	axs[0,1].set_yticklabels(axs[1,0].get_yticklabels(), visible=False)	
	axs[0,1].set_xticklabels(axs[1,0].get_xticklabels(), visible=False)				
		

	#label the axes
	
	axs[1,0].set(xlabel = 'Relative Degree of Unsaturation')
	axs[1,0].set(ylabel = 'Kernel Density Estimate')	
	
	#Do the estimates and plot them
	
	xvals = np.arange(0,1,0.01)
	
	density_carbon = gkde(my_dict['Carbon Star'])
	density_carbon.covariance_factor = lambda : 0.5
	density_carbon._compute_covariance()
	
	density_dark = gkde(my_dict['Dark Cloud'])
	density_dark.covariance_factor = lambda : 0.5
	density_dark._compute_covariance()
	
	density_los = gkde(my_dict['LOS Cloud'])
	density_los.covariance_factor = lambda : 0.5
	density_los._compute_covariance()	
	
	density_sfr = gkde(my_dict['SFR'])
	density_sfr.covariance_factor = lambda : 0.5
	density_sfr._compute_covariance()	
	
	x_ann = 0.97
	y_ann = 0.96
	y_sep = 0.06	
		
	axs[0,0].plot(xvals,density_carbon(xvals),color='darkorange')
	axs[0,0].fill_between(xvals,density_carbon(xvals),0,facecolor='darkorange',alpha=0.25,zorder=4)
	axs[0,0].annotate('Carbon Star',xy=[0.04,0.96],xycoords='axes fraction',ha='left',va='top',size=24,color='darkorange')

	axs[1,0].plot(xvals,density_dark(xvals),color='forestgreen')
	axs[1,0].fill_between(xvals,density_dark(xvals),0,facecolor='forestgreen',alpha=0.25,zorder=4)
	axs[1,0].annotate('Dark Cloud',xy=[0.04,0.96],xycoords='axes fraction',ha='left',va='top',size=24,color='forestgreen')


	axs[0,1].plot(xvals,density_los(xvals),color='red')
	axs[0,1].fill_between(xvals,density_los(xvals),0,facecolor='red',alpha=0.25,zorder=4)
	axs[0,1].annotate('LOS Cloud',xy=[0.04,0.96],xycoords='axes fraction',ha='left',va='top',size=24,color='red')


	axs[1,1].plot(xvals,density_sfr(xvals),color='dodgerblue')
	axs[1,1].fill_between(xvals,density_sfr(xvals),0,facecolor='dodgerblue',alpha=0.25,zorder=4)
	axs[1,1].annotate('SFR',xy=[0.04,0.96],xycoords='axes fraction',ha='left',va='top',size=24,color='dodgerblue')

	
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.show()
	
	plt.savefig('relative_du_by_source_type_kde.pdf',format='pdf',transparent=True,bbox_inches='tight',pad_inches=0)
	
	return		
	
def mass_by_source_type(my_list):

	'''
	Makes a KDE plot of the masses in each source type
	'''
	
	#gather the data

	my_dict = {
	
		'Carbon Star' 	: 	[],
		'Dark Cloud'	:	[],
		'LOS Cloud'	:	[],
		'SFR'	:	[],	
		}
		
	#we have to be a little careful here, because for a given source, there can be two SFRs listed, and we only want to credit it once
		
	masses = []
	
	for mol in my_list:
		
		#we'll make a dictionary here to flag if we've credited things already
		
		credit_dict = {
	
			'Carbon Star' 	: 	False,
			'Dark Cloud'	:	False,
			'LOS Cloud'	:	False,
			'SFR'	:	False,	
			}		
	
		#drop the fullerenes
		
		if mol.fullerene is True:
		
			continue
			
		#add the mass to the list for axis limit purposes
		
		masses.append(mol.mass)	
		
		for source in mol.sources:
		
			if source.type in my_dict:
			
				if credit_dict[source.type] is False:
			
					my_dict[source.type].append(mol.mass)
					
					credit_dict[source.type] = True
	
	plt.close('Mass by Source Type')
	
	fig = plt.figure(num='Mass by Source Type',figsize=(10,8))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':24, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')		

	#axes ticks and limits
	
	ax = plt.subplot(111)
	
	ax.set_xlim([min(masses),max(masses)])
	ax.tick_params(axis='y', which='both', direction='in',length=5,width=1,labelleft='off')
	ax.tick_params(axis='x', which='both', direction='in',length=5,width=1,labelbottom='off')		

	#label the axes
	
	ax.set(xlabel = 'Molecular Mass (amu)')
	ax.set(ylabel = 'Kernel Density Estimate')	
	
	#Do the estimates and plot them
	
	xvals = np.arange(0,max(masses),1)
	
	density_carbon = gkde(my_dict['Carbon Star'])
	density_carbon.covariance_factor = lambda : 0.5
	density_carbon._compute_covariance()
	
	density_dark = gkde(my_dict['Dark Cloud'])
	density_dark.covariance_factor = lambda : 0.5
	density_dark._compute_covariance()
	
	density_los = gkde(my_dict['LOS Cloud'])
	density_los.covariance_factor = lambda : 0.5
	density_los._compute_covariance()	
	
	density_sfr = gkde(my_dict['SFR'])
	density_sfr.covariance_factor = lambda : 0.5
	density_sfr._compute_covariance()	
		
	ax.plot(xvals[min(masses):max(masses)],density_carbon(xvals[min(masses):max(masses)]),color='darkorange')
	ax.fill_between(xvals[min(masses):max(masses)],density_carbon(xvals[min(masses):max(masses)]),0,facecolor='darkorange',alpha=0.25,zorder=4)
	
	ax.plot(xvals[min(masses):max(masses)],density_dark(xvals[min(masses):max(masses)]),color='forestgreen')
	ax.fill_between(xvals[min(masses):max(masses)],density_dark(xvals[min(masses):max(masses)]),0,facecolor='forestgreen',alpha=0.25,zorder=4)
	
	ax.plot(xvals[min(masses):max(masses)],density_los(xvals[min(masses):max(masses)]),color='red')
	ax.fill_between(xvals[min(masses):max(masses)],density_los(xvals[min(masses):max(masses)]),0,facecolor='red',alpha=0.25,zorder=4)

	ax.plot(xvals[min(masses):max(masses)],density_sfr(xvals[min(masses):max(masses)]),color='dodgerblue')
	ax.fill_between(xvals[min(masses):max(masses)],density_sfr(xvals[min(masses):max(masses)]),0,facecolor='dodgerblue',alpha=0.25,zorder=4)
	
	x_ann = 0.97
	y_ann = 0.96
	y_sep = 0.06	
		
	ax.annotate(r'\underline{Source Types}',xy=(x_ann,y_ann-0*y_sep),xycoords='axes fraction',color='black',ha='right',va='top')
	ax.annotate('SFR',xy=(x_ann,y_ann-y_sep),xycoords='axes fraction',color='dodgerblue',ha='right',va='top')
	ax.annotate('Carbon Star',xy=(x_ann,y_ann-2*y_sep),xycoords='axes fraction',color='darkorange',ha='right',va='top')	
	ax.annotate('Dark Cloud',xy=(x_ann,y_ann-3*y_sep),xycoords='axes fraction',color='forestgreen',ha='right',va='top')
	ax.annotate('LOS Cloud',xy=(x_ann,y_ann-4*y_sep),xycoords='axes fraction',color='red',ha='right',va='top')	

	plt.show()
	
	plt.savefig('mass_by_source_type_kde.pdf',format='pdf',transparent=True,bbox_inches='tight',pad_inches=0)
	
	return		

def waves_by_source_type(my_list):

	'''
	Generates four pie charts, one for each generalized source type, with the wedges for the wavelengths used for first detections in those sources
	'''
	#Close an old figure if it exists and initialize a new figure
	
	plt.close('Wavelength by Source Type')
	
	fig,axs = plt.subplots(2,2,num='Wavelength by Source Type',figsize=(15,12))
	
	plt.ion()
	
	#set some font defaults
	
	fontparams = {'size':26, 'family':'sans-serif','sans-serif':['Helvetica']}
	
	plt.rc('font',**fontparams)
	plt.rc('mathtext', fontset='stixsans')	

	#collect the data
	
	type_dict = {
	
		'Carbon Star' 	: 	{'cm': 0, 'mm': 0, 'sub-mm': 0, 'IR': 0, 'UV': 0, 'Vis': 0},
		'Dark Cloud'	:	{'cm': 0, 'mm': 0, 'sub-mm': 0, 'IR': 0, 'UV': 0, 'Vis': 0},
		'LOS Cloud'		:	{'cm': 0, 'mm': 0, 'sub-mm': 0, 'IR': 0, 'UV': 0, 'Vis': 0},
		'SFR'			:	{'cm': 0, 'mm': 0, 'sub-mm': 0, 'IR': 0, 'UV': 0, 'Vis': 0},
		
		}
	
	for mol in my_list:
	
		#we'll make a dictionary here to flag if we've credited things already
		
		credit_dict = {
	
			'Carbon Star' 	: 	False,
			'Dark Cloud'	:	False,
			'LOS Cloud'		:	False,
			'SFR'			:	False,
			
			}		
	
		for source in mol.sources:
		
			if source.type in type_dict:
			
				if credit_dict[source.type] is False:
				
					for wave in mol.wavelengths:
					
						type_dict[source.type][wave] += 1													

					credit_dict[source.type] = True
					
	#make the pie charts
	
	#Carbon Stars
	
	carbon_data = [
	
		type_dict['Carbon Star']['cm'], 
		type_dict['Carbon Star']['mm'],
		type_dict['Carbon Star']['sub-mm'], 
		type_dict['Carbon Star']['IR'], 
		#type_dict['Carbon Star']['UV'] + type_dict['Carbon Star']['Vis'],
		]
		
	carbon_colors = [
	
		'dodgerblue',
		'darkorange',
		'forestgreen',
		'black',
		#'violet',
	
		]		
		
	#Dark Clouds
	
	dark_data = [
	
		type_dict['Dark Cloud']['cm'], 
		type_dict['Dark Cloud']['mm'],
		#type_dict['Dark Cloud']['sub-mm'], 
		#type_dict['Dark Cloud']['IR'], 
		#type_dict['Dark Cloud']['UV'] + type_dict['Dark Cloud']['Vis'],
		]
		
	dark_colors = [
	
		'dodgerblue',
		'darkorange',
		#'forestgreen',
		#'black',
		#'violet',
	
		]			
		
	#LOS Clouds
	
	los_data = [
	
		type_dict['LOS Cloud']['cm'], 
		type_dict['LOS Cloud']['mm'],
		type_dict['LOS Cloud']['sub-mm'], 
		type_dict['LOS Cloud']['IR'], 
		type_dict['LOS Cloud']['UV'] + type_dict['LOS Cloud']['Vis'],
		]
		
	los_colors = [
	
		'dodgerblue',
		'darkorange',
		'forestgreen',
		'black',
		'violet',
	
		]
		
	#SFRs
	
	sfr_data = [
	
		type_dict['SFR']['cm'], 
		type_dict['SFR']['mm'],
		type_dict['SFR']['sub-mm'], 
		#type_dict['SFR']['IR'], 
		#type_dict['SFR']['UV'] + type_dict['SFR']['Vis'],
		]
		
	sfr_colors = [
	
		'dodgerblue',
		'darkorange',
		'forestgreen',
		#'black',
		#'violet',
	
		]
	
	types = [
	
		'cm',
		'mm',
		'sub-mm',
		'IR',
		'UV/Vis',
	
		]
		
	def make_labels(list):
	
		new_labels = []
		
		total = sum(list)
		
		for i in list:
		
			percent = 100*i/total
			
			new_labels.append('{:.1f}\%' .format(percent))
			
		return new_labels			
	
	axs[0,0].pie(carbon_data,labels=make_labels(carbon_data), colors=carbon_colors, labeldistance=1.1, wedgeprops = {'linewidth' : 1.0, 'edgecolor' : 'black', 'alpha' : 0.5})	
	axs[0,0].annotate(r'\textbf{Carbon Stars}',xy=(0.5,1.05),xycoords='axes fraction',color='black', ha='center',va='top',size=30)
	
	axs[0,1].pie(dark_data,labels=make_labels(dark_data),colors=dark_colors, labeldistance=1.1, wedgeprops = {'linewidth' : 1.0, 'edgecolor' : 'black', 'alpha' : 0.5})	
	axs[0,1].annotate(r'\textbf{Dark Clouds}',xy=(0.5,1.05),xycoords='axes fraction',color='black', ha='center',va='top',size=30)

	wedges, texts = axs[1,0].pie(los_data,labels=make_labels(los_data),colors=los_colors, labeldistance=1.1, wedgeprops = {'linewidth' : 1.0, 'edgecolor' : 'black', 'alpha' : 0.5})	
	axs[1,0].annotate(r'\textbf{LOS Clouds}',xy=(0.5,1.05),xycoords='axes fraction',color='black', ha='center',va='top',size=30)
	
	axs[1,1].pie(sfr_data,labels=make_labels(sfr_data),colors=sfr_colors, labeldistance=1.1, wedgeprops = {'linewidth' : 1.0, 'edgecolor' : 'black', 'alpha' : 0.5})		
	axs[1,1].annotate(r'\textbf{SFRs}',xy=(0.5,1.05),xycoords='axes fraction',color='black', ha='center',va='top',size=30)
	
	
	axs[0,1].legend(wedges, types,
          title="Wavelengths",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
	
	fig.tight_layout()
	
	fig.subplots_adjust(wspace=0.1, hspace=0.1)
	
	plt.show()
	
	plt.savefig('waves_by_source_type.pdf',format='pdf',transparent=True,bbox_inches='tight')			

	return 									