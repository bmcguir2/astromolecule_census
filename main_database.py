#!/usr/bin/env python

#############################################################
#						Revision History					#
#############################################################

# generates and writes out a ton of things

# 1.0 - Project start (5/12/2018)
# 1.1 - Cleaned up code
# 1.2 - Added sulfur compounds to DU calculation.

#############################################################
#							Preamble						#
#############################################################

'''

I wrote this code primarily to help me easily generate ascii files for making the plots 
found in the main text (2018 Census of Interstellar [...] Molecules).  I have
attempted to comment the code to some degree to make it more readable and potentially
useful for others.  The utility functions at the end of the code do exactly what I needed
them to do, but may not be optimal python, and indeed depending on when in the process
I wrote each function, some of them may do the same thing in completely different ways.
Depending on when you got this code, some functions may be in development, or may actually
have been found to be in error but not yet fixed.  Moral of the story: read it first =).

A few critical notes for using:

1) I wrote this code intending for it to be loaded and used in an interactive python
environment (I use ipython).  I don't know how it will behave outside of these.

2) It is in Python 3.  I won't be porting it/making it backwards compatible to Python 2.7

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

'''

import os, sys, argparse, math
from operator import itemgetter
from math import ceil

#Python version check

if sys.version_info.major != 3:

	print("This code is written in Python 3.  It will not execute in Python 2.7.  Exiting (sorry).")
	
	quit()

version = 1.2


#############################################################
#						Molecule Class 						#
#############################################################

class Molecule(object):

	def __init__(self,name,formula,year,label,sources,telescopes,wavelengths,other_names='',neutral=False,cation=False,anion=False,radical=False,cyclic=False,mass=0,du=0,natoms=0,H=0,C=0,O=0,N=0,S=0,P=0,Si=0,Cl=0,F=0,Mg=0,Na=0,Al=0,K=0,Fe=0,Ti=0,Ar=0,d_ref=None,lab_ref=None,notes=None,ice=False,ice_d_ref=None,ice_l_ref=None,ppd=None,exgal=None,exo=None,isos=None,isomers=None,ppd_isos=None,ppd_d_ref=None,ppd_l_ref=None,ppd_isos_ref=None,exgal_d_ref=None,exgal_l_ref=None,exo_d_ref=None,exo_l_ref=None,exgal_sources=None,isos_d_ref=None,isos_l_ref=None):
	
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
		self.mass = mass
		self.du = du
		self.natoms = natoms
		self.H = H
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
		
class Source(object):

	def __init__(self,name,type=None,ra=None,dec=None,detects=0,simbad_url=None):
	
		self.name = name
		self.type = type
		self.ra = ra
		self.dec = dec
		self.detects = detects
		self.simbad_url = simbad_url
		
#############################################################
#						Two Atom Molecules					#
#############################################################		

CH = Molecule('methylidyne','CH',1937,'CH','LOS Cloud','Mt Wilson','UV, Vis',neutral=True,H=1,C=1,d_ref='Dunham 1937 PASP 49, 26; Swings & Rosenfeld 1937 ApJ 86, 483; McKellar 1940 PASP 52, 187',lab_ref='Jevons 1932 Phys Soc. pp 177-179; Brazier & Brown 1983 JCP 78, 1608',notes='*First radio in Rydbeck et al. 1973 Nature 246, 466',exgal=True,exgal_d_ref='Whiteoak et al. 1980 MNRAS 190, 17',exgal_sources='LMC, NGC 4945, NGC 5128')
CN = Molecule('cyano radical','CN',1940,'CN','LOS Cloud','Mt Wilson','UV',radical=True,neutral=True,C=1,N=1,d_ref='McKellar 1940 PASP 52, 187',lab_ref='Poletto and Rigutti 1965 Il Nuovo Cimento 39, 519; Dixon & Woods 1977 JCP 67, 3956; Thomas & Dalby 1968 Can. J. Phys. 46, 2815',notes='*First radio in Jefferts et al. 1970 ApJ 161, L87',ppd=True,ppd_d_ref='Kastner et al. 1997 Science 277, 67; Dutrey et al. 1997 A&A 317, L55',ppd_isos='C15N',ppd_isos_ref='[C15N] Hily-Blant et al. 2017 A&A 603, L6',exgal=True,exgal_d_ref='Henkel et al. 1988 A&A 201, L23',exgal_sources='M82, NGC 253, IC 342')
CHp = Molecule('methylidyne cation','CH+',1941,'CH+','LOS Cloud','Mt Wilson','UV, Vis',cation=True,H=1,C=1,d_ref='Douglas & Herzberg 1941 ApJ 94, 381; Dunham 1937 PASP 49, 26',lab_ref='Douglas & Herzberg 1941 ApJ 94, 381',notes=None,ppd=True,ppd_d_ref='Thi et al. 2011 A&A 530, L2',exgal=True,exgal_d_ref='Magain & Gillet 1987 A&A 184, L5',exgal_sources='LMC')
OH = Molecule('hydroxyl radical','OH',1963,'OH','Cas A LOS','Millstone Hill','cm',radical=True,neutral=True,H=1,O=1,d_ref='Weinreb et al. 1963 Nature 200, 829',lab_ref='Ehrenstein et al. 1959 PRL 3, 40',notes=None,ppd=True,ppd_d_ref='Mandell et al. 2008 ApJ 681, L25; Salyk et al. 2008 ApJ 676, L49',exgal=True,exgal_d_ref='Weliachew 1971 ApJ 167, L47',exgal_sources='M82, NGC 253')
CO = Molecule('carbon monoxide','CO',1970,'CO','Orion','NRAO 36-ft','mm',neutral=True,C=1,O=1,d_ref='Wilson et al. 1970 ApJ 161, L43',lab_ref='Cord et al. 1968 Microwave Spectral Tables V5',notes=None,ice=True,ice_d_ref='Soifer et al. 1979 ApJ 232, L53',ice_l_ref='Mantz et al. 1975 JMS 57, 155',ppd=True,ppd_d_ref='Beckwith et al. 1986 ApJ 309, 755',ppd_isos='13CO, C18O, C17O',ppd_isos_ref='[13CO] Sargent & Beckwith 1987 ApJ 323, 294 [C18O] Dutrey et al. 1994 A&A 286, 149 [C17O] Smith et al. 2009 ApJ 701, 163; Guilloteau et al. 2013 A&A 549, A92',exo=True,exo_d_ref='Madhusudhan et al. 2011 Nature 469, 64; Barman et al. 2011 ApJ 733, 65; Lanotte et al. 2014 A&A 572, A73; Barman et al. 2015 ApJ 804, 61',exgal=True,exgal_d_ref='Rickard et al. 1975 ApJ 199, L75',exgal_sources='M82, NGC 253')
H2 = Molecule('hydrogen','H2',1970,'H2','Xi Per LOS','Aerobee-150 Rocket','UV',neutral=True,H=2,d_ref='Carruthers 1970 ApJ 161, L81',lab_ref='Carruthers 1970 ApJ 161, L81',notes=None,ppd=True,ppd_d_ref='Thi et al. 1999 ApJ 521, L63',ppd_isos='HD',ppd_isos_ref='[HD] Bergin et al. 2013 Nature 493, 644',exgal=True,exgal_d_ref='Thompson et al. 1978 ApJ 222, L49',exgal_sources='NGC 1068')
SiO = Molecule('silicon monoxide','SiO',1971,'SiO','Sgr B2','NRAO 36-ft','mm',neutral=True,O=1,Si=1,d_ref='Wilson et al. 1971 ApJ 167, L97',lab_ref='Törring 1968 Z. Naturforschung 23A, 777; Raymonda et al. 1970 JCP 52, 3458',notes=None,exgal=True,exgal_d_ref='Mauersberger & Henkel 1991 A&A 245, 457',exgal_sources='NGC 253')
CS = Molecule('carbon monosulfide','CS',1971,'CS','Orion, W51, IRC+10216, DR 21','NRAO 36-ft','mm',neutral=True,C=1,S=1,d_ref='Penzias et al. 1971 ApJ 168, L53',lab_ref='Mockler & Bird 1955 Phys Rev 98, 1837',notes=None,ppd=True,ppd_d_ref='Ohashi et al. 1991 AJ 102, 2054; Blake et al. 1992 ApJ 391, L99; Guilloteau et al. 2012 A&A 548, A70',ppd_isos='C34S',ppd_isos_ref='[C34S] Artur de la Villarmois et al. 2018 A&A 614, A26',exgal=True,exgal_d_ref='Henkel & Bally 1985 A&A 150, L25',exgal_sources='M82, IC 342')
SO = Molecule('sulfur monoxide','SO',1973,'SO','Orion','NRAO 36-ft','mm',neutral=True,O=1,S=1,d_ref='Gottlieb & Ball 1973 ApJ 184, L59',lab_ref='Winnewisser et al. 1964 JCP 41, 1687',notes=None,ppd=True,ppd_d_ref='Fuente et al. 2010 A&A 524, A19',exgal=True,exgal_d_ref='Johansson 1991 Proc. IAU Symposium 146, 1; Petuchowski & Bennett 1992 ApJ 391, 137',exgal_sources='M82, NGC 253')
SiS = Molecule('silicon monosulfide','SiS',1975,'SiS','IRC+10216','NRAO 36-ft','mm',neutral=True,S=1,Si=1,d_ref='Morris et al. 1975 ApJ 199, L47',lab_ref='Hoeft 1965 Z. fur Naturforschung A, A20, 1327',notes=None)
NS = Molecule('nitrogen monosulfide','NS',1975,'NS','Sgr B2','NRAO 36-ft','mm',neutral=True,N=1,S=1,d_ref='Gottlieb et al. 1975 ApJ 200, L147; Kuiper et al. 1975 ApJ 200, L151',lab_ref='Amano et al. 1969 JMS 32, 97',notes=None,exgal=True,exgal_d_ref='Martin et al. 2003 A&A 411, L465',exgal_sources='NGC 253')
C2 = Molecule('dicarbon','C2',1977,'C2','Cygnus OB2 - 12 LOS','Mt Hopkins','IR',neutral=True,C=2,d_ref='Souza and Lutz 1977 ApJ 216, L49',lab_ref='Phillips 1948 ApJ 107, 389',exgal=True,exgal_d_ref='Welty et al. 2012 MNRAS 428, 1107', exgal_sources='SMC', notes=None)
NO = Molecule('nitric oxide','NO',1978,'NO','Sgr B2','NRAO 36-ft','mm',neutral=True,O=1,N=1,d_ref='Liszt and Turner 1978 ApJ 224, L73',lab_ref='Gallagher & Johnson 1956 Phys Rev 103, 1727',notes=None,exgal=True,exgal_d_ref='Martin et al. 2003 A&A 411, L465',exgal_sources='NGC 253')
HCl = Molecule('hydrogen chloride','HCl',1985,'HCl','Orion','Kuiper','sub-mm',neutral=True,H=1,Cl=1,d_ref='Blake et al. 1985 ApJ 295, 501',lab_ref='de Lucia et al. 1971 Phys Rev A 3, 1849',notes=None)
NaCl = Molecule('sodium chloride','NaCl',1987,'NaCl','IRC+10216','IRAM','mm',neutral=True,Cl=1,Na=1,d_ref='Cernicharo & Guélin 1987 A&A 183, L10',lab_ref='Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609',notes=None)
AlCl = Molecule('aluminum chloride','AlCl',1987,'AlCl','IRC+10216','IRAM','mm',neutral=True,Cl=1,Al=1,d_ref='Cernicharo & Guélin 1987 A&A 183, L10',lab_ref='Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609',notes=None)
KCl = Molecule('potassium chloride','KCl',1987,'KCl','IRC+10216','IRAM','mm',neutral=True,Cl=1,K=1,d_ref='Cernicharo & Guélin 1987 A&A 183, L10',lab_ref='Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609',notes=None)
AlF = Molecule('aluminum fluoride','AlF',1987,'AlF','IRC+10216','IRAM','mm',neutral=True,F=1,Al=1,d_ref='Cernicharo & Guélin 1987 A&A 183, L10',lab_ref='Lovas & Tiemann 1974 J Phys Chem Ref Data 3, 609',notes='*Confirmed in 1994 ApJ 433, 729',isos='26AlF',isos_d_ref='[26AlF] Kamiński et al. 2018 Nature Astronomy 262, 742')
PN = Molecule('phosphorous mononitride','PN',1987,'PN','TMC-1, Orion, W51, Sgr B2','NRAO 12-m, FCRAO 14-m, OVRO','mm',neutral=True,N=1,P=1,d_ref='Sutton et al. 1985 ApJS 58, 341',lab_ref='Wyse et al. 1972 JCP 57, 1106',notes='*Confirmed in Turner & Bally 1987 ApJ 321, L75 and Ziurys 1987 ApJ 321 L81')
SiC = Molecule('silicon carbide','SiC',1989,'SiC','IRC+10216','IRAM','mm',radical=True,neutral=True,C=1,Si=1,d_ref='Cernicharo et al. 1989 ApJ 341, L25',lab_ref='Cernicharo et al. 1989 ApJ 341, L25',notes=None)
CP = Molecule('carbon monophosphide','CP',1990,'CP','IRC+10216','IRAM','mm',radical=True,neutral=True,C=1,P=1,d_ref='Guélin et al. 1990 A&A 230, L9',lab_ref='Saito et al. 1989 ApJ 341, 1114',notes=None)
NH = Molecule('imidogen radical','NH',1991,'NH','Xi Per LOS, HD 27778 LOS','KPNO 4-m, IRAM','UV',radical=True,neutral=True,H=1,N=1,d_ref='Meyer & Roth 1991 ApJ 376, L49',lab_ref='Dixon 1959 Can J. Phys. 37, 1171 and Klaus et al. 1997 A&A 322, L1',notes='*First radio in Cernicharo et al. 2000 ApJ 534, L199',exgal=True,exgal_d_ref='Gonzalez-Alfonso et al. 2004 ApJ 613, 247',exgal_sources='Arp 220')
SiN = Molecule('silicon nitride ','SiN',1992,'SiN','IRC+10216','NRAO 12-m','mm',radical=True,neutral=True,N=1,Si=1,d_ref='Turner 1992 ApJ 388, L35',lab_ref='Saito et al. 1983 JCP 78, 6447',notes=None)
SOp = Molecule('sulfur monoxide cation','SO+',1992,'SO+','IC 443G','NRAO 12-m','mm',cation=True,radical=True,O=1,S=1,d_ref='Turner 1992 ApJ 396, L107',lab_ref='Amano et al. 1991 JMS 146, 519',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS')
COp = Molecule('carbon monoxide cation','CO+',1993,'CO+','M17SW, NGC 7027','NRAO 12-m','mm',cation=True,C=1,O=1,d_ref='Latter et al. 1993 ApJ 419, L97',lab_ref='Sastry et al. 1981 ApJ 250, L91',notes=None,exgal=True,exgal_d_ref='Fuente et al. 2006 ApJ 641, L105',exgal_sources='M82')
HF = Molecule('hydrogen fluoride','HF',1997,'HF','Sgr B2 LOS','ISO','IR',neutral=True,H=1,F=1,d_ref='Neufeld et al. 1997 ApJ 488, L141',lab_ref='Nolt et al. 1987 JMS 125, 274',notes=None,exgal=True,exgal_d_ref='van der Werf et al. 2010 A&A 518, L42; Rangwala et al. 2011 ApJ 743, 94; Monje et al. 2011 ApJL 742, L21',exgal_sources='Mrk 231, Arp 220, Cloverleaf LOS')
N2 = Molecule('nitrogen','N2',2004,'N2','HD 124314 LOS','FUSE','UV',neutral=True,N=2,d_ref='Knauth et al. 2004 Nature 409, 636',lab_ref='Stark et al. 2000 ApJ 531, 321',notes=None)
CFp = Molecule('fluoromethylidynium cation','CF+',2006,'CF+','Orion Bar','IRAM, APEX','mm',cation=True,C=1,F=1,d_ref='Neufeld et al. 2006 A&A 454, L37',lab_ref='Plummer et al. 1986 JCP 84, 2427',notes=None,exgal=True,exgal_d_ref='Muller et al. 2016 A&A 589, L5',exgal_sources='PKS 1830-211 LOS')
PO = Molecule('phosphorous monoxide','PO',2007,'PO','VY Ca Maj','SMT','mm',neutral=True,O=1,P=1,d_ref='Tenenbaum et al. 2007 ApJ 666, L29',lab_ref='Bailleux et al. 2002 JMS 216, 465',notes=None)
O2 = Molecule('oxygen','O2',2007,'O2','Orion, rho Oph A','Odin, Herschel','mm, sub-mm',neutral=True,O=2,d_ref='Goldsmith et al. 2011 ApJ 737, 96',lab_ref='Endo & Mizushima 1982 Jpn J Appl Phys 21, L379; Drouin et al. 2010 J Quant Spec Rad Transf 111, 1167',notes='*Also Larsson et al. 2007 A&A 466, 999; Tentative in Goldsmith 2002 ApJ 576, 814')
AlO = Molecule('aluminum monoxide','AlO',2009,'AlO','VY Ca Maj','SMT','mm',neutral=True,O=1,Al=1,d_ref='Tenenbaum & Ziurys 2009 ApJ 693, L59',lab_ref='Yamada et al. 1990, JCP 92, 2146',notes=None)
CNm = Molecule('cyanide anion','CN-',2010,'CN-','IRC+10216','IRAM','mm',anion=True,C=1,N=1,d_ref='Agúndez et al. 2010 A&A 517, L2',lab_ref='Amano 2008 JCP 129, 244305',notes=None)
OHp = Molecule('hydroxyl cation','OH+',2010,'OH+','Sgr B2 LOS','APEX','sub-mm',cation=True,H=1,O=1,d_ref='Wyrowski et al. 2010 A&A 518, A26; Gerin et al. 2010 A&A 518, L110; Benz et al. 2010 A&A 521, L35',lab_ref='Bekooy et al. 1985 JCP 82, 3868',notes=None,exgal=True,exgal_d_ref='van der Werf et al. 2010 A&A 518, L42; Rangwala et al. 2011 ApJ 743, 94; Gonzalez-Alfonso et al. 2013 A&A 550, A25',exgal_sources='Mrk 231, Arp 220, NGC 4418')
SHp = Molecule('sulfanylium cation','SH+',2011,'SH+','Sgr B2','Herschel','sub-mm',cation=True,H=1,S=1,d_ref='Benz et al. 2010 A&A 521, L35',lab_ref='Brown et al. 2009 JMS 255, 68',notes='*Also in Menten et al. 2011 A&A 525, A77',exgal=True,exgal_d_ref='Muller et al. 2017 A&A 606, A109',exgal_sources='PKS 1830-211')
HClp = Molecule('hydrogen chloride cation','HCl+',2012,'HCl+','W31 LOS, W49 LOS','Herschel','sub-mm',cation=True,H=1,Cl=1,d_ref='de Luca et al. 2012 ApJ 751, L37',lab_ref='Gupta et al. 2012 ApJ 751, L38',notes=None)
SH = Molecule('mercapto radical','SH',2012,'SH','W49 LOS','SOFIA','sub-mm',radical=True,neutral=True,H=1,S=1,d_ref='Neufeld et al. 2012 A&A 542, L6',lab_ref='Morino & Kawaguchi 1995 JMS 170, 172; Klisch et al. 1996 ApJ 473, 1118',notes=None)
TiO = Molecule('titanium monoxide','TiO',2013,'TiO','VY Ca Maj','SMA','mm',neutral=True,O=1,Ti=1,d_ref='Kamiński et al. 2013 A&A 551, A113',lab_ref='Nakimi et al. 1998 JMS 191, 176',notes=None,exo=True,exo_d_ref='Haynes et al. 2015 ApJ 806, 146; Sedaghati et al. 2017 Nature 549, 238; Nugroho et al. 2017 ApJ 154, 221')
ArHp = Molecule('argonium','ArH+',2013,'ArH+','Crab Nebula','Herschel','sub-mm',cation=True,H=1,Ar=1,d_ref='Barlow et al. 2013 Science 342, 1343',lab_ref='Barlow et al. 2013 Science 342, 1343',notes=None,exgal=True,exgal_d_ref='Muller et al. 2015 A&A 582, L4',exgal_sources='PKS 1830-211 LOS')
NSp = Molecule('nitrogen sulfide cation','NS+',2018,'NS+','B1-b, TMC-1, L483','IRAM','mm',cation=True,N=1,S=1,d_ref='Cernicharo et al. 2018 ApJL 853, L22',lab_ref='Cernicharo et al. 2018 ApJL 853, L22',notes=None)

two_atom_list = [CH,CN,CHp,OH,CO,H2,SiO,CS,SO,SiS,NS,C2,NO,HCl,NaCl,AlCl,KCl,AlF,PN,SiC,CP,NH,SiN,SOp,COp,HF,N2,CFp,PO,O2,AlO,CNm,OHp,SHp,HClp,SH,TiO,ArHp,NSp]

#############################################################
#					Three Atom Molecules					#
#############################################################	

H2O = Molecule('water','H2O',1969,'H2O','Sgr B2, Orion, W49','Hat Creek','cm',neutral=True,H=2,O=1,d_ref='Cheung et al. 1969 Nature 221, 626',lab_ref='Golden et al. 1948 Phys Rev 73, 92',notes=None,ice=True,ice_d_ref='Gillett & Forrest 1973 ApJ 179, 483',ice_l_ref='Irvine & Pollack 1968 Icarus 8, 324',isos='HDO',isos_d_ref='[HDO] Turner et al. 1975 ApJ 198, L125',isos_l_ref='[HDO] de Lucia et al. 1974 J Phys Chem Ref Data 3, 211; Erlandsson & Cox 1956 J Chem Phys 25, 778',ppd=True,ppd_d_ref='Carr et al. 2004 ApJ 603, 213; Hogerheijde et al. 2011 Science 344, 338',exo=True,exo_d_ref='Tinetti et al. 2007 Nature 448, 169; Deming et al. 2014 ApJ 774, 95; Kreidberg et al. 2014 ApJL 793, L27; Kreidberg et al. 2015 ApJ 814, 66; Lockwood et al. 2014 ApJ 783, L29',exgal=True,exgal_d_ref='Churchwell et al. 1977 A&A 54, 969',exgal_sources='M33')
HCOp = Molecule('formylium cation','HCO+',1970,'HCO+','W3(OH), Orion, L134, Sgr A, W51','NRAO 36-ft','mm',cation=True,H=1,C=1,O=1,d_ref='Buhl & Snyder 1970 Nature 228, 267',lab_ref='Woods et al. 1975 PRL 35, 1269',notes=None,ppd=True,ppd_d_ref='Kastner et al. 1997 Science 277, 67; Dutrey et al. 1997 A&A 317, L55',ppd_isos='DCO+, H13CO+',ppd_isos_ref='[DCO+] van Dishoeck et al. 2003 A&A 400, L1 [H13CO+] van Zadelhoff et al. 2001 A&A 377, 566; van Dishoeck et al. 2003 A&A 400, L1',exgal=True,exgal_d_ref='Stark et al. 1979 ApJ 229, 118',exgal_sources='M82')
HCN = Molecule('hydrogen cyanide','HCN',1971,'HCN','W3(OH), Orion, Sgr A, W49, W51, DR 21','NRAO 36-ft','mm',neutral=True,H=1,C=1,N=1,d_ref='Snyder et al. 1971 ApJ 163, L47',lab_ref='de Lucia & Gordy 1969 Phys Rev 187, 58',notes=None,ppd=True,ppd_d_ref='Kastner et al. 1997 Science 277, 67; Dutrey et al. 1997 A&A 317, L55',ppd_isos='DCN, H13CN, HC15N',ppd_isos_ref='[DCN] Qi et al. 2008 ApJ 681, 1396 [H13CN] Guzman et al. 2015 ApJ 814, 53 [HC15N] Guzman et al. 2015 ApJ 814, 53',exgal=True,exgal_d_ref='Rickard et al. 1977 ApJ 214, 390',exgal_sources='NGC 253, M82',exo=True,exo_d_ref='Hawker et al. 2018 ApJL 863, L11')
OCS = Molecule('carbonyl sulfide','OCS',1971,'OCS','Sgr B2','NRAO 36-ft','mm',neutral=True,C=1,O=1,S=1,d_ref='Jefferts et al. 1971 ApJ 168, L111',lab_ref='King & Gordy 1954 Phys Rev 93, 407',notes=None,ice=True,ice_d_ref='Palumbo et al. 1995 ApJ 449, 674; Palumbo et al. 1997 ApJ 479, 839',ice_l_ref='Palumbo et al. 1995 ApJ 449, 674',exgal=True,exgal_d_ref='Mauersberger et al. 1995 A&A 294, 23',exgal_sources='NGC 253')
HNC = Molecule('hydrogen isocyanide','HNC',1972,'HNC','W51, NGC 2264','NRAO 36-ft','mm',neutral=True,H=1,C=1,N=1,d_ref='Snyder & Buhl 1972 Annals of the New York Academy of Science 194, 17; Zuckerman et al. 1972 ApJ 173, L125',lab_ref='Blackman et al. 1976 Nature 261, 395',notes=None,ppd=True,ppd_d_ref='Dutrey et al. 1997 A&A 317, L55',exgal=True,exgal_d_ref='Henkel et al. 1988 A&A 201, L23',exgal_sources='IC 342')
H2S = Molecule('hydrogen sulfide','H2S',1972,'H2S','W3, W3(OH), Orion, NGC 2264, Sgr B2, W51, DR 21(OH), NGC 7538','NRAO 36-ft','mm',neutral=True,H=2,S=1,d_ref='Thaddeus et al. 1972 ApJ 176, L73',lab_ref='Cupp et al. 1968 Phys Rev 171, 60',notes=None,exgal=True,exgal_d_ref='Hekkila et al. 1999 A&A 344, 817',exgal_sources='LMC',ppd=True,ppd_d_ref='Phuong et al. 2018 A&A 616, L5')
N2Hp = Molecule('protonated nitrogen','N2H+',1974,'N2H+','Sgr B2, DR 21, NGC 6334, NGC 2264','NRAO 36-ft','mm',cation=True,H=1,N=2,d_ref='Turner 1974 ApJ 193, L83; Green et al. 1974 ApJ 193, L89; Thaddues & Turner 1975 ApJ 201, L25',lab_ref='Saykally et al. 1976 ApJ 205, L101',notes=None,ppd=True,ppd_d_ref='Qi et al. 2003 ApJ 597, 986; Dutrey et al. 2007 A&A 464, 615',ppd_isos='N2D+',ppd_isos_ref='[N2D+] Huang et al. 2015 ApJL 809, L26',exgal=True,exgal_d_ref='Mauersberger & Henkel 1991 A&A 245, 457',exgal_sources='NGC 253, Maffei 2, IC 342, M82, NGC 6946')
C2H = Molecule('ethynyl radical','C2H',1974,'C2H','Orion','NRAO 36-ft','mm',radical=True,neutral=True,H=1,C=2,d_ref='Tucker et al. 1974 ApJ 193, L115',lab_ref='Sastry et al. 1981 ApJ 251, L119',notes=None,ppd=True,ppd_d_ref='Dutrey et al. 1997 A&A 317, L55',exgal=True,exgal_d_ref='Henkel et al. 1988 A&A 201, L23',exgal_sources='M82')
SO2 = Molecule('sulfur dioxide','SO2',1975,'SO2','Orion, Sgr B2','NRAO 36-ft','mm',neutral=True,O=2,S=1,d_ref='Snyder et al. 1975 ApJ 198, L81',lab_ref='Steenbeckeliers 1968 Ann. Soc. Sci. Brux 82, 331',notes=None,exgal=True,exgal_d_ref='Martin et al. 2003 A&A 411, L465',exgal_sources='NGC 253')
HCO = Molecule('formyl radical','HCO',1976,'HCO','W3, NGC 2024, W51, K3-50','NRAO 36-ft','mm',radical=True,neutral=True,H=1,C=1,O=1,d_ref='Snyder et al. 1976 ApJ 208, L91',lab_ref='Saito 1972 ApJ 178, L95',notes=None,exgal=True,exgal_d_ref='Sage & Ziurys 1995 ApJ 447, 625; Garcia-Burillo et al. 2002 ApJ 575, L55',exgal_sources='M82')
HNO = Molecule('nitroxyl radical','HNO',1977,'HNO','Sgr B2, NGC 2024','NRAO 36-ft','mm',neutral=True,H=1,O=1,N=1,d_ref='Ulich et al. 1977 ApJ 217, L105',lab_ref='Saito & Takagi 1973 JMS 47, 99',notes=None)
HCSp = Molecule('protonated carbon monosulfide','HCS+',1981,'HCS+','Orion, Sgr B2','NRAO 36-ft, Bell 7-m','mm',cation=True,H=1,C=1,S=1,d_ref='Thaddeus et al. 1981 ApJ 246, L41',lab_ref='Gudeman et al. 1981 ApJ 246, L47',notes=None,exgal=True,exgal_d_ref='Muller et al. 2013 A&A 551, A109',exgal_sources='PKS 1830-211 LOS')
HOCp = Molecule('hydroxymethyliumylidene','HOC+',1983,'HOC+','Sgr B2','FCRAO 14-m, Onsala','mm',cation=True,H=1,C=1,O=1,d_ref='Woods et al. 1983 ApJ 270, 583',lab_ref='Gudeman et al. 1982 PRL 48, 1344',notes='*Confirmed in 1995 ApJ 455, L73',exgal=True,exgal_d_ref='Usero et al. 2004 A&A 419, 897',exgal_sources='NGC 1068')
SiC2 = Molecule('silacyclopropynylidene','SiC2',1984,'SiC2','IRC+10216','NRAO 36-ft, Bell 7-m','mm',cyclic=True,radical=True,neutral=True,C=2,Si=1,d_ref='Thaddeus et al. 1984 ApJ 283, L45',lab_ref='Michalopoulos et al. 1984 JCP 80, 3556',notes=None)
C2S = Molecule('dicarbon sulfide','C2S',1987,'C2S','TMC-1, IRC+10216, Sgr B2','Nobeyama, IRAM','cm, mm',radical=True,neutral=True,C=2,S=1,d_ref='Saito et al. 1987 ApJ 317, L115',lab_ref='Saito et al. 1987 ApJ 317, L115',notes='*Also Cernicharo et al. 1987 A&A 181, L9',exgal=True,exgal_d_ref='Martin et al. 2006 ApJS 164, 450',exgal_sources='NGC 253')
C3 = Molecule('tricarbon','C3',1988,'C3','IRC+10216','KPNO 4-m','IR',neutral=True,C=3,d_ref='Hinkle et al. 1988 Science 241, 1319',lab_ref='Gausset et al. 1965 ApJ 142, 45',exgal=True,exgal_d_ref='Welty et al. 2012 MNRAS 428, 1107', exgal_sources='SMC',notes=None)
CO2 = Molecule('carbon dioxide','CO2',1989,'CO2','AFGL 961 LOS, AFGL 989 LOS, AFGL 890 LOS','IRAS','IR',neutral=True,C=1,O=2,d_ref='d\'Hendecourt & Jourdain de Muizon 1989 A&A 223, L5; van Dishoeck et al. 1996 A&A 315, L349',lab_ref='d\'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453; Paso et al. 1980 JMS 79, 236; Reichle & Young 1972 Can J Phys 50, 2662',notes='*First detected in ices, then in gas phase',ice=True,ice_d_ref='d\'Hendecourt & Jourdain de Muizon 1989 A&A 223, L5',ice_l_ref='d\'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453',ppd=True,ppd_d_ref='Carr & Najita 2008 Science 319, 1504',exo=True,exo_d_ref='Stevenson et al. 2010 Nature 464, 1161; Madhusudhan et al. 2011 Nature 469, 64; Lanotte et al. 2014 A&A 572, A73')
CH2 = Molecule('methylene','CH2',1989,'CH2','Orion','NRAO 12-m','mm',radical=True,neutral=True,H=2,C=1,d_ref='Hollis et al. 1989 ApJ 346, 794',lab_ref='Lovas et al. 1983 ApJ 267, L131',notes='*Confirmed in 1995 ApJ 438, 259')
C2O = Molecule('dicarbon monoxide','C2O',1991,'C2O','TMC-1','Nobeyama','cm',radical=True,neutral=True,C=2,O=1,d_ref='Ohishi et al. 1991 ApJ 380, L39',lab_ref='Yamada et al. 1985 ApJ 290, L65',notes=None)
MgNC = Molecule('magnesium isocyanide','MgNC',1993,'MgNC','IRC+10216','IRAM','mm',radical=True,neutral=True,C=1,N=1,Mg=1,d_ref='Guélin et al. 1986 A&A 157, L17',lab_ref='Kawaguchi et al. 1993 ApJ 406, L39',notes='*Actually identified in Kawaguchi et al. 1993 ApJ 406, L39 and Guélin et al. 1993 A&A 280, L19')
NH2 = Molecule('amidogen','NH2',1993,'NH2','Sgr B2 LOS','CSO','sub-mm',radical=True,neutral=True,H=2,N=1,d_ref='van Dishoeck et al. 1993 ApJ 416, L83',lab_ref='Charo et al. 1981 ApJ 244, L111',notes=None,exgal=True,exgal_d_ref='Muller et al. 2014 A&A 566, A112',exgal_sources='PKS 1830-211 LOS')
NaCN = Molecule('sodium cyanide','NaCN',1994,'NaCN','IRC+10216','NRAO 12-m','mm',neutral=True,C=1,N=1,Na=1,d_ref='Turner et al. 1994 ApJ 426, L97',lab_ref='van Vaals et al. 1984 Chem Phys 86, 147',notes=None)
N2O = Molecule('nitrous oxide','N2O',1994,'N2O','Sgr B2','NRAO 12-m','mm',neutral=True,O=1,N=2,d_ref='Ziurys et al. 1994 ApJ 436, L181',lab_ref='Lovas 1978 J Phys Chem Ref Data 7, 1445',notes=None)
MgCN = Molecule('magnesium cyanide','MgCN',1995,'MgCN','IRC+10216','NRAO 12-m, IRAM','mm',radical=True,neutral=True,C=1,N=1,Mg=1,d_ref='Ziurys et al. 1995 ApJ 445, L47',lab_ref='Anderson et al. 1994 ApJ 429, L41',notes=None)
H3p = Molecule('','H3+',1996,'H3+','GL2136 LOS, W33 LOS','UKIRT','IR',cation=True,H=3,d_ref='Geballe & Oka 1996 Nature 384, 334',lab_ref='Oka 1980 PRL 45, 531',notes=None,exgal=True,exgal_d_ref='Geballe et al. 2006 ApJ 644, 907',exgal_sources='IRAS 08572+3915')
SiCN = Molecule('silicon monocyanide radical','SiCN',2000,'SiCN','IRC+10216','IRAM','mm',radical=True,neutral=True,C=1,N=1,Si=1,d_ref='Guélin et al. 2000 A&A 363, L9',lab_ref='Apponi et al. 2000 ApJ 536, L55',notes=None)
AlNC = Molecule('aluminum isocyanide','AlNC',2002,'AlNC','IRC+10216','IRAM','mm',neutral=True,C=1,N=1,Al=1,d_ref='Ziurys et al. 2002 ApJ 564, L45',lab_ref='Robinson et al. 1997 Chem Phys Lett 278, 1',notes=None)
SiNC = Molecule('silicon monoisocyanide','SiNC',2004,'SiNC','IRC+10216','IRAM','mm',radical=True,neutral=True,C=1,N=1,Si=1,d_ref='Guélin et al. 2004 A&A 426, L49',lab_ref='Apponi et al. 2000 ApJ 536, L55',notes=None)
HCP = Molecule('phosphaethyne','HCP',2007,'HCP','IRC+10216','IRAM','mm',neutral=True,H=1,C=1,P=1,d_ref='Agúndez et al. 2007 ApJ 662, L91',lab_ref='Bizzocchi et al. 2001 JMS 205, 110',notes='*First attempt 1990 ApJ 365, 59. Confirmed 2008 ApJ 684, 618')
CCP = Molecule('dicarbon phosphide radical','CCP',2008,'CCP','IRC+10216','ARO 12-m','mm',radical=True,neutral=True,C=2,P=1,d_ref='Halfen et al. 2008 ApJ 677, L101',lab_ref='Halfen et al. 2008 ApJ 677, L101',notes=None)
AlOH = Molecule('aluminum hydroxide','AlOH',2010,'AlOH','VY Ca Maj','ARO 12-m, SMT','mm',neutral=True,H=1,O=1,Al=1,d_ref='Tenenbaum & Ziurys 2010 ApJ 712, L93',lab_ref='Apponi et al. 1993 ApJ 414, L129',notes=None)
H2Op = Molecule('oxidaniumyl','H2O+',2010,'H2O+','Sgr B2, Sgr B2 LOS, NGC 6334, DR 21','Herschel','sub-mm',cation=True,H=2,O=1,d_ref='Ossenkopf et al. 2010 A&A 518, L111; Gerin et al. 2010 A&A 518, L110',lab_ref='Strahan et al. 1986 JCP 85, 1252; Murtz et al. 1998 JCP 109, 9744 ',notes=None,exgal=True,exgal_d_ref='Weiss et al. 2010 A&A 521, L1',exgal_sources='M82')
H2Clp = Molecule('chloronium','H2Cl+',2010,'H2Cl+','Sgr B2, Sgr B2 LOS, NGC 6334, NGC 6334 LOS','Herschel','sub-mm',cation=True,H=2,Cl=1,d_ref='Lis et al. 2010 A&A 521, L9',lab_ref='Araki et al. 2001 JMS 210, 132',notes=None,exgal=True,exgal_d_ref='Muller et al. 2014 A&A 566, L6',exgal_sources='PKS 1830-211 LOS')
KCN = Molecule('potassium cyanide','KCN',2010,'KCN','IRC+10216','ARO 12-m, SMT, IRAM','mm',neutral=True,C=1,N=1,K=1,d_ref='Pulliam et al. 2010 ApJ 727, L181',lab_ref='Torring et al. 1980 JCP 73, 4875',notes=None)
FeCN = Molecule('iron cyanide','FeCN',2011,'FeCN','IRC+10216','ARO 12-m','mm',neutral=True,C=1,N=1,Fe=1,d_ref='Zack et al. 2011 ApJ 733, L36',lab_ref='Flory & Ziurys 2011 JCP 135, 184303',notes=None)
HO2 = Molecule('hydroperoxyl radical','HO2',2012,'HO2','rho Oph A','IRAM, APEX','mm',radical=True,neutral=True,H=1,O=2,d_ref='Parise et al. 2012 A&A 541, L11',lab_ref='Beers & Howard 1975 JCP 63, 4212; Saito 1977 JMS 65, 229; Charo & de Lucia 1982 JMS 94, 426',notes=None)
TiO2 = Molecule('titanium dioxide','TiO2',2013,'TiO2','VY Ca Maj','SMA, PdBI','mm',neutral=True,O=1,Ti=1,d_ref='Kamiński et al. 2013 A&A 551, A113',lab_ref='Brunken 2008 APJ 676, 1367; Kania et al. 2011 JMS 268, 173',notes=None)
CCN = Molecule('cyanomethylidyne','CCN',2014,'CCN','IRC+10216','ARO 12-m, SMT','mm',radical=True,neutral=True,C=2,N=1,d_ref='Anderson & Ziurys 2014 ApJ 795, L1',lab_ref='Anderson et al. 2015 JMS 307, 1',notes=None)
SiCSi = Molecule('disilicon carbide','SiCSi',2015,'SiCSi','IRC+10216','IRAM','mm',neutral=True,C=1,Si=2,d_ref='Cernicharo et al. 2015 ApJ 806, L3',lab_ref='McCarthy 2015 JPC Lett 6, 2107',notes=None)
S2H = Molecule('hydrogen disulfide','S2H',2017,'S2H','Horsehead PDR','IRAM','mm',neutral=True,H=1,S=2,d_ref='Fuente et al. 2017 ApJ 851, L49',lab_ref='Tanimoto et al. 2000 JMS 199, 73',notes=None)
HCS = Molecule('thioformyl','HCS',2018,'HCS','L483','IRAM','mm',radical=True,neutral=True,H=1,C=1,S=1,d_ref='Agúndez et al. 2018 A&A 611, L1',lab_ref='Habara et al. 2002 JCP 116, 9232',notes=None)
HSC = Molecule('sulfhydryl carbide','HSC',2018,'HSC','L483','IRAM','mm',radical=True,neutral=True,H=1,C=1,S=1,d_ref='Agúndez et al. 2018 A&A 611, L1',lab_ref='Habara 2000 JCP 112, 10905',notes=None)
NCO = Molecule('isocyanate radical','NCO',2018,'NCO','L483','IRAM','mm',radical=True,neutral=True,C=1,O=1,N=1,d_ref='Marcelino et al. 2018 A&A 612, L10',lab_ref='Kawaguchi et al. 1985 Mol Phys 55, 341; Saito and Amano 1970 JMS34, 383',notes=None)

three_atom_list = [H2O,HCOp,HCN,OCS,HNC,H2S,N2Hp,C2H,SO2,HCO,HNO,HCSp,HOCp,SiC2,C2S,C3,CO2,CH2,C2O,MgNC,NH2,NaCN,N2O,MgCN,H3p,SiCN,AlNC,SiNC,HCP,CCP,AlOH,H2Op,H2Clp,KCN,FeCN,HO2,TiO2,CCN,SiCSi,S2H,HCS,HSC,NCO]

#############################################################
#					Four Atom Molecules						#
#############################################################

NH3 = Molecule('ammonia','NH3',1968,'NH3','Galactic Center','Hat Creek','cm',neutral=True,H=3,N=1,d_ref='Cheung et al. 1968 PRL 25, 1701',lab_ref='Cleeton & Williams 1934 Phys Rev 45, 234',notes=None,ice=True,ice_d_ref='Lacy et al. 1998 ApJ 501, L105',ice_l_ref='d\'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453',ppd=True,ppd_d_ref='Salinas et al. 2016 A&A 591, A122',exgal=True,exgal_d_ref='Martin & Ho 1979 A&A 74, L7',exgal_sources='IC 342, NGC 253')
H2CO = Molecule('formaldehyde','H2CO',1969,'H2CO','M17 LOS, M3 LOS, W49 LOS, NGC 2024 LOS, DR 21 LOS, W43 LOS, W44 LOS, W51 LOS, Sgr A LOS, Sgr B2 LOS, W33 LOS, NGC 6334 LOS, Cas A LOS','NRAO 140-ft','cm',neutral=True,H=2,C=1,O=1,d_ref='Snyder et al. 1969 PRL 22, 679',lab_ref='Shinegari 1967 J Phys Soc Jpn 23, 404',notes=None,ice=True,ice_d_ref='Keane et al. 2001 A&A 376, 254',ice_l_ref='Schutte et al. 1993 Icarus 104, 118',ppd=True,ppd_d_ref='Dutrey et al. 1997 A&A 317, L55',exgal=True,exgal_d_ref='Gardner & Whiteoak 1974 Nature 247, 526',exgal_sources='NGC 253, NGC 4945')
HNCO = Molecule('isocyanic acid','HNCO',1972,'HNCO','Sgr B2','NRAO 36-ft','cm, mm',neutral=True,H=1,C=1,O=1,N=1,d_ref='Snyder & Buhl 1972 ApJ 177, 619',lab_ref='Kewley et al. 1963 JMS 10, 418',notes=None,exgal=True,exgal_d_ref='Nguyen-Q-Rieu et al. 1991 A&A 241, L33',exgal_sources='NGC 253, Maffei 2, IC 342')
H2CS = Molecule('thioformaldehyde','H2CS',1973,'H2CS','Sgr B2 LOS','Parkes 64-m','cm',neutral=True,H=2,C=1,S=1,d_ref='Sinclair et al. 1973 Aust. J. Phys. 26, 85',lab_ref='Johnson & Powell 1970 Science 169, 679',notes=None,exgal=True,exgal_d_ref='Martin et al. 2006 ApJS 164, 450',exgal_sources='NGC 253')
C2H2 = Molecule('acetylene','C2H2',1976,'C2H2','IRC+10216','Mayall 4-m','IR',neutral=True,H=2,C=2,d_ref='Ridgway et al. 1976 Nature 264, 345',lab_ref='Baldacci et al. 1973 JMS 48, 600',notes=None,ppd=True,ppd_d_ref='Lahuis et al. 2006 ApJ 66, L145',exgal=True,exgal_d_ref='Matsuura et al. 2002 ApJ 580, L133',exgal_sources='LMC')
C3N = Molecule('cyanoethynyl radical','C3N',1977,'C3N','IRC+10216, TMC-1','NRAO 36-ft, Onsala','cm, mm',radical=True,neutral=True,C=3,N=1,d_ref='Guelin & Thaddeus 1977 ApJ 212, L81',lab_ref='Gottlieb et al. 1983 ApJ 275, 916',notes='*Confirmed in Friberg et al. 1980 ApJ 241, L99')
HNCS = Molecule('isothiocyanic acid','HNCS',1979,'HNCS','Sgr B2','Bell 7-m, NRAO 36-ft','mm',neutral=True,H=1,C=1,N=1,S=1,d_ref='Frerking et al. 1979 ApJ 234, L143',lab_ref='Kewley et al. 1963 JMS 10, 418',notes=None)
HOCOp = Molecule('protonated carbon dioxide','HOCO+',1981,'HOCO+','Sgr B2','Bell 7-m','mm',cation=True,H=1,C=1,O=2,d_ref='Thaddeus et al. 1981 ApJ 246, L41',lab_ref='Green et al. 1976 Chem Phys 17, 479; Bogey et al. 1984 A&A 138, L11',notes=None,exgal=True,exgal_d_ref='Aladro et al. 2015 A&A 579, A101; Martin et al. 2006 ApJS 164, 450',exgal_sources='NGC 253')
C3O = Molecule('tricarbon monoxide','C3O',1985,'C3O','TMC-1','NRAO 12-m, FCRAO 14-m, Nobeyama','cm, mm',radical=True,neutral=True,C=3,O=1,d_ref='Matthews et al. 1984 Nature 310, 125',lab_ref='Brown et al. 1983 JACS 105, 6496',notes='*Confirmed in Brown et al. 1985 ApJ 297, 302 and Kaifu et al. 2004 PASJ 56, 69')
lC3H = Molecule('propynylidyne radical','l-C3H',1985,'l-C3H','TMC-1, IRC+10216','NRAO 36-ft, Onsala, Bell 7-m, U Mass 14-m','cm, mm',radical=True,neutral=True,H=1,C=3,d_ref='Thaddeus et al. 1985 ApJ 294, L49',lab_ref='Gottlieb et al. 1985 ApJ 294, L55',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS')
HCNHp = Molecule('protonated hydrogen cyanide','HCNH+',1986,'HCNH+','Sgr B2','NRAO 12-m, MWO 4.9-m','mm',cation=True,H=2,C=1,N=1,d_ref='Ziurys & Turner 1986 ApJ 302, L31',lab_ref='Bogey et al. 1985 JCP 83, 3703; Altman et al. 1984 JCP 80, 3911',notes=None)
H3Op = Molecule('hydronium','H3O+',1986,'H3O+','Orion, Sgr B2','NRAO 12-m, MWO 4.9-m','mm',cation=True,H=3,O=1,d_ref='Wootten et al. 1986 A&A 166, L15; Hollis et al. 1986 Nature 322, 524',lab_ref='Plummer et al. 1985 JCP 83, 1428; Bogey et al. 1985 A&A 148, L11; Liu & Oka 1985 PRL 54, 1787',notes='*Confirmed in Wootten et al. 1991 ApJ 390, L79',exgal=True,exgal_d_ref='van der Tak et al. 2007 A&A 477, L5',exgal_sources='M82, Arp 220')
C3S = Molecule('tricarbon monosulfide','C3S',1987,'C3S','TMC-1, IRC+10216','Nobeyama, IRAM','cm, mm',radical=True,neutral=True,C=3,S=1,d_ref='Yamamoto et al. 1987 ApJ 317, L119',lab_ref='Yamamoto et al. 1987 ApJ 317, L119',notes=None)
cC3H = Molecule('cyclopropenylidene radical','c-C3H',1987,'c-C3H','TMC-1','Nobeyama','mm',cyclic=True,radical=True,neutral=True,H=1,C=3,d_ref='Yamamoto et al. 1987 ApJ 322, L55',lab_ref='Yamamoto et al. 1987 ApJ 322, L55',notes=None,exgal='Tentative',exgal_d_ref='Martin et al. 2006 ApJS 164, 450',exgal_sources='NGC 253')
HC2N = Molecule('cyanocarbene radical','HC2N',1991,'HC2N','IRC+10216','IRAM','mm',radical=True,neutral=True,H=1,C=2,N=1,d_ref='Guélin & Cernicharo 1991 A&A 244, L21',lab_ref='Saito et al. 1984 JCP 80, 1427; Brown et al. 1990 JMS 143, 203',notes=None)
H2CN = Molecule('methylene amidogen radical','H2CN',1994,'H2CN','TMC-1','NRAO 12-m','cm',radical=True,neutral=True,H=2,C=1,N=1,d_ref='Ohishi et al. 1994 ApJ 427, L51',lab_ref='Yamamoto & Saito 1992 JCP 96, 4157',notes=None)
SiC3 = Molecule('silicon tricarbide','SiC3',1999,'SiC3','IRC+10216','NRAO 12-m','mm',neutral=True,cyclic=True,C=3,Si=1,d_ref='Apponi et al. 1999 ApJ 516, L103',lab_ref='Apponi et al. 1999 JCP 111, 3911; McCarthy et al. JCP 110, 1064',notes=None)
CH3 = Molecule('methyl radical','CH3',2000,'CH3','Sgr A LOS','ISO','IR',radical=True,neutral=True,H=3,C=1,d_ref='Feuchtgruber et al. 2000 ApJ 535, L111',lab_ref='Yamada et al. 1981 JCP 75, 5256',notes=None)
C3Nm = Molecule('cyanoethynyl anion','C3N-',2008,'C3N-','IRC+10216','IRAM','mm',anion=True,C=3,N=1,d_ref='Thaddeus et al. 2008 ApJ 677, 1132',lab_ref='Thaddeus et al. 2008 ApJ 677, 1132',notes=None)
PH3 = Molecule('phosphine','PH3',2008,'PH3','IRC+10216, CRL 2688','IRAM, Herschel, SMT, CSO','mm, sub-mm',neutral=True,H=3,P=1,d_ref='Agúndez et al. 2008 A&A 485, L33',lab_ref='Cazzoli & Puzzarini 2006 JMS 239, 64; Sousa-Silva et al. 2013 JMS 288, 28; Muller 2013 JQSRT 130, 335',notes='*Confirmed in Agúndez et al. 2014 ApJL 790, L27')
HCNO = Molecule('fulminic acid','HCNO',2009,'HCNO','B1-b, L1544, L183, L1527','IRAM','mm',neutral=True,H=1,C=1,O=1,N=1,d_ref='Marcelino et al. 2009 ApJ 690, L27',lab_ref='Winnewisser & Winnewisser 1971 Z Naturforsch 26, 128',notes=None)
HOCN = Molecule('cyanic acid','HOCN',2009,'HOCN','Sgr B2','Bell 7-m, NRAO 36-ft, NRAO 12-m','mm',neutral=True,H=1,C=1,O=1,N=1,d_ref='Brünken et al. 2009 ApJ 697, 880',lab_ref='Brünken et al. 2009 ApJ 697, 880',notes='*Confirmed in Brünken et al. 2010 A&A 516, A109')
HSCN = Molecule('thiocyanic acid','HSCN',2009,'HSCN','Sgr B2','ARO 12-m','mm',neutral=True,H=1,C=1,N=1,S=1,d_ref='Halfen et al. 2009 ApJ 702, L124',lab_ref='Brunken et al. 2009 ApJ 706, 1588',notes=None)
HOOH = Molecule('hydrogen peroxide','HOOH',2011,'HOOH','rho Oph A','APEX','mm',neutral=True,H=2,O=2,d_ref='Bergman et al. 2011 A&A 531, L8',lab_ref='Petkie et al. 1995 JMS 171, 145; Helminger et al. 1981 JMS 85, 120',notes=None)
lC3Hp = Molecule('cyclopropynylidynium cation','l-C3H+',2012,'l-C3H+','Horsehead PDR','IRAM','mm',cation=True,H=1,C=3,d_ref='Pety et al. 2012 A&A 549, A68',lab_ref='Brunken et al. 2014 ApJ 783, L4',notes=None)
HMgNC = Molecule('hydromagnesium isocyanide','HMgNC',2013,'HMgNC','IRC+10216','IRAM','mm',neutral=True,H=1,C=1,N=1,Mg=1,d_ref='Cabezas et al. 2013 ApJ 75, 133',lab_ref='Cabezas et al. 2013 ApJ 75, 133',notes=None)
HCCO = Molecule('ketenyl radical','HCCO',2015,'HCCO','Lupus-1A, L483','IRAM','mm',radical=True,neutral=True,H=1,C=2,O=1,d_ref='Agúndez et al. 2015 A&A 577, L5',lab_ref='Endo & Hirota 1987 JCP 86, 4319; Oshima & Endo 1993 JMS 159, 458',notes=None)
CNCN = Molecule('isocyanogen','CNCN',2018,'CNCN','L483','IRAM','mm',neutral=True,C=2,N=2,d_ref='Agundez et al. 2018 ApJL 861, L22',lab_ref='Gerry et al. 1990 JMS 140, 147; Winnewisser et al. 1992 JMS 153, 635',notes=None)

four_atom_list = [NH3,H2CO,HNCO,H2CS,C2H2,C3N,HNCS,HOCOp,C3O,lC3H,HCNHp,H3Op,C3S,cC3H,HC2N,H2CN,SiC3,CH3,C3Nm,PH3,HCNO,HOCN,HSCN,HOOH,lC3Hp,HMgNC,HCCO,CNCN]

#############################################################
#					Five Atom Molecules						#
#############################################################

HC3N = Molecule('cyanoacetylene','HC3N',1971,'HC3N','Sgr B2','NRAO 140-ft','cm',neutral=True,H=1,C=3,N=1,d_ref='Turner 1971 ApJ 163, L35',lab_ref='Tyler & Sheridan 1963 Trans Faraday Soc 59, 2661',notes='*Confirmed in Dickinson 1972 AL 12, 235',ppd=True,ppd_d_ref='Chapillon et al. 2012 ApJ 756, 58',exgal=True,exgal_d_ref='Mauersberger et al. 1990 A&A 236, 63; Henkel et al. 1988 A&A 201, L23',exgal_sources='NGC 253')
HCOOH = Molecule('formic acid','HCOOH',1971,'HCOOH','Sgr B2','NRAO 140-ft','cm',neutral=True,H=2,C=1,O=2,d_ref='Zukerman et al. 1971 ApJ 163, L41',lab_ref='Zukerman et al. 1971 ApJ 163, L41; Bellet et al. 1971 J Mol Struct 9, 49; Bellet et al. 1971 J Mol Struct 9, 65',notes='*Confirmed in Winnewisser & Churchwell 1975 ApJ 200, L33',ice=True,ice_d_ref='Schutte et al. 1999 A&A 343, 966',ice_l_ref='Schutte et al. 1999 A&A 343, 966',ppd=True,ppd_d_ref='Favre et al. 2018 ApJL 862, L2')
CH2NH = Molecule('methanimine','CH2NH',1973,'CH2NH','Sgr B2','Parkes 64-m','cm',neutral=True,H=3,C=1,N=1,d_ref='Godfrey et al. 1973 ApL 13, 119',lab_ref='Godfrey et al. 1973 ApL 13, 119; Johnson & Lovas 1972 CPL 15, 65',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS')
NH2CN = Molecule('cyanamide','NH2CN',1975,'NH2CN','Sgr B2','NRAO 36-ft','mm',neutral=True,H=2,C=1,N=2,d_ref='Turner et al. 1975 ApJ 201, L149',lab_ref='Tyler et al. 1972 JMS 43, 248; Miller et al. 1962 JMS 8, 153; Lide 1962 JMS 8, 142; Johnson & Suenram 1976 ApJ 208, 245',notes=None,exgal=True,exgal_d_ref='Martin et al. 2006 ApJS 164, 450',exgal_sources='NGC 253')
H2CCO = Molecule('ketene','H2CCO',1977,'H2CCO','Sgr B2','NRAO 36-ft','mm',neutral=True,H=2,C=2,O=1,d_ref='Turner 1977 ApJ 213, L75',lab_ref='Johnson & Strandberg 1952 JCP 20, 687; Johns et al. 1972 JMS 42, 523',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS')
C4H = Molecule('butadiynyl radical','C4H',1978,'C4H','IRC+10216','NRAO 36-ft','mm',radical=True,neutral=True,H=1,C=4,d_ref='Guélin et al. 1978 ApJ 224, L27',lab_ref='Gottlieb et al. 1983 ApJ 275, 916',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS')
SiH4 = Molecule('silane','SiH4',1984,'SiH4','IRC+10216','IRTF','IR',neutral=True,H=4,Si=1,d_ref='Goldhaber and Betz 1977 ApJ 279, L55',lab_ref='Goldhaber and Betz 1977 ApJ 279, L55',notes=None)
cC3H2 = Molecule('cyclopropenylidene','c-C3H2',1985,'c-C3H2','Sgr B2, Orion, TMC-1','Bell 7-m','cm, mm',neutral=True,cyclic=True,H=2,C=3,d_ref='Thaddeus et al. 1985 ApJ 299, L63',lab_ref='Thaddeus et al. 1985 ApJ 299, L63',notes='*See also Vrtilek et al. 1987 ApJ 314, 716',ppd=True,ppd_d_ref='Qi et al. 2013 ApJL 765, L14',exgal=True,exgal_d_ref='Seaquist & Bell 1986 ApJ 303, L67',exgal_sources='NGC 5128')
CH2CN = Molecule('cyanomethyl radical','CH2CN',1988,'CH2CN','TMC-1, Sgr B2','FCRAO 14-m, NRAO 140-ft, Onsala, Nobeyama','cm',radical=True,neutral=True,H=2,C=2,N=1,d_ref='Irvine et al. 1988 ApJ 334, L107',lab_ref='Saito et al. 1988 ApJ 334, L113',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS')
C5 = Molecule('pentacarbon','C5',1989,'C5','IRC+10216','KPNO 4-m','IR',neutral=True,C=5,d_ref='Bernath et al. 1989 Science 244, 562',lab_ref='Vala et al. 1989 JCP 90, 595',notes=None)
SiC4 = Molecule('silicon tetracarbide','SiC4',1989,'SiC4','IRC+10216','Nobeyama','cm, mm',neutral=True,C=4,Si=1,d_ref='Ohishi et al. 1989 ApJ 345, L83',lab_ref='Ohishi et al. 1989 ApJ 345, L83',notes=None)
H2CCC = Molecule('propadienylidene','H2CCC',1991,'H2CCC','TMC-1','IRAM, Effelsberg 100-m','cm, mm',neutral=True,H=2,C=3,d_ref='Cernicharo et al. 1991 ApJ 368, L39',lab_ref='Vrtilek et al. 1990 ApJ 364, L53',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS')
CH4 = Molecule('methane','CH4',1991,'CH4','NGC 7538 LOS','IRTF','IR',neutral=True,H=4,C=1,d_ref='Lacy et al. 1991 ApJ 376, 556',lab_ref='Champion et al. 1989 JMS 133, 256; d\'Hendecourt & Allamandola 1986 A&A Supp Ser. 64, 453 ',notes=None,ice=True,ice_d_ref='Lacy et al. 1991 ApJ 376, 556',ice_l_ref='d\'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453',ppd=True,ppd_d_ref='Gibb et al. 2013 ApJL 776, L28',exo=True,exo_d_ref='Swain et al. 2008 Nature 452, 329; Barman et al. 2011 ApJ 733, 65; Stevenson et al. 2014 ApJ 791, 36; Barman et al. 2015 ApJ 804, 61')
HCCNC = Molecule('isocyanoacetylene','HCCNC',1992,'HCCNC','TMC-1','Nobeyama','cm, mm',neutral=True,H=1,C=3,N=1,d_ref='Kawaguchi et al. 1992 ApJ 386, L51',lab_ref='Kruger et al. 2010 Ang. Chem. 23, 1644',notes=None)
HNCCC = Molecule('','HNCCC',1992,'HNCCC','TMC-1','Nobeyama','cm',neutral=True,H=1,C=3,N=1,d_ref='Kawaguchi et al. 1992 ApJ 396, L49',lab_ref='Kawaguchi et al. 1992 ApJ 396, L49',notes=None)
H2COHp = Molecule('protonated formaldehyde','H2COH+',1996,'H2COH+','Sgr B2, Orion, W51','Nobeyama, NRAO 12-m','cm, mm',cation=True,H=3,C=1,O=1,d_ref='Ohishi et al. 1996 ApJ 471, L61',lab_ref='Chomiak et al. 1994 Can J Phys 72, 1078',notes=None)
C4Hm = Molecule('butadiynyl anion','C4H-',2007,'C4H-','IRC+10216','IRAM','mm',anion=True,H=1,C=4,d_ref='Cernicharo et al. 2007 A&A 467, L37',lab_ref='Gupta et al. 2007 ApJ 655, L57',notes=None)
CNCHO = Molecule('cyanoformaldehyde','CNCHO',2007,'CNCHO','Sgr B2','GBT','cm',neutral=True,H=1,C=2,O=1,N=1,d_ref='Remijan et al. 2008 ApJ 675, L85',lab_ref='Bogey et al. 1988 CPL 146, 227; Bogey et al. 1995 JMS 172, 344',notes=None)
HNCNH = Molecule('carbodiimide','HNCNH',2012,'HNCNH','Sgr B2','GBT','cm',neutral=True,H=2,C=1,N=2,d_ref='McGuire et al. 2012 ApJ 758, L33',lab_ref='Birk et al. 1989 JMS 135, 402; Wagener et al. 1995 JMS 170, 323; Jabs et al. 1997 Chem Phys 225, 77',notes=None)
CH3O = Molecule('methoxy radical','CH3O',2012,'CH3O','B1-b','IRAM','mm',radical=True,neutral=True,H=3,C=1,O=1,d_ref='Cernicharo et al. 2012 ApJ 759, L43',lab_ref='Momose et al. 1988 JCP 88, 5338; Endo et al. 1984 JCP 81, 122',notes=None)
NH3Dp = Molecule('ammonium ion','NH3D+',2013,'NH3D+','Orion, B1-b','IRAM','mm',cation=True,H=4,N=1,d_ref='Gupta et al. 2013 ApJ 778, L1',lab_ref='Gupta et al. 2013 ApJ 778, L1',notes='*Confirmed in Marcelino et al. 2018 A&A 612, L10')
H2NCOp = Molecule('protonated isocyanic acid','H2NCO+',2013,'H2NCO+','Sgr B2, L483','GBT','cm',cation=True,H=2,C=1,O=1,N=1,d_ref='Cernicharo et al. 2013 ApJ 771, L10',lab_ref='Cernicharo et al. 2013 ApJ 771, L10',notes='*See also Doménech et al. 2013 ApJ 77, L11')
NCCNHp = Molecule('protonated cyanogen','NCCNH+',2015,'NCCNH+','TMC-1, L483','Yebes 40-m, IRAM','cm, mm',cation=True,H=1,C=2,N=2,d_ref='Agúndez et al. 2015 A&A 579, L10',lab_ref='Amano & Scappini 1991 JCP 95, 2280; Gottlieb et al. 200 JCP 113, 1910',notes=None)
CH3Cl = Molecule('chloromethane','CH3Cl',2017,'CH3Cl','IRAS 16293','ALMA','mm',neutral=True,H=3,C=1,Cl=1,d_ref='Fayolle et al. 2017 Nature Astron. 1, 702',lab_ref='Wlodarczak et al. 1986 JMS 116, 251',notes=None)

five_atom_list = [HC3N,HCOOH,CH2NH,NH2CN,H2CCO,C4H,SiH4,cC3H2,CH2CN,C5,SiC4,H2CCC,CH4,HCCNC,HNCCC,H2COHp,C4Hm,CNCHO,HNCNH,CH3O,NH3Dp,H2NCOp,NCCNHp,CH3Cl]


#############################################################
#					Six Atom Molecules						#
#############################################################


CH3OH = Molecule('methanol','CH3OH',1970,'CH3OH','Sgr A, Sgr B2','NRAO 140-ft','cm',neutral=True,H=4,C=1,O=1,d_ref='Ball et al. 1970 ApJ 162, L203',lab_ref='Ball et al. 1970 ApJ 162, L203',notes=None,ice=True,ice_d_ref='Grim et al. 1991 A&A 243, 473',ice_l_ref='d\'Hendecourt & Allamandola 1986 A&A Sup. Ser. 64, 453',ppd=True,ppd_d_ref='Walsh et al. 2016 ApJL 823, L10',exgal=True,exgal_d_ref='Henkel et al. 1987 A&A 188, L1',exgal_sources='NGC 253, IC 342')
CH3CN = Molecule('methyl cyanide','CH3CN',1971,'CH3CN','Sgr A, Sgr B2','NRAO 36-ft','mm',neutral=True,H=3,C=2,N=1,d_ref='Solomon et al. 1971 ApJ 168, L107',lab_ref='Cord et al. 1968 Microwave Spectral Tables V5; Kessler et al. Phys Rev 79, 54',notes=None,ppd=True,ppd_d_ref='Oberg et al. 2015 Nature 520, 198',exgal=True,exgal_d_ref='Mauersberger et al. 1991 A&A 247, 307',exgal_sources='NGC 253')
NH2CHO = Molecule('formamide','NH2CHO',1971,'NH2CHO','Sgr B2','NRAO 140-ft','cm',neutral=True,H=3,C=1,O=1,N=1,d_ref='Rubin et al. 1971 ApJ 169, L39',lab_ref='Rubin et al. 1971 ApJ 169, L39',notes=None,exgal=True,exgal_d_ref='Muller et al. 2013 A&A 551, A109',exgal_sources='PKS 1830-211 LOS')
CH3SH = Molecule('methyl mercaptan','CH3SH',1979,'CH3SH','Sgr B2','Bell 7-m','mm',neutral=True,H=4,C=1,S=1,d_ref='Linke et al. 1979 ApJ 234, L139',lab_ref='Kilb 1955 JCP 23, 1736',notes=None)
C2H4 = Molecule('ethylene','C2H4',1981,'C2H4','IRC+10216','McMath Solar Telescope','IR',neutral=True,H=4,C=2,d_ref='Betz 1981 ApJ 244, L103',lab_ref='Lambeau et al. 1980 JMS 81, 227',notes=None)
C5H = Molecule('pentynylidyne radical','C5H',1986,'C5H','IRC+10216','IRAM','cm',radical=True,neutral=True,H=1,C=5,d_ref='Cernicharo et al. 1986 A&A 164, L1',lab_ref='Gottlieb et al. 1986 A&A 164, L5',notes='*See also Cernicharo et al. 1986 A&A 167, L5 and Cernicharo et al. 1987 A&A 172, L5')
CH3NC = Molecule('methyl isocyanide','CH3NC',1988,'CH3NC','Sgr B2','IRAM','mm',neutral=True,H=3,C=2,N=1,d_ref='Cernicharo et al. 1988 A&A 189, L1',lab_ref='Kukolich 1972 JCP 57, 869; Ring et al. 1947 Phys Rev 72, 1262',notes='*Confirmed in Remijan et al. 2005 ApJ 632, 333 and Gratier et al. 2013 557, A101')
HC2CHO = Molecule('propynal','HC2CHO',1988,'HC2CHO','TMC-1','NRAO 140-ft, Nobeyama','cm',neutral=True,H=2,C=3,O=1,d_ref='Irvine et al. 1988 ApJ 335, L89',lab_ref='Winnewisser 1973 JMS 46, 16',notes=None)
H2C4 = Molecule('butatrienylidene','H2C4',1991,'H2C4','IRC+10216','IRAM','mm',neutral=True,H=2,C=4,d_ref='Cernicharo et al. 1991 ApJ 368, L43',lab_ref='Killian et al. 1990 ApJ 365, L89',notes=None)
C5S = Molecule('pentacarbon monosulfide radical','C5S',1993,'C5S','IRC+10216','NRAO 140-ft','cm',radical=True,neutral=True,C=5,S=1,d_ref='Bell et al. 1993 ApJ 417, L37',lab_ref='Kasai et al. 1993 ApJ 410, L45; Gordon et al. 2001 ApJS 134, 311',notes='*Confirmed in Agúndez et al. 2014 A&A 570, A45')
HC3NHp = Molecule('protonated cyanoacetylene','HC3NH+',1994,'HC3NH+','TMC-1','Nobeyama','cm',cation=True,H=2,C=3,N=1,d_ref='Kawaguchi et al. 1994 ApJ 420, L95',lab_ref='Lee & Amano 1987 ApJ 323',notes=None)
C5N = Molecule('cyanobutadiynyl radical','C5N',1998,'C5N','TMC-1','Effelsberg 100-m, IRAM','cm',radical=True,neutral=True,C=5,N=1,d_ref='Guélin et al. 1998 A&A 355, L1',lab_ref='Kasai et al. 1997 ApJ 477, L65',notes=None)
HC4H = Molecule('diacetylene','HC4H',2001,'HC4H','CRL 618','ISO','IR',neutral=True,H=2,C=4,d_ref='Cernicharo et al. 2001 ApJ 546, L123',lab_ref='Arie & Johns 1992 JMS 155, 195',notes='*Confirmed in 2018 ApJ 852, 80',exgal=True,exgal_d_ref='Bernard-Salas et al. 2006 ApJ 652, L29',exgal_sources='SMP LMC 11')
HC4N = Molecule('','HC4N',2004,'HC4N','IRC+10216','IRAM','mm',radical=True,neutral=True,H=1,C=4,N=1,d_ref='Cernicharo et al. 2004 ApJ 615, L145',lab_ref='Tang et al. 1999 CPL 315, 69',notes=None)
cH2C3O = Molecule('cyclopropenone','c-H2C3O',2006,'c-H2C3O','Sgr B2','GBT','cm',neutral=True,cyclic=True,H=2,C=3,O=1,d_ref='Hollis et al. 2006 ApJ 642, 933',lab_ref='Benson et al. 1973 JACS 95, 2772; Guillemin et al. 1990 JMS 140, 190',notes=None)
CH2CNH = Molecule('ketenimine','CH2CNH',2006,'CH2CNH','Sgr B2','GBT','cm',neutral=True,H=3,C=2,N=1,d_ref='Lovas et al. 2006 ApJ 645, L137',lab_ref='Rodler et al. 1984 CPL 110, 447; Rodler et al. 1986 JMS 118, 267',notes=None)
C5Nm = Molecule('cyanobutadiynyl anion','C5N-',2008,'C5N-','IRC+10216','IRAM','mm',anion=True,C=5,N=1,d_ref='Cernicharo et al. 2008 ApJ 688, L83',lab_ref='Botschwina & Oswald 2008 JCP 129, 044305',notes=None)
HNCHCN = Molecule('E-cyanomethanimine','HNCHCN',2013,'HNCHCN','Sgr B2','GBT','cm',neutral=True,H=2,C=2,N=2,d_ref='Zaleski et al. 2013 ApJ 765, L9',lab_ref='Zaleski et al. 2013 ApJ 765, L9',notes=None)
SiH3CN = Molecule('silyl cyanide','SiH3CN',2014,'SiH3CN','IRC+10216','IRAM','mm',neutral=True,H=3,C=1,N=1,Si=1,d_ref='Agúndez et al. 2014 A&A 570, A45',lab_ref='Priem et al. 1998 JMS 191, 183',notes='*Confirmed in Cernicharo et al. 2017 A&A 606, L5')

six_atom_list = [CH3OH,CH3CN,NH2CHO,CH3SH,C2H4,C5H,CH3NC,HC2CHO,H2C4,C5S,HC3NHp,C5N,HC4H,HC4N,cH2C3O,CH2CNH,C5Nm,HNCHCN,SiH3CN]


#############################################################
#					Seven Atom Molecules					#
#############################################################

CH3CHO = Molecule('acetaldehyde','CH3CHO',1973,'CH3CHO','Sgr B2','NRAO 140-ft','cm',neutral=True,H=4,C=2,O=1,d_ref='Gottlieb 1973 Molecules in the Galactic Environment 181; Fourikis et al. 1974 Aust J Phys 27, 425; Gilmore et al. 1976 ApJ 204, 43',lab_ref='Kilb et al. 1957 JCP 26, 1695; Souter & Wood 1970 JCP 52, 674',notes=None,ice='Tentative',ice_d_ref='Schutte et al. 1999 A&A 343, 966',ice_l_ref='Schutte et al. 1999 A&A 343, 966',exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS')
CH3CCH = Molecule('methylacetylene','CH3CCH',1973,'CH3CCH','Sgr B2','NRAO 36-ft','mm',neutral=True,H=4,C=3,d_ref='Buhl & Snyder 1973 Molecules in the Galactic Environment 187',lab_ref='Trambarulo et al. 1950 JCP 18, 1613',notes=None,exgal=True,exgal_d_ref='Mauersberger et al. 1991 A&A 247, 307',exgal_sources='NGC 253, M82')
CH3NH2 = Molecule('methylamine','CH3NH2',1974,'CH3NH2','Sgr B2, Orion','Mitaka 6-m, NRAO 36-ft, Parkes 64-m','cm, mm',neutral=True,H=5,C=1,N=1,d_ref='Fourikis et al. 1974 ApJ 191, L139; Kaifu et al. 1974 ApJ 191, L135',lab_ref='Takagi & Kojima 1973 ApJ 181, L91',notes=None,exgal=True,exgal_d_ref='Muller et al. 2011 A&A 535, A103',exgal_sources='PKS 1830-211 LOS')
CH2CHCN = Molecule('vinylcyanide','CH2CHCN',1975,'CH2CHCN','Sgr B2','Parkes 64-m','cm',neutral=True,H=3,C=3,N=1,d_ref='Gardner & Winnewisser 1975 ApJ 195, L127',lab_ref='Gerry & Winnewisser 1973 JMS 48, 1',notes=None)
HC5N = Molecule('cyanodiacetylene','HC5N',1976,'HC5N','Sgr B2','Algonquin 46-m','cm',neutral=True,H=1,C=5,N=1,d_ref='Broten et al. 1976 ApJ 209, L143; Avery et al. 1976 ApJ 205 L173',lab_ref='Alexander et al. 1976 JMS 62, 175',notes=None,exgal='Tentative',exgal_d_ref='Aladro et al. 2015 A&A 579, A101',exgal_sources='NGC 253')
C6H = Molecule('hexatriynyl radical','C6H',1986,'C6H','TMC-1','Nobeyama','cm',radical=True,neutral=True,H=1,C=6,d_ref='Suzuki et al. 1986 PASJ 38, 911',lab_ref='Pearson et al. 1988 A&A 189, L13',notes=None)
cC2H4O = Molecule('ethylene oxide','c-C2H4O',1997,'c-C2H4O','Sgr B2','Haystack, Nobeyama, SEST 15-m','cm, mm',neutral=True,cyclic=True,H=4,C=2,O=1,d_ref='Dickens et al. 1997 ApJ 489, 753',lab_ref='Hirose 1974 ApJ 189, L145',notes=None)
CH2CHOH = Molecule('vinyl alcohol','CH2CHOH',2001,'CH2CHOH','Sgr B2','NRAO 12-m','mm',neutral=True,H=4,C=2,O=1,d_ref='Turner & Apponi 2001 ApJ 561, L207',lab_ref='Rodler 1985 JMS 114, 23; Kaushik 1977 CPL 49, 90',notes=None)
C6Hm = Molecule('hexatriynyl anion','C6H-',2006,'C6H-','TMC-1, IRC+10216','GBT','cm',anion=True,H=1,C=6,d_ref='McCarthy et al. 2006 ApJ 652, L141',lab_ref='McCarthy et al. 2006 ApJ 652, L141',notes='*First gas-phase molecular anion')
CH3NCO = Molecule('methyl isocyanate','CH3NCO',2015,'CH3NCO','Sgr B2, Orion','ARO 12-m, SMT','mm',neutral=True,H=3,C=2,O=1,N=1,d_ref='Halfen et al. 2015 ApJ 812, L5',lab_ref='Halfen et al. 2015 ApJ 812, L5',notes='*see also Cernicharo et al. 2016 A&A 587, L4')
HC5O = Molecule('butadiynylformyl radical','HC5O',2017,'HC5O','TMC-1','GBT','cm',radical=True,neutral=True,H=1,C=5,O=1,d_ref='McGuire et al. 2017 ApJ 843, L28',lab_ref='Mohamed et al. 2005 JCP 123, 234301',notes=None)

seven_atom_list = [CH3CHO,CH3CCH,CH3NH2,CH2CHCN,HC5N,C6H,cC2H4O,CH2CHOH,C6Hm,CH3NCO,HC5O]


#############################################################
#					Eight Atom Molecules					#
#############################################################

HCOOCH3 = Molecule('methyl formate','HCOOCH3',1975,'HCOOCH3','Sgr B2','Parkes 64-m, Effelsberg 100-m','cm',neutral=True,H=4,C=2,O=2,d_ref='Churchwell & Winnewisser 1975 A&A 45, 229; Brown et al. 1975 ApJ 197, L29',lab_ref='Brown et al. 1975 ApJ 197, L29',notes='*t-mf detected 2012 ApJ 755, 143',exgal=True,exgal_d_ref='Sewiło et al. 2018 ApJL 853, L19',exgal_sources='LMC')
CH3C3N = Molecule('methylcyanoacetylene','CH3C3N',1984,'CH3C3N','TMC-1','NRAO 140-ft','cm',neutral=True,H=3,C=4,N=1,d_ref='Broten et al. 1984 ApJ 276, L25',lab_ref='Moises et al. 1982 JMS 92, 497',notes=None)
C7H = Molecule('heptatriynylidyne radical','C7H',1997,'C7H','IRC+10216','IRAM','mm',radical=True,neutral=True,H=1,C=7,d_ref='Guélin et al. 1997 A&A 317, L1',lab_ref='Travers et al. 1996 ApJ 465, L77',notes=None)
CH3COOH = Molecule('acetic acid','CH3COOH',1997,'CH3COOH','Sgr B2','BIMA, OVRO','mm',neutral=True,H=4,C=2,O=2,d_ref='Mehringer et al. 1997 ApJ 480, L71',lab_ref='Tabor 1957 JCP 27, 974',notes=None)
H2C6 = Molecule('hexapentaenylidene','H2C6',1997,'H2C6','TMC-1','Goldstone 70-m','cm',neutral=True,H=2,C=6,d_ref='Langer et al. 1997 ApJ 480, L63',lab_ref='McCarthy et al. 1997 Science 275, 518',notes=None)
CH2OHCHO = Molecule('glycolaldehyde','CH2OHCHO',2000,'CH2OHCHO','Sgr B2','NRAO 12-m','mm',neutral=True,H=4,C=2,O=2,d_ref='Hollis et al. 2000 ApJ 540, L107',lab_ref='Marstokk & Mollendal 1973 J Mol Struct 16, 259',notes=None)
HC6H = Molecule('triacetylene','HC6H',2001,'HC6H','CRL 618','ISO','IR',neutral=True,H=2,C=6,d_ref='Cernicharo et al. 2001 ApJ 546, L123',lab_ref='Haas etal. 1994 JMS 167, 176',notes=None,exgal=True,exgal_d_ref='Bernard-Salas et al. 2006 ApJ 652, L29',exgal_sources='SMP LMC 11')
CH2CHCHO = Molecule('propenal','CH2CHCHO',2004,'CH2CHCHO','Sgr B2, G327.3-0.6','NRAO 12-m, SEST 15-m','cm',neutral=True,H=4,C=3,O=1,d_ref='Hollis et al. 2004 ApJ 610, L21',lab_ref='Winnewisser et al. 1975 Z Naturforsch 30, 1001',notes=None)
CH2CCHCN = Molecule('cyanoallene','CH2CCHCN',2006,'CH2CCHCN','TMC-1','GBT','cm',neutral=True,H=3,C=4,N=1,d_ref='Lovas et al. 2006 ApJ 637, L37',lab_ref='Bouche et al. 1973 J Mol Struct 18, 211',notes='*Also Chin et al. 2006 AIP Conf. Proc. 855, 149')
NH2CH2CN = Molecule('aminoacetonitrile','NH2CH2CN',2008,'NH2CH2CN','Sgr B2','IRAM, PdBI,  ATCA','mm',neutral=True,H=4,C=2,N=2,d_ref='Belloche et al. 2008 A&A 482, 179',lab_ref='Bogey et al. 1990 JMS 143, 180',notes=None)
CH3CHNH = Molecule('ethanimine','CH3CHNH',2013,'CH3CHNH','Sgr B2','GBT','cm',neutral=True,H=5,C=2,N=1,d_ref='Loomis et al. 2013 ApJL 765, L10',lab_ref='Loomis et al. 2013 ApJL 765, L10',notes=None)
CH3SiH3 = Molecule('methyl silane','CH3SiH3',2017,'CH3SiH3','IRC+10216','IRAM','mm',neutral=True,H=6,C=1,Si=1,d_ref='Cernicharo et al. 2017 A&A 606, L5',lab_ref='Wong et al. 1983 JMS 102, 89',notes=None)

eight_atom_list = [HCOOCH3,CH3C3N,C7H,CH3COOH,H2C6,CH2OHCHO,HC6H,CH2CHCHO,CH2CCHCN,NH2CH2CN,CH3CHNH,CH3SiH3]


#############################################################
#					Nine Atom Molecules						#
#############################################################

CH3OCH3 = Molecule('dimethyl ether','CH3OCH3',1974,'CH3OCH3','Orion','NRAO 36-ft, NRL 85-ft','cm, mm',neutral=True,H=6,C=2,O=1,d_ref='Snyder et al. 1974 ApJ 191, L79',lab_ref='Kasai & Myers JCP 30, 1096; Blukis et al. 1963 JCP 38, 2753',notes=None,exgal=True,exgal_d_ref='Qiu et al. 2018 A&A 613, A3; Sewiło et al. 2018 ApJL 853, L19',exgal_sources='NGC 1068, LMC')
CH3CH2OH = Molecule('ethanol','CH3CH2OH',1975,'CH3CH2OH','Sgr B2','NRAO 36-ft','mm',neutral=True,H=6,C=2,O=1,d_ref='Zukerman et al. 1975 ApJ 196, L99',lab_ref='Takano et al. 1986 JMS 26, 157',notes='*g-ethanol detected 1997 ApJ 480, 420')
CH3CH2CN = Molecule('ethyl cyanide','CH3CH2CN',1977,'CH3CH2CN','Orion, Sgr B2','NRAO 36-ft','mm',neutral=True,H=5,C=3,N=1,d_ref='Johnson et al. 1977 ApJ 218, 370',lab_ref='Johnson et al. 1977 ApJ 218, 370',notes=None)
HC7N = Molecule('cyanotriacetylene','HC7N',1977,'HC7N','TMC-1','Algonquin 46-m, Haystack','cm',neutral=True,H=1,C=7,N=1,d_ref='Kroto et al. 1977 Bull. Am. As. Soc. 9, 303',lab_ref='Kirby et al. 1980 JMS 83, 261',notes=None)
CH3C4H = Molecule('methyldiacetylene','CH3C4H',1984,'CH3C4H','TMC-1','Effelsberg 100-m, Haystack, NRAO 140-ft','cm',neutral=True,H=4,C=5,d_ref='Walmsley et al. 1984 A&A 134, L11',lab_ref='Heath et al. 1955 Faraday Discuss. 19, 38',notes=None)
C8H = Molecule('octatriynyl radical','C8H',1996,'C8H','IRC+10216','IRAM, Nobeyama','mm',radical=True,neutral=True,H=1,C=8,d_ref='Cernicharo & Guélin 1996 A&A 309, L27',lab_ref='Pauzat et al. 1991 ApJ 369, L13',notes=None)
CH3CONH2 = Molecule('acetamide','CH3CONH2',2006,'CH3CONH2','Sgr B2','GBT','cm',neutral=True,H=5,C=2,O=1,N=1,d_ref='Hollis et al. 2006 ApJ 643, L25',lab_ref='Suenram et al. 2001 JMS 208, 188',notes=None)
C8Hm = Molecule('octatriynyl anion','C8H-',2007,'C8H-','TMC-1, IRC+10216','GBT','cm',anion=True,H=1,C=8,d_ref='Brünken et al. 2007 ApJ 664, L43; Remijan et al. 2007 ApJ 664, L47',lab_ref='Gupta et al. 2007 ApJ 655, L57',notes=None)
CH2CHCH3 = Molecule('propylene','CH2CHCH3',2007,'CH2CHCH3','TMC-1','IRAM','mm',neutral=True,H=6,C=3,d_ref='Marcelino et al. 2007 ApJ 665, L127',lab_ref='Pearson et al. 1994 JMS 166, 120; Wlodarczak et al. 1994 JMS 167, 239',notes=None)
CH3CH2SH = Molecule('ethyl mercaptan','CH3CH2SH',2014,'CH3CH2SH','Orion','IRAM','mm',neutral=True,H=6,C=2,S=1,d_ref='Kolesniková et al. 2014 ApJ 784, L7',lab_ref='Kolesniková et al. 2014 ApJ 784, L7',notes=None)
HC7O = Molecule('hexadiynylformyl radical','HC7O',2017,'HC7O','TMC-1','GBT','cm',radical=True,neutral=True,H=1,C=7,O=1,d_ref='McGuire et al. 2017 ApJ 843, L28',lab_ref='Mohamed et al. 2005 JCP 123, 234301',notes='*Confirmed in Cordiner et al. 2017 ApJ 850, 194')

nine_atom_list = [CH3OCH3,CH3CH2OH,CH3CH2CN,HC7N,CH3C4H,C8H,CH3CONH2,C8Hm,CH2CHCH3,CH3CH2SH,HC7O]


#############################################################
#					Ten Atom Molecules						#
#############################################################

acetone = Molecule('acetone','(CH3)2CO',1987,'acetone','Sgr B2','IRAM, NRAO 140-ft, NRAO 12-m','cm, mm',neutral=True,H=6,C=3,O=1,d_ref='Combes et al. 1987 A&A 180, L13',lab_ref='Vacherand et al. 1986 JMS 118, 355',notes='*Confirmed in 2002 ApJ 578, 245')
HOCH2CH2OH = Molecule('ethylene glycol','HOCH2CH2OH',2002,'HOCH2CH2OH','Sgr B2','NRAO 12-m','mm',neutral=True,H=6,C=2,O=2,d_ref='Hollis et al. 2002 ApJ 571, L59',lab_ref='Christen et al. 1995 JMS 172, 57',notes='*aGg\' conformer in 2017 A&A 598, A59')
CH2CH2CHO = Molecule('propanal','CH2CH2CHO',2004,'CH2CH2CHO','Sgr B2','GBT','cm',neutral=True,H=6,C=3,O=1,d_ref='Hollis et al. 2004 ApJ 610, L21',lab_ref='Butcher & Wilson 1964 JCP 40, 1671',notes=None)
CH3C5N = Molecule('methylcyanodiacetylene','CH3C5N',2006,'CH3C5N','TMC-1','GBT','cm',neutral=True,H=4,C=6,N=1,d_ref='Snyder et al. 2006 ApJ 647, 412',lab_ref='Chen et al. 1998 JMS 192, 1',notes=None)
CH3CHCH2O = Molecule('propylene oxide','CH3CHCH2O',2016,'CH3CHCH2O','Sgr B2','GBT','cm',neutral=True,cyclic=True,H=6,C=3,O=1,d_ref='McGuire & Carroll et al. 2016 Science 352, 1449',lab_ref='McGuire & Carroll et al. 2016 Science 352, 1449',notes='*First chiral molecule')
CH3OCH2OH = Molecule('methoxymethanol','CH3OCH2OH',2017,'CH3OCH2OH','NGC 6334','ALMA','mm',neutral=True,H=6,C=2,O=2,d_ref='McGuire et al. 2017 ApJ 851, L46',lab_ref='Motiyenko et al. 2018 PCCP 20, 5509',notes=None)

ten_atom_list = [acetone,HOCH2CH2OH,CH2CH2CHO,CH3C5N,CH3CHCH2O,CH3OCH2OH]


#############################################################
#					Eleven Atom Molecules					#
#############################################################

HC9N = Molecule('cyanotetraacetylene','HC9N',1978,'HC9N','TMC-1','Algonquin 46-m, NRAO 140-ft','cm',neutral=True,H=1,C=9,N=1,d_ref='Broten et al. 1978 ApJ 223, L105',lab_ref='Iida et al. 1991 ApJ 371, L45',notes=None)
CH3C6H = Molecule('methyltriacetylene','CH3C6H',2006,'CH3C6H','TMC-1','GBT','cm',neutral=True,H=4,C=7,d_ref='Remijan et al. 2006 ApJ 643, L37',lab_ref='Alexander et al. 1978 JMS 70, 84',notes=None)
C2H5OCHO = Molecule('ethyl formate','C2H5OCHO',2009,'C2H5OCHO','Sgr B2','IRAM','mm',neutral=True,H=6,C=3,O=2,d_ref='Belloche et al. 2009 A&A 499, 215',lab_ref='Medvedev et al. 2009 ApJS 181, 433',notes=None)
CH3COOCH3 = Molecule('methyl acetate','CH3COOCH3',2013,'CH3COOCH3','Orion','IRAM','mm',neutral=True,H=6,C=3,O=2,d_ref='Tercero et al. 2013 ApJ 770, L13',lab_ref='Tudorie et al. 2011 JMS 269, 211',notes=None)

eleven_atom_list = [HC9N,CH3C6H,C2H5OCHO,CH3COOCH3]


#############################################################
#					Twelve Atom Molecules					#
#############################################################

C6H6 = Molecule('benzene','C6H6',2001,'C6H6','CRL 618','ISO','IR',neutral=True,cyclic=True,H=6,C=6,d_ref='Cernicharo et al. 2001 ApJ 546, L123',lab_ref='Lindenmayer et al. 1988 JMS 128 172',notes=None,exgal=True,exgal_d_ref='Bernard-Salas et al. 2006 ApJ 652, L29',exgal_sources='SMP LMC 11')
nC3H7CN = Molecule('n-propyl cyanide','n-C3H7CN',2009,'n-C3H7CN','Sgr B2','IRAM','mm',neutral=True,H=7,C=4,N=1,d_ref='Belloche et al. 2009 A&A 499, 215',lab_ref='Belloche et al. 2009 A&A 499, 215',notes=None)
iC3H7CN = Molecule('isopropyl cyanide','i-C3H7CN',2014,'i-C3H7CN','Sgr B2','ALMA','mm',neutral=True,H=7,C=4,N=1,d_ref='Belloche et al. 2014 Science 345, 1584',lab_ref='Muller et al. 2011 JMS 267, 100',notes=None)

twelve_atom_list = [C6H6,nC3H7CN,iC3H7CN]


#############################################################
#					Thirteen Atom Molecules					#
#############################################################

cC6H5CN = Molecule('benzonitrile','c-C6H5CN',2018,'C6H5CN','TMC-1','GBT','cm',neutral=True,cyclic=True,H=5,C=7,N=1,d_ref='McGuire et al. 2018 Science 359, 202',lab_ref='Wohlfart et al. 2008 JMS 247, 119',notes=None)

thirteen_atom_list = [cC6H5CN]


#############################################################
#					Fullerene Molecules						#
#############################################################

C60 = Molecule('buckminsterfullerene','C60',2010,'C60','TC 1, NGC 7023','Spitzer','IR',neutral=True,C=60,d_ref='Cami et al. 2010 Science 329, 1180',lab_ref='Nemes et al. 1994 CPL 218, 295',notes='*See also Sellgren et al. 2010 ApJ 722, L54 and Werner 2004b, Sellgren 2007 therein')
C60p = Molecule('buckminsterfullerene cation','C60+',2013,'C60+','NGC 7023','Spitzer','IR',cation=True,C=60,d_ref='Berné et al. 2013 A&A 550, L4',lab_ref='Kern et al. 2013 JPCA 117, 8251',notes='*See also Campbell et al. 2015 Nature 523, 322')
C70 = Molecule('rugbyballene','C70',2010,'C70','TC 1','Spitzer','IR',neutral=True,C=70,d_ref='Cami et al. 2010 Science 329, 1180',lab_ref='Nemes et al. 1994 CPL 218, 295',notes=None)

fullerene_list = [C60,C60p,C70]


#############################################################
#							Full List						#
#############################################################


full_list = [CH,CN,CHp,OH,CO,H2,SiO,CS,SO,SiS,NS,C2,NO,HCl,NaCl,AlCl,KCl,AlF,PN,SiC,CP,NH,SiN,SOp,COp,HF,N2,CFp,PO,O2,AlO,CNm,OHp,SHp,HClp,SH,TiO,ArHp,NSp,H2O,HCOp,HCN,OCS,HNC,H2S,N2Hp,C2H,SO2,HCO,HNO,HCSp,HOCp,SiC2,C2S,C3,CO2,CH2,C2O,MgNC,NH2,NaCN,N2O,MgCN,H3p,SiCN,AlNC,SiNC,HCP,CCP,AlOH,H2Op,H2Clp,KCN,FeCN,HO2,TiO2,CCN,SiCSi,S2H,HCS,HSC,NCO,NH3,H2CO,HNCO,H2CS,C2H2,C3N,HNCS,HOCOp,C3O,lC3H,HCNHp,H3Op,C3S,cC3H,HC2N,H2CN,SiC3,CH3,C3Nm,PH3,HCNO,HOCN,HSCN,HOOH,lC3Hp,HMgNC,HCCO,CNCN,HC3N,HCOOH,CH2NH,NH2CN,H2CCO,C4H,SiH4,cC3H2,CH2CN,C5,SiC4,H2CCC,CH4,HCCNC,HNCCC,H2COHp,C4Hm,CNCHO,HNCNH,CH3O,NH3Dp,H2NCOp,NCCNHp,CH3Cl,CH3OH,CH3CN,NH2CHO,CH3SH,C2H4,C5H,CH3NC,HC2CHO,H2C4,C5S,HC3NHp,C5N,HC4H,HC4N,cH2C3O,CH2CNH,C5Nm,HNCHCN,SiH3CN,CH3CHO,CH3CCH,CH3NH2,CH2CHCN,HC5N,C6H,cC2H4O,CH2CHOH,C6Hm,CH3NCO,HC5O,HCOOCH3,CH3C3N,C7H,CH3COOH,H2C6,CH2OHCHO,HC6H,CH2CHCHO,CH2CCHCN,NH2CH2CN,CH3CHNH,CH3SiH3,CH3OCH3,CH3CH2OH,CH3CH2CN,HC7N,CH3C4H,C8H,CH3CONH2,C8Hm,CH2CHCH3,CH3CH2SH,HC7O,acetone,HOCH2CH2OH,CH2CH2CHO,CH3C5N,CH3CHCH2O,CH3OCH2OH,HC9N,CH3C6H,C2H5OCHO,CH3COOCH3,C6H6,nC3H7CN,iC3H7CN,cC6H5CN,C60,C60p,C70]

natoms_list = [two_atom_list,three_atom_list,four_atom_list,five_atom_list,six_atom_list,seven_atom_list,eight_atom_list,nine_atom_list,ten_atom_list,eleven_atom_list,twelve_atom_list,thirteen_atom_list,fullerene_list]

elements_list = ['H','C','O','N','S','P','Si','Cl','F','Mg','Na','Al','K','Fe','Ti','Ar']

#############################################################
#					 	Update Stats			 			#
#############################################################

for x in full_list:

	x.natoms = x.H + x.C + x.O + x.N + x.S + x.P + x.Si + x.Cl + x.F + x.Mg + x.Na + x.Al + x.K + x.Fe + x.Ti + x.Ar
	x.mass = x.H + x.C*12 + x.O*16 + x.N*14 + x.S*32 + x.P*31 + x.Si*28 + x.Cl*35 + x.F*19 + x.Mg*24 + x.Na*23 + x.Al*27 + x.K*39 + x.Fe*56 + x.Ti*48 + x.Ar*36
	
	#for molecules containing C, and only H,O,N, and halogens, calculate the degree of unsaturation
	
	if all(int(i) is 0 for i in [x.Si,x.Mg,x.Na,x.Al,x.K,x.Fe,x.Ti,x.Ar,x.P,x.Si]) and x.C != 0:
	
		x.du = 1 + 0.5*(x.H*-1 + x.C*2 + x.N*1 + x.Cl*-1 + x.F*-1)
		
	else:
	
		x.du = None

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

source_tag_list = [AFGL890LOS,AFGL961LOS,AFGL989LOS,B1b,CRL2688,CRL618,CasALOS,CrabNebula,CygnusOB212LOS,DR21,DR21LOS,DR21OH,GL2136LOS,GalacticCenter,HD124314LOS,HD27778LOS,HorseheadPDR,IC443G,IRAS16293,IRC10216,K350,L134,L1527,L1544,L183,L483,LOSCloud,Lupus1A,M17LOS,M17SW,M3LOS,NGC2024,NGC2024LOS,NGC2264,NGC6334,NGC6334LOS,NGC7023,NGC7027,NGC7538,NGC7538LOS,Orion,OrionBar,SgrA,SgrALOS,SgrB2,SgrB2LOS,TC1,TMC1,VYCaMaj,W3,W3OH,W31LOS,W33LOS,W43LOS,W44LOS,W49,W49LOS,W51,W51LOS,XiPerLOS,rhoOphA]

#############################################################
#					 	Update Stats			 			#
#############################################################

for x in source_tag_list:

	for y in full_list:
	
		mol_sources = y.sources.split(',')
		
		for source in mol_sources:
		
			if x.name == source.strip():
		
				x.detects +=1			
	
#############################################################
#						Functions	 						#
#############################################################	

#find_tag_name(y) suggests a list of possible entries based on the user-input molecular formula y.  This is still buggy and not implemented.

# def find_tag_name(y):
# 
# 	elements_dict = dict(
# 	
# 	H = 0,
# 	C = 0,
# 	O = 0,
# 	N = 0,
# 	S = 0,
# 	P = 0,
# 	Si = 0,
# 	Cl = 0,
# 	F = 0, 
# 	Mg = 0,
# 	Na = 0,
# 	Al = 0,
# 	K = 0,
# 	Fe = 0,
# 	Ti = 0,
# 	Ar = 0
# 	
# 	)
# 	
# 	for x in range(len(y)):
# 	
# 		#if it's a number, we move on:
# 		
# 		if y[x].isalpha() == False:
# 		
# 			continue
# 			
# 		#if it's a lower case letter, we move on:
# 		
# 		if y[x].islower() == True:
# 		
# 			continue
# 	
# 		#if it's the last character and it's a capital letter, then there's only one of that element, and we proceed to add it to the dictionary
# 		
# 		elif x == len(y)-1 and y[x].isalpha() == True and y[x].islower() == False:
# 			
# 			element = y[x]
# 			natoms = 1
# 	
# 		#if it's not the last character and it's a capital letter...
# 		
# 		elif y[x].isalpha() == True and y[x].islower() == False:
# 				
# 			#if the next character is a lower case, then we've got a two letter element
# 			
# 			if y[x+1].islower() == True:
# 			
# 				#make the element name out of the two letters
# 				
# 				element = y[x]+y[x+1]
# 				
# 				#if we aren't at the end of the string
# 				
# 				if (x+2 in range(len(y))):
# 				
# 					#if the letter after the lowercase is another letter, then it's just one of that element
# 					
# 					if y[x+2].isalpha() == True:
# 				
# 						natoms = 1
# 					
# 					#otherwise it's a number and that's how many we have
# 					
# 					else:
# 				
# 						natoms = int(y[x+2])
# 				
# 				#if we are at the end of the string, then it's just one of that element
# 						
# 				else:
# 				
# 					natoms = 1
# 					
# 			#Otherwise, if it was just a capital letter, and the next character is a number
# 				
# 			elif y[x+1].isalpha() == False:
# 				
# 				#the element is that capital letter
# 				
# 				element = y[x]
# 				
# 				#and we get the number from the next character, unless the character after that is also a number
# 				
# 				if (x+2 in range(len(y))):
# 				
# 					#let's check to make sure something with 100+ atoms hasn't been entered:
# 					
# 					if (x+3 in range(len(y))):
# 					
# 						if y[x+3].isalpha() == False:
# 						
# 							natoms = int(y[x+1:x+4])
# 						
# 							print('Sorry, there are no molecules with {}(!) {} atoms in the database.' .format(natoms,element))
# 							return
# 					
# 					if y[x+2].isalpha() == False:
# 					
# 						natoms = int(y[x+1:x+3])
# 						
# 				else:
# 				
# 					natoms = int(y[x+1])
# 			
# 			#if the next character is not a number, then we just have one of that element and we move on
# 				
# 			elif y[x+1].isalpha() == True:
# 			
# 				element = y[x]
# 				
# 				natoms = 1
# 		
# 		if element not in elements_dict:
# 			
# 			print('Sorry, there are no molecules containing {} in the database.' .format(element))
# 			return
# 				
# 		else:
# 			
# 			elements_dict[element] += natoms
# 	
# 	possible_hits = []
# 	
# 	for x in full_list:
# 	
# 		test_array = [
# 		
# 		x.H == elements_dict['H'],
# 		x.C == elements_dict['C'],
# 		x.O == elements_dict['O'],
# 		x.N == elements_dict['N'],
# 		x.S == elements_dict['S'],
# 		x.P == elements_dict['P'],
# 		x.Si == elements_dict['Si'],
# 		x.Cl == elements_dict['Cl'],
# 		x.F == elements_dict['F'],
# 		x.Mg == elements_dict['Mg'],
# 		x.Na == elements_dict['Na'],
# 		x.Al == elements_dict['Al'],
# 		x.K == elements_dict['K'],
# 		x.Fe == elements_dict['Fe'],
# 		x.Ti == elements_dict['Ti'],
# 		x.Ar == elements_dict['Ar'],
# 		
# 		]		
# 		
# 		if all(test_array) == True:
# 		
# 			possible_hits.append(x.label)
# 				
# 	print('\n[Tag] Molecule Name (Formula)')
# 	print('-----------------------------')
# 
# 	for x in possible_hits:
# 	
# 		tag = x
# 		
# 		if tag[-1] == '-':
# 		
# 			tag = tag[:-1] + 'm'
# 			
# 		if tag[-1] == '+':
# 		
# 			tag = tag[:-1] + 'p'
# 			
# 		tag = tag.replace('-','')
# 	
# 		print('[{}] {} ({})' .format(tag,globals()[tag].name,globals()[tag].formula))

#summary(y) prints a nicely-formatted summary of the database entry for that molecule or source

def summary(y):

	if isinstance(y,Molecule):
	
		n_dash = len(y.name) + len(y.formula) +3 
	
		dashes = '-' * n_dash

		print('\n' + dashes)
		print('{} ({})' .format(y.name,y.formula))
		print(dashes + '\n')
		print('Atoms:\t{}' .format(y.natoms))
		print('Mass:\t{} amu' .format(y.mass))
		print('Year Detected:\t{}' .format(y.year))
		print('Source(s):\t{}' .format(y.sources))
		print('Telescope(s) Used:\t{}' .format(y.telescopes))
	
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
	
#refs(y) prints a nicely-formatted list of references and notes for molecule y

def refs(y):

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
		
		#to be enabled later if lab references are cataloged..
		#print('[Lab] {}' .format(y.exo_l_ref))		

#output_summary(y,filename=None) writes out an ascii file containing the output of summary(y) for molecule y or a list of molecules y.  A filename can optionally be specified.

def output_summary(y,filename=None):

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

#mols_with_elements generates a tab-delimited ascii file with the number of molecules detected containing each element

def mols_with_elements():

	elements_dict = dict((element,0) for element in elements_list)
	
	for x in full_list:
	
		if x.H != 0:
		
			elements_dict['H'] += 1
			
		if x.C != 0:
		
			elements_dict['C'] += 1
			
		if x.O != 0:
		
			elements_dict['O'] += 1
			
		if x.N != 0:
		
			elements_dict['N'] += 1
			
		if x.S != 0:
		
			elements_dict['S'] += 1
			
		if x.P != 0:
		
			elements_dict['P'] += 1
			
		if x.Si != 0:
		
			elements_dict['Si'] += 1
			
		if x.Cl != 0:
		
			elements_dict['Cl'] += 1
			
		if x.F != 0:
		
			elements_dict['F'] += 1
			
		if x.Mg != 0:
		
			elements_dict['Mg'] += 1
			
		if x.Al != 0:
		
			elements_dict['Al'] += 1
			
		if x.K != 0:
		
			elements_dict['K'] += 1
			
		if x.Fe != 0:
		
			elements_dict['Fe'] += 1
			
		if x.Ti != 0:
		
			elements_dict['Ti'] += 1													
			
		if x.Ar != 0:
		
			elements_dict['Ar'] += 1	
			
		if x.Na != 0:
		
			elements_dict['Na'] += 1			
			
	with open('mols_with_elements.txt','w') as output:
	
		output.write('Element\tMolecules\n')
		
		for x in sorted(elements_dict):
		
			output.write('{}\t{}\n' .format(x,elements_dict[x]))		

#molecule_types prints the number of molecules of each type that have been detected, in total.  Set output=True to have it return a dictionary with these values.

def molecule_types(output=False):

	types_dict = {}
	
	types_dict['neutral'] = 0
	types_dict['cation'] = 0
	types_dict['anion'] = 0
	types_dict['radical'] = 0
	types_dict['cyclic'] = 0

	for x in full_list:
	
		types_dict['neutral'] += x.neutral
		types_dict['cation'] += x.cation
		types_dict['anion'] += x.anion
		types_dict['radical'] += x.radical
		types_dict['cyclic'] += x.cyclic
	
	for x in sorted(types_dict):
	
		print('{}: {}' .format(x,types_dict[x]))
		
	if output == True:	
	
		return types_dict	
		
#make_source_types_list finds all the unique source types among cataloged sources and returns it as a list

def make_source_types_list():

	types_list = []
	
	for source in source_tag_list:
	
		if source.type not in types_list:
		
			types_list.append(source.type)
			
	return types_list

#molecule_types_by_source_type produces a tab delimited ascii file where the number of each molecule detected in a given source type, with a given attribute, is totaled.

def molecule_types_by_source_type():

	types_list = make_source_types_list()
	
	types_dict = dict((type,[0,0,0,0,0]) for type in types_list)
	
	for x in full_list:
	
		#get the sources for a molecule
	
		sources = x.sources.split(',')
		
		#initialize a list of source types so we don't overcount some source types
		
		types_tracker = make_source_types_list()
		
		for y in sources:
		
			z = y.strip()
			
			for a in source_tag_list:
			
				if a.name == z:
				
					#check if that type has yet to be used, if it hasn't, go ahead and credit it with a detection
				
					if a.type in types_tracker:					
								
						types_dict[a.type][0] += x.neutral
						types_dict[a.type][1] += x.cation
						types_dict[a.type][2] += x.anion
						types_dict[a.type][3] += x.radical
						types_dict[a.type][4] += x.cyclic
						
						#remove that source type from the list, so it doesn't get credited again for this molecule
						
						types_tracker.remove(a.type)							
					
	with open('molecule_types_by_source_type.txt','w') as output:
	
		output.write('SourceType\tNeutral\tCation\tAnion\tRadical\tCyclic\n')
	
		for x in sorted(types_dict):
		
			output.write('{}\t{}\t{}\t{}\t{}\t{}\n' .format(x,types_dict[x][0],types_dict[x][1],types_dict[x][2],types_dict[x][3],types_dict[x][4]))
		
#mols_in_source_type generates a tab-delimited ascii file of generalized source type and total molecules detected there

def mols_in_source_type():

	types_list = make_source_types_list()
	
	types_dict = dict((type,0) for type in types_list)

	for x in full_list:
	
		#get the sources for a molecule
	
		sources = x.sources.split(',')
		
		#initialize a list of source types so we don't overcount some source types
		
		types_tracker = make_source_types_list()
		
		for y in sources:
		
			z = y.strip()
			
			for a in source_tag_list:
			
				if a.name == z:
					
					#check if that type has yet to be used, if it hasn't, go ahead and credit it with a detection
				
					if a.type in types_tracker:					
			
						types_dict[a.type] += 1

						#remove that source type from the list, so it doesn't get credited again for this molecule
						
						types_tracker.remove(a.type)						
		
	with open('mols_in_source_type.txt','w') as output:
	
		output.write('Source\tDetections\n')
		
		for x in sorted(types_dict):
		
			output.write('{}\t{}\n' .format(x,types_dict[x]))

#mols_in_sources generates a tab-delimited ascii file of detections in individual sources, with the exception of LOS clouds, which it lumps together

def mols_in_sources():

	sources_dict = dict((tag,0) for tag in source_tag_list if tag.type != 'LOS Cloud')
	
	los_count = 0
	
	for x in source_tag_list:
	
		if x.type != 'LOS Cloud':
		
			sources_dict[x] += x.detects

	for x in full_list:
	
		#get the sources for a molecule
	
		sources = x.sources.split(',')
		
		#check variable to prevent double counting LOS for a molecule
		
		los_counted = False
		
		for y in sources:
		
			z = y.strip()
			
			for a in source_tag_list:
			
				if a.name == z:
				
					if a.type == 'LOS Cloud' and los_counted == False:					
			
						los_count += 1

						los_counted = True				
			
	with open('mols_in_sources.txt','w') as output:
	
		output.write('Source\tDetections\n')
	
		for x in sources_dict: 
		
			output.write('{}\t{}\n' .format(x.name,sources_dict[x]))
				
		output.write('LOS Clouds\t{}' .format(los_count))

#mols_in_source_latex generates a latex table for the number of molecules detected in each source

def mols_in_source_latex():

	sources_dict = dict((tag,0) for tag in source_tag_list if tag.type != 'LOS Cloud')
	
	los_count = 0
	
	for x in source_tag_list:
		
		if x.type != 'LOS Cloud':
		
			sources_dict[x] += x.detects
			
	for x in full_list:
	
		#get the sources for a molecule
	
		sources = x.sources.split(',')
		
		#check variable to prevent double counting LOS for a molecule
		
		los_counted = False
		
		for y in sources:
		
			z = y.strip()
			
			for a in source_tag_list:
			
				if a.name == z:
				
					if a.type == 'LOS Cloud' and los_counted == False:					
			
						los_count += 1

						los_counted = True		
			
	sources_list = []
	
	for x in sources_dict:
	
		sources_list.append([x.name,sources_dict[x]])
		
	sources_list.append(['LOS Clouds',los_count])
	
	sources_list = sorted(sources_list, key=itemgetter(1), reverse=True)			
			
	with open('mols_in_sources.tex','w') as output:
	
		output.write('\\begin{table}[h!]\n\\centering\n')
		output.write('\\caption{Total number of detections for each source listed in \\S\\ref{known}.}\n')
		output.write('\\begin{tabular*}{\\columnwidth}{l @{\extracolsep{\\fill}}  c @{\extracolsep{\\fill}}  l @{\extracolsep{\\fill}}  c }\n')
		output.write('\\hline\\hline\n')
		output.write('Source\t&\t\\#\t&Source\t&\t\\# \\\\\n')
		output.write('\\hline\n')
	
		x = len(sources_list)
		
		z = int(0.5*x)
		
		for y in range(z):
		
				source_1 = sources_list[y][0]
				source_2 = sources_list[y+z][0]
				
				count_1 = sources_list[y][1]
				count_2 = sources_list[y+z][1]
		
				output.write('{}\t&\t{}\t&\t{}\t&\t{}\\\\\n' .format(source_1,count_1,source_2,count_2))
				
		output.write('\\hline\n\\end{tabular*}\n\\label{detects_by_source}\n\\end{table}')												

#du_by_source_type returns a tab-delimited ascii file with the degrees of unsaturation for each generalized source type, disregarding the fullerenes, and dropping the top and bottom 10% of values (or highest and lowest, in case of small numbers)

def du_by_source_type():

	types_list = make_source_types_list()
	
	types_dict = dict((type,[]) for type in types_list)
	
	for x in full_list:
	
		#get the sources for a molecule
	
		sources = x.sources.split(',')
		
		#initialize a list of source types so we don't overcount some source types
		
		types_tracker = make_source_types_list()
		
		for y in sources:
		
			z = y.strip()
			
			for a in source_tag_list:
			
				if a.name == z:
				
					if x.du != None and x.du < 20:
					
						#check if that type has yet to be used, if it hasn't, go ahead and credit it with a detection
					
						if a.type in types_tracker:					
				
							types_dict[a.type].append(x.du)

							#remove that source type from the list, so it doesn't get credited again for this molecule
							
							types_tracker.remove(a.type)					
					
	for x in types_dict:
		
		du_values = sorted(types_dict[x])
		
		n = len(du_values)
		
		#outliers = ceil(n*0.1)
		
		#types_dict[x] = du_values[outliers:n-outliers]		
		
		types_dict[x] = du_values
	
	write_dict = {}
	
	for x in types_dict:
	
		du_str = ''
		
		for y in types_dict[x]:
		
			du_str += '{}\t' .format(y)
			
		du_str.strip('\t')
		
		du_str += '\n'
		
		write_dict[x] = du_str
		
						
					
	with open('du_by_source_type.txt','w') as output:
	
		output.write('SourceType\tDU_vals\n')
	
		for x in sorted(write_dict):
		
			output.write('{}\t{}' .format(x,write_dict[x]))	

#mass_by_source_type returns a tab-delimited ascii file with the masses for each generalized source type, disregarding the fullerenes, and dropping the top and bottom 10% of values (or highest and lowest, in case of small numbers)

def mass_by_source_type():

	types_list = make_source_types_list()
	
	types_dict = dict((type,[]) for type in types_list)
	
	for x in full_list:
	
		#get the sources for a molecule
	
		sources = x.sources.split(',')
		
		#initialize a list of source types so we don't overcount some source types
		
		types_tracker = make_source_types_list()
		
		for y in sources:
		
			z = y.strip()
			
			for a in source_tag_list:
			
				if a.name == z:
				
					if x.C < 20:
				
						#check if that type has yet to be used, if it hasn't, go ahead and credit it with a detection
					
						if a.type in types_tracker:

							types_dict[a.type].append(x.mass)
							
							#remove that source type from the list, so it doesn't get credited again for this molecule
							
							types_tracker.remove(a.type)							
					
					
	for x in types_dict:
		
		mass_values = sorted(types_dict[x])
		
		n = len(mass_values)
		
		outliers = ceil(n*0.1)
		
		types_dict[x] = mass_values[outliers:n-outliers]		
	
	write_dict = {}
	
	for x in types_dict:
	
		mass_str = ''
		
		for y in types_dict[x]:
		
			mass_str += '{}\t' .format(y)
			
		mass_str.strip('\t')
		
		mass_str += '\n'
		
		write_dict[x] = mass_str						
					
	with open('mass_by_source_type.txt','w') as output:
	
		output.write('SourceType\tMasses\n')
	
		for x in sorted(write_dict):
		
			output.write('{}\t{}' .format(x,write_dict[x]))	

#atoms_by_source_type returns a tab-delimited ascii file with the number of atoms in each molecule for each generalized source type, disregarding the fullerenes, and dropping the top and bottom 10% of values (or highest and lowest, in case of small numbers)

def atoms_by_source_type():

	types_list = make_source_types_list()
	
	types_dict = dict((type,[]) for type in types_list)
	
	for x in full_list:
	
		sources = x.sources.split(',')
		
		for y in sources:
		
			z = y.strip()
			
			for a in source_tag_list:
			
				if a.name == z:
				
					if x.natoms < 60:
				
						types_dict[a.type].append(x.natoms)
					
					
	for x in types_dict:
		
		natoms_values = sorted(types_dict[x])
		
		n = len(natoms_values)
		
		outliers = ceil(n*0.1)
		
		types_dict[x] = natoms_values[outliers:n-outliers]		
	
	write_dict = {}
	
	for x in types_dict:
	
		natoms_str = ''
		
		for y in types_dict[x]:
		
			natoms_str += '{}\t' .format(y)
			
		natoms_str.strip('\t')
		
		natoms_str += '\n'
		
		write_dict[x] = natoms_str						
					
	with open('atoms_by_source_type.txt','w') as output:
	
		output.write('SourceType\tNAtoms\n')
	
		for x in sorted(write_dict):
		
			output.write('{}\t{}' .format(x,write_dict[x]))	
			
#mass_per_atom_by_source_type returns a tab-delimited ascii file with the mass per atom in each molecule for each generalized source type, disregarding the fullerenes, and dropping the top and bottom 10% of values (or highest and lowest, in case of small numbers)

def mass_per_atom_by_source_type():

	types_list = make_source_types_list()
	
	types_dict = dict((type,[]) for type in types_list)
	
	for x in full_list:
	
		sources = x.sources.split(',')
		
		for y in sources:
		
			z = y.strip()
			
			for a in source_tag_list:
			
				if a.name == z:
				
					if x.natoms < 60:
				
						mass_per = x.mass/x.natoms
				
						types_dict[a.type].append(mass_per)
					
					
	for x in types_dict:
		
		natoms_values = sorted(types_dict[x])
		
		n = len(natoms_values)
		
		outliers = ceil(n*0.1)
		
		types_dict[x] = natoms_values[outliers:n-outliers]		
	
	write_dict = {}
	
	for x in types_dict:
	
		natoms_str = ''
		
		for y in types_dict[x]:
		
			natoms_str += '{}\t' .format(y)
			
		natoms_str.strip('\t')
		
		natoms_str += '\n'
		
		write_dict[x] = natoms_str						
					
	with open('mass_per_atom_by_source_type.txt','w') as output:
	
		output.write('SourceType\tMassperAtom\n')
	
		for x in sorted(write_dict):
		
			output.write('{}\t{}' .format(x,write_dict[x]))				
			
#wave_by_source_type returns a tab-delimited ascii file with the number of molecules detected in each wavelength range for each generalized source type

def wave_by_source_type():

	types_list = make_source_types_list()
	
	'''
	The list for each source type is, in order:
	
	[0]cm
	[1]mm
	[2]sub-mm
	[3]IR
	[4]Vis
	[5]UV
	
	'''
	
	types_dict = dict((type,[0,0,0,0,0,0]) for type in types_list)
	
	for x in full_list:
	
		#get the sources for a molecule
	
		sources = x.sources.split(',')
		
		#initialize a list of source types so we don't overcount some source types
		
		cm_types_tracker = make_source_types_list()
		mm_types_tracker = make_source_types_list()
		sub_mm_types_tracker = make_source_types_list()
		ir_types_tracker = make_source_types_list()
		vis_types_tracker = make_source_types_list()
		uv_types_tracker = make_source_types_list()
		
		for y in sources:
		
			z = y.strip()
			
			for a in source_tag_list:
			
				if a.name == z:
					
					waves = x.wavelengths.split(',')
					
					for b in waves:
					
						if b.strip() == 'cm':
						
							if a.type in cm_types_tracker:
							
								types_dict[a.type][0] += 1
								
								cm_types_tracker.remove(a.type)
						
						if b.strip() == 'mm':
						
							if a.type in mm_types_tracker:
							
								types_dict[a.type][1] += 1
								
								mm_types_tracker.remove(a.type)
							
						if b.strip() == 'sub-mm':
						
							if a.type in sub_mm_types_tracker:
							
								types_dict[a.type][2] += 1
								
								sub_mm_types_tracker.remove(a.type)	
							
						if b.strip() == 'IR':
						
							if a.type in ir_types_tracker:
							
								types_dict[a.type][3] += 1
								
								ir_types_tracker.remove(a.type)
							
						if b.strip() == 'Vis':
						
							if a.type in vis_types_tracker:
							
								types_dict[a.type][4] += 1
								
								vis_types_tracker.remove(a.type)
							
						if b.strip() == 'UV':
						
							if a.type in uv_types_tracker:
							
								types_dict[a.type][5] += 1
								
								uv_types_tracker.remove(a.type)
										
	with open('wave_by_source_type.txt','w') as output:
	
		output.write('SourceType\tcm\tmm\tsub-mm\tIR\tVis\tUV\n')
	
		for x in sorted(types_dict):
		
			output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n' .format(x,types_dict[x][0],types_dict[x][1],types_dict[x][2],types_dict[x][3],types_dict[x][4],types_dict[x][5]))	

#mols_in_source returns all the molecules detected in source y

def mols_in_source(y):

	mol_list = []

	for x in full_list:
	
		if y in x.sources:
		
			mol_list.append(x.formula)
			
	return mol_list

def count(x,y):

	i = 0

	for z in full_list:
	
		if getattr(z,x) == y:
		
			i += 1
			
	return i
	
def make_eight_more_table():

	nlines = max(count('natoms',8),count('natoms',9))
	
	table_list = []
	
	eight_more_list = [eight_atom_list,nine_atom_list,ten_atom_list,eleven_atom_list,twelve_atom_list,thirteen_atom_list,fullerene_list]

	for x in range(nlines):
	
		table_line = ''
		
		for r in eight_more_list:
		
			if x < len(r):
			
				label = getattr(r[x],'label')
				formula = getattr(r[x],'formula')
				
				table_line += '\\hyperref[{}]{{\ce{{{}}}}}\t&\t' .format(label,formula)
				
			else:
			
				table_line += '\t&\t'
				
		table_line = table_line[:-2]
		table_line += '\\\\\n'		
		
		table_list.append(table_line)
		
	with open('latex_eight_more_table.txt','w') as output:
	
		for x in table_list:
		
			output.write(x)
			
def make_two_seven_table():

	table_list = []
	
	i = int(len(two_atom_list)/2)+1
	
	two_atom_list_1 = two_atom_list[:i]
	two_atom_list_2 = two_atom_list[i:]
	
	j = int(len(three_atom_list)/2)+1
	
	three_atom_list_1 = three_atom_list[:j]
	three_atom_list_2 = three_atom_list[j:]	
	
	nlines = max(len(two_atom_list_1),len(two_atom_list_2),len(three_atom_list_1),len(three_atom_list_2),len(four_atom_list),len(five_atom_list),len(six_atom_list),len(seven_atom_list))
	
	two_seven_list = [two_atom_list_1,two_atom_list_2,three_atom_list_1,three_atom_list_2,four_atom_list,five_atom_list,six_atom_list,seven_atom_list]

	for x in range(nlines):
	
		table_line = ''
		
		for r in two_seven_list:
		
			if x < len(r):
			
				label = getattr(r[x],'label')
				formula = getattr(r[x],'formula')
				
				table_line += '\\hyperref[{}]{{\ce{{{}}}}}\t&\t' .format(label,formula)
				
			else:
			
				table_line += '\t&\t'
				
		table_line = table_line[:-2]
		table_line += '\\\\\n'		
		
		table_list.append(table_line)
		
	with open('latex_two_seven_table.txt','w') as output:
	
		output.write('\\begin{table*}\n')
		output.write('\\centering\n')
		output.write('\\caption{List of detected interstellar molecules with two to seven atoms, categorized by number of atoms, and vertically ordered by detection year.  Column headers and molecule formulas are in-document hyperlinks in most PDF viewers.}\n')
		output.write('\\begin{tabular*}{\\textwidth}{l l @{\\extracolsep{\\fill}}  l @{\\extracolsep{\\fill}}  l  l   l @{\\extracolsep{\\fill}}  l @{\\extracolsep{\\fill}} l}\n')
		output.write('\\hline\\hline\n')
		output.write('\\multicolumn{2}{l}{\\hyperref[2atoms]{2 Atoms}} &\multicolumn{2}{l}{\\hyperref[3atoms]{3 Atoms}}& \\hyperref[4atoms]{4 Atoms} & \\hyperref[5atoms]{5 Atoms} & \\hyperref[6atoms]{6 Atoms} & \\hyperref[7atoms]{7 Atoms} \\\\\n')
		output.write('\\hline\n')
	
		for x in table_list:
		
			output.write(x)			
			
		output.write('\\hline\n\\end{tabular*}\n\\label{two_seven}\n\\end{table*}')
			
def generate_years(write=False,cumulative=False):

	years_list = []

	for x in full_list:
	
		year = getattr(x,'year')
	
		years_list.append(year)
	
	years_dict = {}
	
	for x in years_list:
	
		if x in years_dict:
		
			years_dict[x] += 1
			
		else:
		
			years_dict[x] = 1	
			
	if write == True:
	
		with open('detects_by_year.txt','w') as output:
		
			for x in sorted(years_dict):
			
				output.write('{}\t{}\n' .format(x, years_dict[x]))		
		
	if cumulative == True:
	
		with open('cumulative_by_year.txt','w') as output:
		
			i = 0
		
			for x in sorted(years_dict):
			
				i += years_dict[x]
			
				output.write('{}\t{}\n' .format(x, i))				
	
	return
	
def generate_years_by_atoms():

	final_list = [['Year', 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 'Fullerenes']]
	
	years_list = []
	
	for x in full_list:
	
		year = getattr(x,'year')
	
		years_list.append(year)
	
	years_dict = {}
	
	for x in years_list:
	
		if x in years_dict:
		
			years_dict[x] += 1
			
		else:
		
			years_dict[x] = 1	
	
	for x in sorted(years_dict):
	
		final_entry = [x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		
		for y in full_list:
		
			if y.year == x:
		
				natoms_tmp = y.natoms
			
				if natoms_tmp > 13:
			
					final_entry[13] += 1
				
				else:
			
					final_entry[natoms_tmp-1] += 1	
				
		final_list.append(final_entry)	
		
	with open('years_by_atoms.txt', 'w') as output:		
	
		output.write('Year\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\tfullerenes\n')
		
		i2 = 0
		i3 = 0
		i4 = 0
		i5 = 0
		i6 = 0
		i7 = 0
		i8 = 0
		i9 = 0
		i10 = 0
		i11 = 0
		i12 = 0
		i13 = 0
		ifull = 0
		
		for x in final_list[1:]:
			
			i2 += x[1]
			i3 += x[2]
			i4 += x[3]
			i5 += x[4]
			i6 += x[5]
			i7 += x[6]
			i8 += x[7]
			i9 += x[8]
			i10 += x[9]
			i11 += x[10]
			i12 += x[11]
			i13 += x[12]
			ifull += x[13]
			
			output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n' .format(x[0],i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,ifull))

def count_string(x,y):

	i = 0
	
	for z in full_list:
	
		if x in getattr(z,y):
		
			i += 1
			
	return i	
	
def wavelength_by_mass():

	mass_dict = {}
	
	mass_dict['0-10'] = [0,0,0,0,0,0]
	mass_dict['10-20'] = [0,0,0,0,0,0]
	mass_dict['20-30'] = [0,0,0,0,0,0]
	mass_dict['30-40'] = [0,0,0,0,0,0]
	mass_dict['40-50'] = [0,0,0,0,0,0]
	mass_dict['50-60'] = [0,0,0,0,0,0]
	mass_dict['60-70'] = [0,0,0,0,0,0]
	mass_dict['70-80'] = [0,0,0,0,0,0]
	mass_dict['80-90'] = [0,0,0,0,0,0]
	mass_dict['90-100'] = [0,0,0,0,0,0]
	mass_dict['100-110'] = [0,0,0,0,0,0]
	mass_dict['110-120'] = [0,0,0,0,0,0]
	mass_dict['120-130'] = [0,0,0,0,0,0]
	mass_dict['130+'] = [0,0,0,0,0,0]	
	
	for x in full_list:
	
		if x.mass < 10:
		
			dict_entry = '0-10'
			
		elif x.mass < 20:
		
			dict_entry = '10-20'
			
		elif x.mass < 30:
		
			dict_entry = '20-30'
			
		elif x.mass < 40:
		
			dict_entry = '30-40'
			
		elif x.mass < 50:
		
			dict_entry = '40-50'
			
		elif x.mass < 60:
		
			dict_entry = '50-60'
			
		elif x.mass < 70:
		
			dict_entry = '60-70'
			
		elif x.mass < 80:
		
			dict_entry = '70-80'
			
		elif x.mass < 90:
		
			dict_entry = '80-90'
			
		elif x.mass < 100:
		
			dict_entry = '90-100'
			
		elif x.mass < 110:
		
			dict_entry = '100-110'
			
		elif x.mass < 120:
		
			dict_entry = '110-120'
			
		elif x.mass < 130:
		
			dict_entry = '120-130'
			
		else:
		
			dict_entry = '130+'
	
		waves = x.wavelengths.split(',')
		
		for y in range(len(waves)):
		
			if waves[y].strip() == 'cm':
			
				mass_dict[dict_entry][0] += 1

			if waves[y].strip() == 'mm':
			
				mass_dict[dict_entry][1] += 1	
				
			if waves[y].strip() == 'sub-mm':
			
				mass_dict[dict_entry][2] += 1	
				
			if waves[y].strip() == 'IR':
			
				mass_dict[dict_entry][3] += 1					
				
			if waves[y].strip() == 'Vis':
			
				mass_dict[dict_entry][4] += 1	
				
			if waves[y].strip() == 'UV':
			
				mass_dict[dict_entry][5] += 1	
	
	with open('waves_by_mass.txt','w') as output:
	
		output.write('Mass Range\tcm\tmm\tsub-mm\tIR\tVis\tUV\n')
				
		for x in mass_dict:		
		
			output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n' .format(x,mass_dict[x][0],mass_dict[x][1],mass_dict[x][2],mass_dict[x][3],mass_dict[x][4],mass_dict[x][5]))
			
def mass_by_wavelength():

	waves_dict = {}
	
	waves_dict['cm'] = []
	waves_dict['mm'] = []
	waves_dict['sub-mm'] = []
	waves_dict['IR'] = []
	waves_dict['Vis'] = []
	waves_dict['UV'] = []
	
	for x in full_list:
	
		waves = x.wavelengths.split(',')
		
		for y in range(len(waves)):
		
			if waves[y].strip() == 'cm':
			
				waves_dict['cm'].append(x.mass)

			if waves[y].strip() == 'mm':
			
				waves_dict['mm'].append(x.mass)
				
			if waves[y].strip() == 'sub-mm':
			
				waves_dict['sub-mm'].append(x.mass)
				
			if waves[y].strip() == 'IR':
			
				waves_dict['IR'].append(x.mass)				
				
			if waves[y].strip() == 'Vis':
			
				waves_dict['Vis'].append(x.mass)
				
			if waves[y].strip() == 'UV':
			
				waves_dict['UV'].append(x.mass)
	
	with open('cm_masses.txt','w') as output:
	
		for x in range(len(waves_dict['cm'])):		
		
			output.write('{}\n' .format(waves_dict['cm'][x]))
			
	with open('mm_masses.txt','w') as output:
	
		for x in range(len(waves_dict['mm'])):		
		
			output.write('{}\n' .format(waves_dict['mm'][x]))
			
	with open('sub-mm_masses.txt','w') as output:
	
		for x in range(len(waves_dict['sub-mm'])):		
		
			output.write('{}\n' .format(waves_dict['sub-mm'][x]))
			
	with open('IR_masses.txt','w') as output:
	
		for x in range(len(waves_dict['IR'])):		
		
			output.write('{}\n' .format(waves_dict['IR'][x]))
			
	with open('Vis_masses.txt','w') as output:
	
		for x in range(len(waves_dict['Vis'])):		
		
			output.write('{}\n' .format(waves_dict['Vis'][x]))
			
	with open('UV_masses.txt','w') as output:
	
		for x in range(len(waves_dict['UV'])):		
		
			output.write('{}\n' .format(waves_dict['UV'][x]))
			
def years_by_atoms():

	atoms_dict = {}
	
	atoms_dict['2'] = []
	atoms_dict['3'] = []
	atoms_dict['4'] = []
	atoms_dict['5'] = []
	atoms_dict['6'] = []
	atoms_dict['7'] = []
	atoms_dict['8'] = []
	atoms_dict['9'] = []
	atoms_dict['10'] = []
	atoms_dict['11'] = []
	atoms_dict['12'] = []				
	atoms_dict['13'] = []
	atoms_dict['fullerenes'] = []
				
	for x in full_list:
	
		atoms_str = ''
	
		if x.natoms > 13:
		
			atoms_str = 'fullerenes'
			
		else:
		
			atoms_str = '{}' .format(x.natoms)
			
		atoms_dict[atoms_str].append(x.year)	
	
	with open('years_by_atoms.txt', 'w') as output:
	
		for x in atoms_dict:
		
			output.write('{} Atoms\n' .format(x))
			
			for y in range(len(atoms_dict[x])):
			
				output.write('{}\n' .format(atoms_dict[x][y]))

#make_sources_list runs through the data and returns a full list of unique sources

def make_sources_list():

	sources_list = []
	
	for x in full_list:
	
		for y in x.sources.split(','):
	
			if y.strip() not in sources_list:
		
				sources_list.append(y.strip())
				
	return sorted(sources_list)	
			
#make_du_list generates an ascii file, two column, tab delimited, with the formula and degree of unsaturation for all hydrocarbon molecules (H, C, O, N, S, F, Cl) except fullerenes.
			
def make_du_list():

	with open('du_list.txt','w') as output:

		for x in full_list:
	
			if x.du != None and x.natoms < 60:
		
				output.write('{}\t{}\n' .format(x.formula,x.du))

#count_du_gt runs through the data and counts all degree of unsaturation values greater than 'y'
				
def count_du_gt(y):

	i = 0
	
	for x in full_list:
	
		if x.du != None:
		
			if x.du > y:
			
				i += 1
				
	return i				
	
#make_scopes_list runs through the data and returns a full list of unique telescopes

def make_scopes_list():

	scopes_list = []
	
	for x in full_list:
	
		for y in x.telescopes.split(','):
	
			if y.strip() not in scopes_list:
		
				scopes_list.append(y.strip())
				
	return sorted(scopes_list)

#detects_by_scope_by_year generates an ascii file, tab delimited columns that give the cumulative number of detections for each telescope in each year, for all telescopes with greater than z detections.

def detects_by_scope_by_year(z=4,get_dict=False):

	scopes_dict = {}
	
	scopes_list = make_scopes_list()
	
	for scope in scopes_list:
	
		count = 0
		
		for x in full_list:
		
			if scope in x.telescopes:
			
				count += 1
				
		if count > z:
		
			scopes_dict[scope] = []	
			
	for scope in scopes_dict:
	
		for x in full_list:
		
			if scope in x.telescopes:
			
				scopes_dict[scope].append(x.year)
				
	years = list(range(1963,2019))
	
	
	output_table = []
	
	output_table.append(years)
	
	for scope in scopes_dict:
	
		count = 0
		
		count_list = [scope]
		
		for year in years:
		
			count += scopes_dict[scope].count(year)
			
			count_list.append(count)
			
		output_table.append(count_list)
		
	output_table[0].insert(0,'year')
		
	with open('detects_by_scope_by_year.txt','w') as output:
	
		for x in range(len(output_table[0])):
		
			line = ''
		
			for y in range(len(output_table)):
			
				line += str(output_table[y][x]) + '\t'
				
			output.write(line + '\n')
				
	if get_dict == True:
	
		return scopes_dict		
				
#make_scope_detect_latex makes a latex table for the number of detections per telescope.

def make_scope_detect_latex():

	scopes_list = make_scopes_list()
	
	scopes_dict = {}
	
	for scope in scopes_list:
	
		count = 0
	
		for x in full_list:
		
			if scope in x.telescopes:
			
				count += 1
			
		scopes_dict[scope] = count
		
	with open('detects_scope_table.tex','w') as output:
	
		output.write('\\begin{table}[h!]\n\\centering\n\\footnotesize\n')
		output.write('\\caption{Total number of detections for each facility listed in \\S\\ref{known}.}\n')
		output.write('\\begin{tabular*}{\\columnwidth}{l @{\extracolsep{\\fill}}  c @{\extracolsep{\\fill}}  l @{\extracolsep{\\fill}}  c }\n')
		output.write('\\hline\\hline\n')
		output.write('Facility\t&\t\\#\t&\tFacility\t&\t\\# \\\\\n')
		output.write('\\hline\n')
	
		x = len(scopes_list)
		
		z = int(0.5*x)
		
		for y in range(z):
		
				scope_1 = scopes_list[y]
				scope_2 = scopes_list[y+z]
		
				output.write('{}\t&\t{}\t&\t{}\t&\t{}\\\\\n' .format(scope_1,scopes_dict[scope_1],scope_2,scopes_dict[scope_2]))
				
		output.write('\\hline\n\\end{tabular*}\n\\label{detects_by_scope}\n\\end{table}')								
		
#refs_list(y) is a utility function for getting ascii output of references to feed into summary_list(y)

def refs_list(y):

	refs_list = []

	lab_refs = y.lab_ref.split(';')
	d_refs = y.d_ref.split(';')
	
	if y.notes != None:
		notes = y.notes.strip('*')
	
	refs_list.append('Detection Reference(s)')
	
	for x in range(len(d_refs)):
	
		refs_list.append('[{}] {}' .format(x+1,d_refs[x].strip()))
		
	refs_list.append('\nLaboratory Reference(s)')	
	
	for x in range(len(lab_refs)):
	
		refs_list.append('[{}] {}' .format(x+1,lab_refs[x].strip()))	
		
	if y.notes!= None:
		refs_list.append('\nNotes')
		refs_list.append(notes)
		
	if y.isos != None:
	
		iso_d_refs = y.isos_d_ref.split('[')
		#iso_l_refs = y.isos_l_ref.split('[')
		
		del iso_d_refs[0]
		#del iso_l_refs[0]
	
		refs_list.append('\nIsotopologue Detection Reference(s)')
		
		for x in iso_d_refs:
		
			refs_list.append('[' + x.strip())
		
		#refs_list.append('\nIsotopologue Laboratory Reference(s)')	
			
		#for x in iso_l_refs:
		
		#	refs_list.append('[' + x.strip())	

	if y.ice == True or y.ice == 'Tentative':
	
		refs_list.append('\nIce Reference(s)')
		
		refs_list.append('[Det] {}' .format(y.ice_d_ref))
		refs_list.append('[Lab] {}' .format(y.ice_l_ref))

	if y.ppd == True or y.ppd == 'Tentative':
	
		refs_list.append('\nProtoplanetary Disks Reference(s)')
		
		refs_list.append('[Det] {}' .format(y.ppd_d_ref))
		#refs_list.append('[Lab] {}' .format(y.ppd_l_ref))	
		
	if y.exgal == True or y.exgal == 'Tentative':
	
		refs_list.append('\nExternal Galaxies Reference(s)')
		
		refs_list.append('[Det] {}' .format(y.exgal_d_ref))
		#refs_list.append('[Lab] {}' .format(y.exgal_l_ref))		
		
	if y.exo == True or y.exo == 'Tentative':
	
		refs_list.append('\nExoplanetary Atmospheres Reference(s)')
		
		refs_list.append('[Det] {}' .format(y.exo_d_ref))
		#refs_list.append('[Lab] {}' .format(y.exo_l_ref))		
		
	return refs_list	

#summary_list(y) is a utility function for getting ascii output of summary for printing to a file with output_summary(y)

def summary_list(y):

	summary_list = []

	if isinstance(y,Molecule):

		n_dash = len(y.name) + len(y.formula) +3 

		dashes = '-' * n_dash

		summary_list.append('\n' + dashes)
		summary_list.append('{} ({})' .format(y.name,y.formula))
		summary_list.append(dashes + '\n')
		summary_list.append('Atoms:\t{}' .format(y.natoms))
		summary_list.append('Mass:\t{} amu' .format(y.mass))
		summary_list.append('Year Detected:\t{}' .format(y.year))
		summary_list.append('Source(s):\t{}' .format(y.sources))
		summary_list.append('Telescope(s) Used:\t{}' .format(y.telescopes))

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

		summary_list.append('Attributes:\t{}\n' .format(attr_str))
	
		if y.isos != None:
	
			summary_list.append('Known Isotopologues: {}\n' .format(y.isos))
	
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
	
			summary_list.append('Also Detected In:\t{}' .format(other_str))

		refs_list_tmp = refs_list(y)
		
		for x in refs_list_tmp:
		
			summary_list.append(x)
	
	elif isinstance(y,Source):

		n_dash = len(y.name)
	
		dashes = '-' * n_dash
	
		summary_list.append('\n' + dashes)
		summary_list.append('{}' .format(y.name))
		summary_list.append(dashes + '\n')

		summary_list.append('RA (J2000):\t{}' .format(y.ra))
		summary_list.append('DEC (J2000):\t{}\n' .format(y.dec))
	
		summary_list.append('Generalized Type:\t{}\n' .format(y.type))
	
		summary_list.append('Number of Detections:\t{}\n' .format(y.detects))
	
		summary_list.append('Simbad URL:\t{}' .format(y.simbad_url))
		
		mol_str = 'Molecules Detected in {}' .format(y.name)
	
		dashes = '-' * len(mol_str)
		
		summary_list.append('\n' + dashes)	
		summary_list.append(mol_str)
		summary_list.append(dashes)
	
		summary_list.append('{}' .format(mols_in_source(y.name)).replace("'",'').strip(']').strip('['))
		
	return summary_list	

#calc_max_du(y) calculates the maximum value of the degree of unsaturation a molecule could have

def calc_max_du(y):

	max_du = 1 + 0.5*(y.C*2 + y.N*1)
	
	return max_du
	
#calc_rel_du(y) calculates the relative level of unsaturation for a molecule		

def calc_rel_du(y):

	rel_du = y.du / calc_max_du(y)
	
	return rel_du
	
#make_rel_du_list generates an ascii file, two column, tab delimited, with the formula and relative degree of unsaturation for all hydrocarbon molecules
			
def make_rel_du_list():

	with open('rel_du_list.txt','w') as output:

		for x in full_list:
	
			if x.du != None:
			
				rel_du = calc_rel_du(x)
		
				output.write('{}\t{}\n' .format(x.formula,rel_du))	

#rel_du_by_source_type() returns a tab-delimited ascii file with the relative degrees of unsaturation for each generalized source type, disregarding the fullerenes.
				
def rel_du_by_source_type():

	types_list = make_source_types_list()
	
	types_dict = dict((type,[]) for type in types_list)
	
	for x in full_list:
	
		#get the sources for a molecule
	
		sources = x.sources.split(',')
		
		#initialize a list of source types so we don't overcount some source types
		
		types_tracker = make_source_types_list()
		
		for y in sources:
		
			z = y.strip()
			
			for a in source_tag_list:
			
				if a.name == z:
				
					if x.du != None and x.du < 20:
					
						#check if that type has yet to be used, if it hasn't, go ahead and credit it with a detection
					
						if a.type in types_tracker:
					
							rel_du = calc_rel_du(x)
				
							types_dict[a.type].append(rel_du)
							
							#remove that source type from the list, so it doesn't get credited again for this molecule
							
							types_tracker.remove(a.type)
					
					
	for x in types_dict:
		
		du_values = sorted(types_dict[x])
		
		types_dict[x] = du_values
	
	write_dict = {}
	
	for x in types_dict:
	
		du_str = ''
		
		for y in types_dict[x]:
		
			du_str += '{}\t' .format(y)
			
		du_str.strip('\t')
		
		du_str += '\n'
		
		write_dict[x] = du_str
		
						
					
	with open('rel_du_by_source_type.txt','w') as output:
	
		output.write('SourceType\tDU_vals\n')
	
		for x in sorted(write_dict):
		
			output.write('{}\t{}' .format(x,write_dict[x]))		
			
#make_rel_du_list generates an ascii file, three column, tab delimited, with the formula, mass, and relative degree of unsaturation for all hydrocarbon molecules
			
def make_rel_du_mass_list():

	with open('rel_du_mass_list.txt','w') as output:

		for x in full_list:
	
			if x.du != None:
			
				rel_du = calc_rel_du(x)
		
				output.write('{}\t{}\t{}\n' .format(x.formula,x.mass,rel_du))	

#make_galaxies_list generates a list of unique external galaxies in which detections have been made

def make_galaxies_list():

	raw_list = []
	
	for x in full_list:
	
		if x.exgal != None:
		
			sources_tmp = x.exgal_sources.split(',')
			
			for y in sources_tmp:
			
				raw_list.append(y.strip())
				
	galaxies_list = sorted(list(set(raw_list)))
	
	return galaxies_list

#make_exgal_gals_list generates an ascii file, two column, tab delimited, with the galaxies in which detections have been made, and the number of detections.

def make_exgal_gals_list():

	galaxies_list = make_galaxies_list()
	
	galaxies_dict = dict((galaxy,0) for galaxy in galaxies_list)
	
	for x in full_list:		
	
		if x.exgal != None:
		
			sources_tmp = x.exgal_sources.split(',')
		
			for y in sources_tmp:
		
				source = y.strip()
			
				galaxies_dict[source] += 1
			
	with open('exgal_galaxies_list.txt','w') as output:
	
		for x in galaxies_dict:
		
			output.write('{}\t{}\n' .format(x,galaxies_dict[x]))
			
#atoms_by_wavelength generates an ascii file, tab delimited, of the number of atoms in each molecule detected at each wavelength.

def atoms_by_wavelength():

	waves_dict = {}

	waves_dict['cm'] = []
	waves_dict['mm'] = []
	waves_dict['sub-mm'] = []
	waves_dict['IR'] = []
	waves_dict['Vis'] = []
	waves_dict['UV'] = []
	
	for x in full_list:
	
		for y in x.wavelengths.split(','):
		
			waves_dict[y.strip()].append(x.natoms)
			
	lengths = []
	
	for x in waves_dict:
	
		lengths.append(len(waves_dict[x]))
		
	rows = max(lengths)
	
	with open('atoms_by_wavelength.txt','w') as output:
	
		output.write('cm\tmm\tsub-mm\tIR\tVis\tUV\n')
		
		for x in range(rows):
		
			cm = ''
			mm = ''
			submm = ''
			ir = ''
			vis = ''
			uv = ''
		
			if x < len(waves_dict['cm']):		
				
				cm = waves_dict['cm'][x]
				
			if x < len(waves_dict['mm']):		
				
				mm = waves_dict['mm'][x]
				
			if x < len(waves_dict['sub-mm']):		
				
				submm = waves_dict['sub-mm'][x]
				
			if x < len(waves_dict['IR']):		
				
				ir = waves_dict['IR'][x]
				
			if x < len(waves_dict['Vis']):		
				
				vis = waves_dict['Vis'][x]
				
			if x < len(waves_dict['UV']):		
				
				uv = waves_dict['UV'][x]																				
				
			output.write('{}\t{}\t{}\t{}\t{}\t{}\n' .format(cm,mm,submm,ir,vis,uv))	
				
				
				
				
				
				
				
				
				
				
				
				
			